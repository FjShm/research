#include "crystallinity-A.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_dump_path = param["output_cry_dump_path"].as<std::string>();
    const std::string out_cry_path = param["output_cry_text_path"].as<std::string>();
    const std::vector<bool> special_bonds = param["special_bonds"].as< std::vector<bool> >(std::vector<bool>(3, false));
    const int k = param["k"].as<int>();
    const int beta = param["beta"].as<int>();
    const int frames = param["dump_frames"].as<int>(-1);
    const double lath = param["lambda_threshold"].as<double>();
    const double dith = param["neighbor_stem_distance_threshold"].as<double>();
    const double thth = param["theta_deg_threshold"].as<double>();
    int N = param["N"].as<int>(-1);
    int M = param["M"].as<int>(-1);

    const double cos_thth = std::cos(thth * M_PI/180.);

    // -------------------------------
    // max loop
    if (frames == -1){
        std::cout << "It is recommended to put 'dump_frames' in the input file.\n"
            << "This will automatically start counting the number of frames,\n"
            << "but this process puts a heavy load on the I/O of the system.\n";
    }
    int max_loop = frames == -1 ?
        std::count_rows(dump_path, "TIMESTEP"):
        frames;
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    ReadDump::ExtraReadDump rd(dump_path);
    std::ifstream rotxt{rot_path};
    WriteDump::WriteDump wd(out_dump_path);
    std::ofstream out_cry{out_cry_path, std::ios::out | std::ios::trunc};

    std::string rotxt_row;
    bool f1rst_loop = true;

    // skip header of rotation.txt
    for (int i = 0; i < 2; i++) std::getline(rotxt, rotxt_row);

    while(rd.read_1frame()){
        std::getline(rotxt, rotxt_row);
        // read rotation.txt
        int timestep;
        Eigen::Matrix3d rot;
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);

        // check timestep
        if (timestep != rd.timestep){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep << std::endl
                << "dump file:    " << rd.timestep << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // validation
        rd.header_validation("id", "xu", "yu", "zu", "x", "y", "z");
        rd.add_column_if_not_exist("mol", N, M);
        int mol = rd.header_map->at("mol");

        // set unwrapped_coordination, N, M
        std::vector<Eigen::Vector3d> u_coordinations;
        rd.join_3columns(u_coordinations, "xu", "yu", "zu");
        M = M != -1 ? M : rd.max_of_col("mol");
        N = N != -1 ? N : rd.num_atoms/M;

        // calc center of gravity of each polymers
        std::vector<Eigen::Vector3d> cogs(M);
        Eigen::Vector3d sum_tmp = Eigen::Vector3d::Zero();
        for (int i = 0; i < rd.num_atoms; i++){
            sum_tmp += u_coordinations[i];
            if (
                i == rd.num_atoms - 1 ||
                rd.atoms_all_data->coeff(i, mol) < rd.atoms_all_data->coeff(i+1, mol)
            ){
                cogs[(int)rd.atoms_all_data->coeff(i, mol)-1] = sum_tmp / (double)N;
                sum_tmp << Eigen::Vector3d::Zero();
            }
        }

        // move vector
        // 各鎖の重心をセルボックス内へ移動させるために使うベクトル
        std::vector<Eigen::Vector3d> movecs(M);
        for (int i = 0; i < M; i++){
            cogs[i] -= rd.cellbox_origin;
            int xi, yi, zi;
            zi = std::floor(cogs[i](2) / rd.cellbox_c(2));
            yi = std::floor((cogs[i](1) - (double)zi*rd.cellbox_c(1)) / rd.cellbox_b(1));
            xi = std::floor(
                (cogs[i](0) - (double)yi*rd.cellbox_b(0) - (double)zi*rd.cellbox_c(0)) / rd.cellbox_a(0)
            );
            movecs[i] = (double)xi*rd.cellbox_a + (double)yi*rd.cellbox_b + (double)zi*rd.cellbox_c;
        }

        // move polymer into simulation box
        std::vector<Eigen::Vector3d> position_fixed_coordinations(rd.num_atoms);
        for (int i = 0; i < rd.num_atoms; i++){
            position_fixed_coordinations[i]
                = u_coordinations[i] - movecs[(int)rd.atoms_all_data->coeff(i, mol)-1];
        }

        // rotate coordinations
        Eigen::MatrixXd output_coordinations(rd.num_atoms, 3);
        for (int i = 0; i < rd.num_atoms; i++)
            output_coordinations.row(i) << position_fixed_coordinations[i].transpose() * rot;
        rd.append_columns(output_coordinations, true, "xu", "yu", "zu");

        // store lamda
        std::vector<double> lambdas(rd.num_atoms);
        std::vector<Eigen::Vector3d> stem_vecs(rd.num_atoms);
        std::vector<double> stem_vec_norms(rd.num_atoms);
        std::vector<bool> isStem(rd.num_atoms, false);
        int lb = k + beta; // 1本の鎖内でstemを定義出来るビーズindexの下限
        int ub = N - k - beta;
        for (int i = 0; i < rd.num_atoms; i++){
            int n = i % N;
            // stemを定義出来ない範囲
            if (n < lb || ub <= n){
                lambdas[i] = 0.;
                stem_vecs[i].setZero();
                stem_vec_norms[i] = 0.;
            // stem定義可能
            } else {
                Eigen::Vector3d sum_d = Eigen::Vector3d::Zero();
                isStem[i] = true;
                for (int kk = -k; kk <= k; kk++){
                    Eigen::Vector3d d;
                    d = position_fixed_coordinations[i+kk+beta]
                        - position_fixed_coordinations[i+kk-beta];
                    sum_d += d/d.norm();
                }
                double sum_d_norm = sum_d.norm();
                lambdas[i] = sum_d_norm / (2.*(double)k + 1.);
                stem_vecs[i] = sum_d;
                stem_vec_norms[i] = sum_d_norm;
            }
        }

        // 周期境界条件も考えてKD-tree作成
        std::vector<Eigen::Vector3d> wrapped_coordinations;
        rd.join_3columns(wrapped_coordinations, "x", "y", "z");
        const Eigen::Vector3d zeros = Eigen::Vector3d::Zero();
        const std::vector<Eigen::Vector3d> _a_ = {-rd.cellbox_a, zeros, rd.cellbox_a};
        const std::vector<Eigen::Vector3d> _b_ = {-rd.cellbox_b, zeros, rd.cellbox_b};
        const std::vector<Eigen::Vector3d> _c_ = {-rd.cellbox_c, zeros, rd.cellbox_c};

        std::vector<Point> points(27*M*(ub-lb));
        std::vector<int> kdtree_idx_to_molcule_id(27*M*(ub-lb));
        int counter = 0;
        for (int i = 0; i < rd.num_atoms; i++){
            if (!isStem[i]) continue;
            for (int xi = 0; xi < 3; xi++)
                for (int yi = 0; yi < 3; yi++)
                    for (int zi = 0; zi < 3; zi++){
                        kdtree_idx_to_molcule_id[counter] = i;
                        points[counter++] = Point(
                            wrapped_coordinations[i]
                            + _a_[xi] + _b_[yi] + _c_[zi]
                        );
                    }
        }
        kdt::KDTree<Point> kdtree(points);

        // 隣接stemの総数とcrystal条件を満たすstemの数
        Eigen::VectorXd total_neighbor_stems = Eigen::VectorXd::Zero(rd.num_atoms);
        Eigen::VectorXd cry_neighbor_stems = Eigen::VectorXd::Zero(rd.num_atoms);
        Point query;
        Eigen::Vector3d stem_i, stem_j;
        for (int i = 0; i < rd.num_atoms; i++){
            if (!isStem[i]) continue;
            stem_i = stem_vecs[i];
            
            // stem_i の中央のみについてquery作成
            query = Point(wrapped_coordinations[i]);
            
            // idx == i, isStem[idx] == false, special_bonds == false は除外する
            std::vector<int> idxes;
            for (int idx : kdtree.radiusSearch(query, dith)){
                int molecule_idx = kdtree_idx_to_molcule_id[idx]; // 0 <= id < N
                int mol_i = rd.atoms_all_data->coeff(i, mol);
                int mol_j = rd.atoms_all_data->coeff(molecule_idx, mol);
                if (
                    molecule_idx != i &&
                    isStem[molecule_idx] &&
                    is_special_bond(special_bonds, i, molecule_idx, mol_i, mol_j)
                ){
                    idxes.emplace_back(molecule_idx);
                }
            }
            total_neighbor_stems(i) = idxes.size();
            for (auto j : idxes){
                stem_j = stem_vecs[j];
                if (lambdas[i] >= lath && lambdas[j] >= lath){
                    double cos =
                        std::abs(stem_i.dot(stem_j)) / (stem_vec_norms[i]*stem_vec_norms[j]);
                    if (cos >= cos_thth) cry_neighbor_stems(i) += 1;
                }
            }
        }

        // calculate ratio
        Eigen::VectorXd ratio(rd.num_atoms);
        for (size_t i = 0; i < rd.num_atoms; i++)
            ratio(i) =
                (int)total_neighbor_stems(i) == 0 ? 0 : cry_neighbor_stems(i)/total_neighbor_stems(i);

        // output cry text
        if (f1rst_loop){
            out_cry << "TimeStep";
            for (size_t i = 0; i < rd.num_atoms; i++){
                if (!isStem[i]) continue;
                out_cry << " " << i;
            }
            out_cry << std::endl;
            f1rst_loop = false;
        }
        out_cry << rd.timestep;
        for (size_t i = 0; i < rd.num_atoms; i++){
            if (!isStem[i]) continue;
            out_cry << " " << ratio(i);
        }
        out_cry << std::endl;

        // output dump
        rd.append_column(cry_neighbor_stems, "cry_neighbor_stems");
        rd.append_column(total_neighbor_stems, "total_neighbor_stems");
        rd.append_column(ratio, "cry_neighbor_stems(ratio)");
        wd.set_by_ReadDump(rd);
        wd.set_headers("id", "mol", "xu", "yu", "zu",
                "total_neighbor_stems", "cry_neighbor_stems", "cry_neighbor_stems(ratio)");
        wd.write_1frame();

        // update progress bar
        ++show_progress;
    }
}


void rotationtxt2rotmatrix(std::string &row, Eigen::Matrix3d &rot, int &timestep){
    std::vector<std::string> row_split = std::split(row, ' ');
    timestep = std::stoi(row_split[0]);
    rot(0, 0) = std::stod(row_split[1]);
    rot(0, 1) = std::stod(row_split[2]);
    rot(0, 2) = std::stod(row_split[3]);
    rot(1, 0) = std::stod(row_split[4]);
    rot(1, 1) = std::stod(row_split[5]);
    rot(1, 2) = std::stod(row_split[6]);
    rot(2, 0) = std::stod(row_split[7]);
    rot(2, 1) = std::stod(row_split[8]);
    rot(2, 2) = std::stod(row_split[9]);
}

// special_bond かどうか判定する
// special_bond ... 同一鎖内, nつ隣のビーズであり, 結晶化度の計算に含めるもの
bool is_special_bond(
    const std::vector<bool> &special_bonds,
    const int &i,
    const int &j,
    const int &mol_i,
    const int &mol_j
){
    if (mol_i != mol_j) return true;
    int diff = std::abs(i - j);
    if (diff > special_bonds.size()) return true;
    return special_bonds[diff - 1];
}

