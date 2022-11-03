#include "crystallinity-A.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_dump_path = param["output_cry_dump_path"].as<std::string>();
    const std::string out_cry_path = param["output_cry_text_path"].as<std::string>();
    const int k = param["k"].as<int>();
    const int beta = param["beta"].as<int>();
    const double lath = param["lambda_threshold"].as<double>();
    const double dith = param["neighbor_stem_distance_threshold"].as<double>();
    const double thth = param["theta_deg_threshold"].as<double>();
    const double sqrdDith = dith*dith;
    const double cos_thth = std::cos(thth * M_PI/180.);

    // -------------------------------
    // max loop
    int max_loop = std::count_rows(dump_path, "TIMESTEP");
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
        std::getline(rotxt, rotxt_row)
        // read rotation.txt
        int timestep;
        Eigen::Matrix3d rot;
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);

        // check timestep
        if (timestep != rd.timestep){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep << std::endl
                << "dump file:    " << timestep_dump << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // validation
        rd.header_validation("id", "xu", "yu", "zu");
        rd.add_column_if_not_exist("mol", N, M);

        std::vector<Eigen::Vector3d> coordinations;
        rd.join_3columns(coordinations, "xu", "yu", "zu");
        // calc center of gravity of each polymers
        std::vector<Eigen::Vector3d> cogs(mol_max);
        Eigen::Vector3d sum_tmp;
        sum_tmp << 0., 0., 0.;
        for (int i = 0; i < num_atoms; i++){
            sum_tmp += coordinations[i];
            if (i == num_atoms - 1 || mols[i] < mols[i+1]){
                cogs[mols[i]-1] = sum_tmp / (double)bead_per_chain;
                sum_tmp << 0., 0., 0.;
            }
        }

        // move vector
        std::vector<Eigen::Vector3d> movecs(mol_max);
        for (int i = 0; i < mol_max; i++){
            cogs[i] -= cell_origin;
            int xi, yi, zi;
            zi = std::floor(cogs[i](2) / c(2));
            yi = std::floor((cogs[i](1) - (double)zi*c(1)) / b(1));
            xi = std::floor((cogs[i](0) - (double)yi*b(0) - (double)zi*c(0)) / a(0));
            movecs[i] = (double)xi*a + (double)yi*b + (double)zi*c;
        }

        // move polymer into simulation box
        std::vector<Eigen::Vector3d> position_fixed_coordinations(num_atoms);
        for (int i = 0; i < num_atoms; i++){
            position_fixed_coordinations[i] = coordinations[i] - movecs[mols[i]-1];
        }

        // rotate coordinations
        std::vector<Eigen::Vector3d> output_coordinations(num_atoms);
        for (int i = 0; i < num_atoms; i++){
            output_coordinations[i] = position_fixed_coordinations[i].transpose() * rot;
        }

        // store alpha (bead index)
        int num_stems = mol_max*(bead_per_chain - 2*k - 2*beta);
        int idx = 0;
        std::vector<int> alphas(num_stems);
        for (int i = 0; i < mol_max; i++)
            for (int alpha = k+beta; alpha < bead_per_chain-k-beta; alpha++)
                alphas[idx++] = i*bead_per_chain + alpha;

        // store lamda
        std::vector<double> lambdas(num_stems);
        std::vector<Eigen::Vector3d> stem_vecs(num_stems);
        std::vector<double> stem_vec_norms(num_stems);
        for (int i = 0; i < num_stems; i++){
            Eigen::Vector3d sum_d;
            sum_d << 0., 0., 0.;
            for (int kk = -k; kk <= k; kk++){
                Eigen::Vector3d d;
                d = position_fixed_coordinations[alphas[i]+kk+beta]
                    - position_fixed_coordinations[alphas[i]+kk-beta];
                d /= d.norm();
                sum_d += d;
            }
            double sum_d_norm = sum_d.norm();
            lambdas[i] = sum_d_norm / (2.*(double)k + 1.);
            stem_vecs[i] = sum_d;
            stem_vec_norms[i] = sum_d_norm;
        }

        // 周期境界条件も考えてKD-tree作成
        Eigen::Vector3d zeros = {0., 0., 0.};
        std::vector<Eigen::Vector3d> _a_ = {-a, zeros, a};
        std::vector<Eigen::Vector3d> _b_ = {-b, zeros, b};
        std::vector<Eigen::Vector3d> _c_ = {-c, zeros, c};

        std::vector<Point> points(27*num_stems);
        int counter = 0;
        for (int i = 0; i < num_stems; i++)
            for (int xi = 0; xi < 3; xi++)
                for (int yi = 0; yi < 3; yi++)
                    for (int zi = 0; zi < 3; zi++)
                        points[counter++] = Point(
                                position_fixed_coordinations[alphas[i]]
                                + _a_[xi]
                                + _b_[yi]
                                + _c_[zi]
                                );
        kdt::KDTree<Point> kdtree(points);

        // 隣接stemの総数とcrystal条件を満たすstemの数
        std::vector<int> total_neighbor_stems(num_stems, 0);
        std::vector<int> cry_neighbor_stems(num_stems, 0);
        Point query;
        Eigen::Vector3d stem_i, stem_j;
        for (int i = 0; i < num_stems; i++){
            stem_i = stem_vecs[i];
            query = Point(stem_i);
            std::vector<int> idxes = kdtree.radiusSearch(query, dith);
            total_neighbor_stems[i] = idxes.size();
            for (int j = 0; j < idxes.size(); j++){
                int jj = idxes[j] % num_stems;
                stem_j = stem_vecs[jj];
                if (lambdas[i] >= lath && lambdas[jj] >= lath){
                    double cos =
                        std::abs(stem_i.dot(stem_j)) / (stem_vec_norms[i]*stem_vec_norms[jj]);
                    if (cos >= cos_thth) cry_neighbor_stems[i]++;
                }
            }
        }

        // output cry text
        if (f1rst_loop){
            out_cry << "TimeStep";
            for (size_t i = 0; i < num_stems; i++)
                out_cry << " " << alphas[i];
            out_cry << std::endl;
            f1rst_loop = false;
        }
        out_cry << timestep;
        for (size_t i = 0; i < num_stems; i++){
            if (total_neighbor_stems[i] != 0){
                out_cry << " " << (double)cry_neighbor_stems[i] / (double)total_neighbor_stems[i];
            } else {
                out_cry << " 0";
            }
        }
        out_cry << std::endl;


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

