#include "crystallinity-A.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_dump_path = param["output_cry_dump_path"].as<std::string>();
    const int k = param["k"].as<int>();
    const int beta = param["beta"].as<int>();
    const double lath = param["lambda_threshold"].as<double>();
    const double dith = param["neighbor_stem_distance_threshold"].as<double>();
    const double thth = param["theta_deg_threshold"].as<double>();
    const double sqrdDith = dith*dith;
    const double cos_thth = std::cos(thth * M_PI/180.);

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(rot_path, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    std::ifstream dump{dump_path};
    std::ifstream rotxt{rot_path};
    std::ofstream out_dump{out_dump_path, std::ios::out | std::ios::trunc};

    std::string rotxt_row;;

    // skip header of rotation.txt
    for (int i = 0; i < 2; i++) std::getline(rotxt, rotxt_row);


    while(std::getline(rotxt, rotxt_row)){
        // read rotation.txt
        int timestep;
        Eigen::Matrix3d rot;
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);


        // read dump one timestep
        // coordinations, molsは
        // read_one_timestep_of_dump内でresize
        Eigen::Vector3d a, b, c, cell_origin;
        std::vector<Eigen::Vector3d> coordinations;
        std::vector<int> mols;
        int timestep_dump, num_atoms;
        read_one_timestep_of_dump(
                dump, a, b, c, cell_origin, coordinations, mols, timestep_dump, num_atoms);

        int mol_max = *std::max_element(mols.begin(), mols.end());
        int bead_per_chain = num_atoms / mol_max;

        // check timestep
        if (timestep != timestep_dump){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep
                << "dump file:    " << timestep_dump << std::endl;
            return 1;
        }

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

        // 周期境界条件も考えてKB-tree作成
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
            //for (int j = 1; j < idxes.size(); j++){
                int jj = idxes[j] / 27;
                stem_j = stem_vecs[jj];
                if (lambdas[i] >= lath && lambdas[jj] >= lath){
                    double cos =
                        std::abs(stem_i.dot(stem_j)) / (stem_vec_norms[i]*stem_vec_norms[jj]);
                    if (cos >= cos_thth) cry_neighbor_stems[i]++;
                }
            }
        }


        // output new dump
        write_to_newdump(
                out_dump, timestep, num_atoms, a, b, c,
                cell_origin, output_coordinations, mols,
                total_neighbor_stems, cry_neighbor_stems, alphas
                );

        // update progress bar
        ++show_progress;
    }
}

void dumpcell_to_vector_converter(
        double &xlo_b, double &xhi_b, double &xy,
        double &ylo_b, double &yhi_b, double &xz,
        double &zlo_b, double &zhi_b, double &yz,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin
        ){
    double xlo,xhi, ylo,yhi, zlo,zhi;
    xlo = xlo_b - std::min({0., xy, xz, xy + xz});
    xhi = xhi_b - std::max({0., xy, xz, xy + xz});
    ylo = ylo_b - std::min({0., yz});
    yhi = yhi_b - std::max({0., yz});
    zlo = zlo_b;
    zhi = zhi_b;
    cell_origin << xlo, ylo, zlo;
    a << xhi - xlo, 0., 0.;
    b << xy, yhi - ylo, 0.;
    c << xz, yz, zhi - zlo;
}

void read_one_timestep_of_dump(
        std::ifstream &dump,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin,
        std::vector<Eigen::Vector3d> &coordinations,
        std::vector<int> &mols,
        int &timestep_dump,
        int &num_atoms
        ){
    std::string row;

    // TIMESTEP
    std::getline(dump, row);
    std::getline(dump, row);
    timestep_dump = std::stoi(row);

    // NUMBER OF ATOMS
    std::getline(dump, row);
    std::getline(dump, row);
    num_atoms = std::stoi(row);
    coordinations.resize(num_atoms);
    mols.resize(num_atoms);

    // BOX BOUNDS xy xz yz pp pp pp
    double xlo_b,xhi_b,xy, ylo_b,yhi_b,xz, zlo_b,zhi_b,yz;
    std::getline(dump, row);
    dump >> xlo_b >> xhi_b >> xy >> ylo_b >> yhi_b >> xz >> zlo_b >> zhi_b >> yz;
    dumpcell_to_vector_converter(
            xlo_b, xhi_b, xy, ylo_b, yhi_b, xz, zlo_b, zhi_b, yz, a, b, c, cell_origin);
    std::getline(dump, row);

    // ATOMS id mol xu yu zu
    std::getline(dump, row);
    int id, mol;
    double xu, yu, zu;
    for (int i = 0; i < num_atoms; i++){
        dump >> id >> mol >> xu >> yu >> zu;
        mols[id - 1] = mol;
        coordinations[id - 1] << xu, yu, zu;
    }
    std::getline(dump, row);
}

std::vector<std::string> split(const std::string &s, char delim){
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) {
            elems.push_back(item);
        }
    }
    return elems;
}

void rotationtxt2rotmatrix(std::string &row, Eigen::Matrix3d &rot, int &timestep){
    std::vector<std::string> row_split = split(row, ' ');
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

void vector_to_dumpcell_converter(
        Eigen::Matrix3d &dumpcell,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin
        ){
    double xlo,xhi,xy, ylo,yhi,xz, zlo,zhi,yz;
    xlo = cell_origin(0);
    ylo = cell_origin(1);
    zlo = cell_origin(2);
    xy = b(0);
    xz = c(0);
    yz = c(1);
    xhi = a(0) + xlo;
    yhi = b(1) + ylo;
    zhi = c(2) + zlo;

    double xlo_b,xhi_b, ylo_b,yhi_b, zlo_b,zhi_b;
    xlo_b = xlo + std::min({0., xy, xz, xy + xz});
    xhi_b = xhi + std::max({0., xy, xz, xy + xz});
    ylo_b = ylo + std::min({0., yz});
    yhi_b = yhi + std::max({0., yz});
    zlo_b = zlo;
    zhi_b = zhi;

    dumpcell << xlo_b, xhi_b, xy,
             ylo_b, yhi_b, xz,
             zlo_b, zhi_b, yz;
}

void write_to_newdump(
        std::ofstream &out,
        int &timestep,
        int& num_atoms,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin,
        std::vector<Eigen::Vector3d> &coordinations,
        std::vector<int> &mols,
        std::vector<int> &total_neighbor_stems,
        std::vector<int> &cry_neighbor_stems,
        std::vector<int> &alphas
        ){
    out << "ITEM: TIMESTEP\n" << timestep << std::endl;
    out << "ITEM: NUMBER OF ATOMS\n" << num_atoms << std::endl;
    out << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";

    // calculate triclinic box
    Eigen::Matrix3d dumpcell;
    vector_to_dumpcell_converter(dumpcell, a, b, c, cell_origin);
    out << dumpcell(0, 0) << " " << dumpcell(0, 1) << " " << dumpcell(0, 2) << std::endl;
    out << dumpcell(1, 0) << " " << dumpcell(1, 1) << " " << dumpcell(1, 2) << std::endl;
    out << dumpcell(2, 0) << " " << dumpcell(2, 1) << " " << dumpcell(2, 2) << std::endl;

    out << "ITEM: ATOMS id mol xu yu zu cry_stems ";
    out << "total_neighbor_stems ratio(critical) ratio(whole)\n";
    int al = 0;
    for (int i = 0; i < num_atoms; i++){
        out << i + 1 << " " << mols[i] << " " << coordinations[i](0) << " "
            << coordinations[i](1) << " " << coordinations[i](2) << " ";
        if (i == alphas[al] && total_neighbor_stems[al] != 0){
            out << cry_neighbor_stems[al] << " " << total_neighbor_stems[al] << " "
                << (double)cry_neighbor_stems[al]/(double)total_neighbor_stems[al] << " "
                << (double)cry_neighbor_stems[al]/(double)cry_neighbor_stems.size() << std::endl;
            al++;
        } else if (i == alphas[al]){
            out << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
            al++;
        } else {
            out << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
        }
    }
}


void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    for (int i = 0; i < 2; i++) std::getline(in, row);
    while(std::getline(in, row)) max_loop++;
}

