#include "legendre_p2.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string out_cry_path = param["output_cry_text_path"].as<std::string>();
    const int k = param["k"].as<int>();
    const double dith = param["neighbor_stem_distance_threshold"].as<double>();

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(dump_path, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    std::ifstream dump{dump_path};
    std::ofstream out_cry{out_cry_path, std::ios::out | std::ios::trunc};

    std::string rotxt_row;
    bool f1rst_loop = true;


    int looper = 0;
    while(looper < max_loop){
        // read dump one timestep
        // coordinations, molsは
        // read_one_timestep_of_dump内でresize
        Eigen::Vector3d a, b, c, cell_origin;
        std::vector<Eigen::Vector3d> coordinations;
        std::vector<int> mols;
        int timestep, num_atoms;
        read_one_timestep_of_dump(
                dump, a, b, c, cell_origin, coordinations, mols, timestep, num_atoms);

        int mol_max = *std::max_element(mols.begin(), mols.end());
        int bead_per_chain = num_atoms / mol_max;

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

        // 周期境界条件も考えてKD-tree作成
        Eigen::Vector3d zeros = {0., 0., 0.};
        std::vector<Eigen::Vector3d> _a_ = {-a, zeros, a};
        std::vector<Eigen::Vector3d> _b_ = {-b, zeros, b};
        std::vector<Eigen::Vector3d> _c_ = {-c, zeros, c};

        std::vector<Point> points(27*position_fixed_coordinations.size());
        int counter = 0;
        for (int i = 0; i < position_fixed_coordinations.size(); i++)
            for (int xi = 0; xi < 3; xi++)
                for (int yi = 0; yi < 3; yi++)
                    for (int zi = 0; zi < 3; zi++)
                        points[counter++] = Point(
                                position_fixed_coordinations[i]
                                + _a_[xi]
                                + _b_[yi]
                                + _c_[zi]
                                );
        kdt::KDTree<Point> kdtree(points);

        // 近傍bond(i,j)について考える
        // i<jについてのみ計算
        double cos2_ave = 0;
        counter = 0;
        Point query;
        Eigen::Vector3d stem_i, stem_j;
        for (size_t i = 0; i < position_fixed_coordinations.size()-k; i++){
            if (mols[i] < mols[i+k]) continue;
            Eigen::Vector3d bond_i =
                position_fixed_coordinations[i+k] - position_fixed_coordinations[i];
            query = Point(bond_i);
            std::vector<int> idxes = kdtree.radiusSearch(query, dith);
            for (size_t j = 0; j < idxes.size(); j++){
                int idxesj_ = idxes[j] % position_fixed_coordinations.size();
                if (
                        (idxes[j] < position_fixed_coordinations.size()) &&
                        (i < idxes[j]) &&
                        (mols[idxesj_] < mols[idxesj_+k])
                    ) continue;
                Eigen::Vector3d bond_j =
                    position_fixed_coordinations[idxesj_+k] - position_fixed_coordinations[idxesj_];
                double dot = bond_i.dot(bond_j);
                cos2_ave += dot*dot  / (bond_i.squaredNorm() * bond_j.squaredNorm());
                counter++;
            }
        }
        cos2_ave /= counter;
        double P2 = 0.5 * (3*cos2_ave - 1.);

        // output cry text
        out_cry << timestep << "," << P2 << std::endl;


        // update progress bar
        ++show_progress;
        ++looper;
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


void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    while(std::getline(in, row)){
        if ("ITEM: TIMESTEP" == row)
            max_loop++;
    }
}

