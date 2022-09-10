#include "rotate_dump.h"


int main(int argc, char *argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_dump_path = param["output_rotated_dump_path"].as<std::string>();

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


    bool f1rst_loop = true;
    Eigen::MatrixXd cell_origin0, a0, b0, c0;
    while(std::getline(rotxt, rotxt_row)){
        // read rotation.txt
        int timestep;
        Eigen::Matrix3d rot;
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);

        // read dump one timestep
        Eigen::MatrixXd a(1, 3), b(1, 3), c(1, 3), cell_origin(1, 3);
        std::vector<Eigen::MatrixXd> coordinations(1, Eigen::MatrixXd (1, 3));
        std::vector<int> mols;
        int timestep_dump, num_atoms;
        read_one_timestep_of_dump(
                dump, a, b, c, cell_origin, coordinations, mols, timestep_dump, num_atoms);
        if (f1rst_loop){
            cell_origin0 = cell_origin;
            a0 = a;
            b0 = b;
            c0 = c;
            f1rst_loop = false;
        }

        // check timestep
        if (timestep != timestep_dump){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep
                << "dump file:    " << timestep_dump << std::endl;
            return 1;
        }

        // rotate simulation cell
        Eigen::MatrixXd new_a(1, 3), new_b(1, 3), new_c(1, 3), new_cell_origin(1, 3);
        new_cell_origin = cell_origin * rot;
        new_a = a * rot;
        new_b = b * rot;
        new_c = c * rot;

        // fix position of cell
        new_cell_origin = cell_origin0;

        // rotate coordinations
        std::vector<Eigen::MatrixXd> new_coordinations(num_atoms, Eigen::MatrixXd (1, 3));
        for (int i = 0; i < num_atoms; i++){
            new_coordinations[i] = coordinations[i] * rot;
        }

        // move the center of gravity of the polymer's coordinates into the cell
        // converte basis vector
        Eigen::MatrixXd ex(1, 3), ey(1, 3), ez(1, 3);
        ex << 1., 0., 0.;
        ey << 0., 1., 0.;
        ez << 0., 0., 1.;
        Eigen::Matrix3d M;
        // M << new_a(0, 0), new_b(0, 0), new_c(0, 0),
        //   new_a(0, 1), new_b(0, 1), new_c(0, 1),
        //   new_a(0, 2), new_b(0, 2), new_c(0, 2);
        M << a0(0, 0), b0(0, 0), c0(0, 0),
          a0(0, 1), b0(0, 1), c0(0, 1),
          a0(0, 2), b0(0, 2), c0(0, 2);
        Eigen::Matrix3d M_inv = M.inverse();

        // calculate center of gravity of each polymer
        Eigen::MatrixXd cg(1, 3);
        cg << 0., 0., 0.;
        int num_beads_per_chain = 0;
        for (int i = 0; i < num_atoms; i++){
            num_beads_per_chain++;
            cg += new_coordinations[i];
            if ((i == num_atoms - 1) || (mols[i + 1] > mols[i])){
                cg /= num_beads_per_chain;
                // cg -= new_cell_origin;
                cg -= cell_origin0;
                double p, q, r;
                p = cg(0, 0); q = cg(0, 1); r = cg(0, 2);
                double coeff_a = p*M_inv(0, 0) + q*M_inv(0, 1) + r*M_inv(0, 2);
                double coeff_b = p*M_inv(1, 0) + q*M_inv(1, 1) + r*M_inv(1, 2);
                double coeff_c = p*M_inv(2, 0) + q*M_inv(2, 1) + r*M_inv(2, 2);
                if (coeff_a < 0) coeff_a -= 1.;
                if (coeff_b < 0) coeff_b -= 1.;
                if (coeff_c < 0) coeff_c -= 1.;
                for (int j = 0; j < num_beads_per_chain; j++){
                    // new_coordinations[i - j] -= (double)((int)coeff_a) * new_a;
                    // new_coordinations[i - j] -= (double)((int)coeff_b) * new_b;
                    // new_coordinations[i - j] -= (double)((int)coeff_c) * new_c;
                    new_coordinations[i - j] -= (double)((int)coeff_a) * a0;
                    new_coordinations[i - j] -= (double)((int)coeff_b) * b0;
                    new_coordinations[i - j] -= (double)((int)coeff_c) * c0;
                }
                num_beads_per_chain = 0;
                cg << 0., 0., 0.;
            }
        }

        // output new dump
        write_to_newdump(
                out_dump, timestep, num_atoms, new_a, new_b, new_c,
                new_cell_origin, new_coordinations, mols
                );

        // update progress bar
        ++show_progress;
    }
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


void dumpcell_to_vector_converter(
        double &xlo_b, double &xhi_b, double &xy,
        double &ylo_b, double &yhi_b, double &xz,
        double &zlo_b, double &zhi_b, double &yz,
        Eigen::MatrixXd &a,
        Eigen::MatrixXd &b,
        Eigen::MatrixXd &c,
        Eigen::MatrixXd &cell_origin
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

void vector_to_dumpcell_converter(
        Eigen::Matrix3d &dumpcell,
        Eigen::MatrixXd &a,
        Eigen::MatrixXd &b,
        Eigen::MatrixXd &c,
        Eigen::MatrixXd &cell_origin
    ){
    double xlo,xhi,xy, ylo,yhi,xz, zlo,zhi,yz;
    xlo = cell_origin(0, 0);
    ylo = cell_origin(0, 1);
    zlo = cell_origin(0, 2);
    xy = b(0, 0);
    xz = c(0, 0);
    yz = c(0, 1);
    xhi = a(0, 0) + xlo;
    yhi = b(0, 1) + ylo;
    zhi = c(0, 2) + zlo;

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

void read_one_timestep_of_dump(
        std::ifstream &dump,
        Eigen::MatrixXd &a,
        Eigen::MatrixXd &b,
        Eigen::MatrixXd &c,
        Eigen::MatrixXd &cell_origin,
        std::vector<Eigen::MatrixXd> &coordinations,
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
    coordinations.resize(num_atoms, Eigen::MatrixXd (1, 3));
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
 
void write_to_newdump(
        std::ofstream &out,
        int &timestep,
        int& num_atoms,
        Eigen::MatrixXd &a,
        Eigen::MatrixXd &b,
        Eigen::MatrixXd &c,
        Eigen::MatrixXd &cell_origin,
        std::vector<Eigen::MatrixXd> &coordinations,
        std::vector<int> &mols
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

    out << "ITEM: ATOMS id mol xu yu zu\n";
    for (int i = 0; i < num_atoms; i++){
        out << i + 1 << " " << mols[i] << " " << coordinations[i](0, 0) << " "
            << coordinations[i](0, 1) << " " << coordinations[i](0, 2) << std::endl;
    }
}
 
void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    for (int i = 0; i < 2; i++) std::getline(in, row);
    while(std::getline(in, row)) max_loop++;
}
