#include "rotate_dump.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_dump_path = param["output_rotated_dump_path"].as<std::string>();
    const bool ifMove = param["move_to_cell"].as<bool>();
    int N = param["N"].as<int>(-1);
    int M = param["M"].as<int>(-1);

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(rot_path, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    ReadDump::ReadDump rd(dump_path);
    std::ifstream rotxt{rot_path};
    std::ofstream out_dump{out_dump_path, std::ios::out | std::ios::trunc};

    std::string rotxt_row;;

    // skip header of rotation.txt
    for (int i = 0; i < 2; i++) std::getline(rotxt, rotxt_row);


    while(rd.read_1frame()){
        // read rotation.txt
        int timestep;
        Eigen::Matrix3d rot;
        std::getline(rotxt, rotxt_row);
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);


        // read dump one timestep
        std::vector<Eigen::Vector3d> coordinations;
        rd.join_3columns(coordinations, "xu", "yu", "zu");
        int mol, id;
        if (rd.header_map.count("mol") == 0){
            id = rd.header_map.at("id");
            mol = rd.atoms_all_data.cols();
            rd.header_map.insert(std::make_pair("mol", mol));
            rd.atoms_all_data.conservativeResize(rd.num_atoms, mol+1);
            if (N != -1){
                for (int i = 0; i < rd.num_atoms; i++)
                    rd.atoms_all_data(i, mol) = ((int)rd.atoms_all_data(i, id) - 1) / N + 1;
            } else if (M != -1){
                N = rd.num_atoms / M;
                for (int i = 0; i < rd.num_atoms; i++)
                    rd.atoms_all_data(i, mol) = ((int)rd.atoms_all_data(i, id) - 1) / N + 1;
            } else {
                std::cout << "Since there is no 'mol' in ATOMS in the dump file,"
                    << " it is necessary to write N or M in the "
                    << argv[1] << ".\n";
                std::exit(EXIT_FAILURE);
            }
        } else {
            mol = rd.header_map.at("mol");
        }

        int mol_max = rd.atoms_all_data.col(mol).array().maxCoeff();
        int bead_per_chain = rd.num_atoms / mol_max;

        // check timestep
        if (timestep != rd.timestep){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep << std::endl
                << "dump file:    " << rd.timestep << std::endl;
            return 1;
        }

        if (ifMove){
            // calc center of gravity of each polymers
            std::vector<Eigen::Vector3d> cogs(mol_max);
            Eigen::Vector3d sum_tmp;
            sum_tmp << 0., 0., 0.;
            for (int i = 0; i < rd.num_atoms; i++){
                sum_tmp += coordinations[i];
                if (i == rd.num_atoms - 1 || rd.atoms_all_data(i, mol) < rd.atoms_all_data(i+1, mol)){
                    cogs[rd.atoms_all_data(i, mol)-1] = sum_tmp / (double)bead_per_chain;
                    sum_tmp << 0., 0., 0.;
                }
            }

            // move vector
            std::vector<Eigen::Vector3d> movecs(mol_max);
            for (int i = 0; i < mol_max; i++){
                cogs[i] -= rd.cellbox_origin;
                int xi, yi, zi;
                zi = std::floor(cogs[i](2) / rd.cellbox_c(2));
                yi = std::floor((cogs[i](1) - (double)zi*rd.cellbox_c(1)) / rd.cellbox_b(1));
                xi = std::floor(
                    (cogs[i](0) - (double)yi*rd.cellbox_b(0) - (double)zi*rd.cellbox_c(0))
                    / rd.cellbox_a(0)
                );
                movecs[i] =
                    (double)xi*rd.cellbox_a + (double)yi*rd.cellbox_b + (double)zi*rd.cellbox_c;
            }

            // move polymer into simulation box
            for (int i = 0; i < rd.num_atoms; i++){
                coordinations[i] -= movecs[rd.atoms_all_data(i, mol)-1];
            }
        }

        // rotate coordinations
        for (int i = 0; i < rd.num_atoms; i++){
            coordinations[i] = coordinations[i].transpose() * rot;
        }


        // output new dump
        write_to_newdump(out_dump, rd, coordinations);

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
        ReadDump::ReadDump &rd,
        std::vector<Eigen::Vector3d> &coordinations
        ){
    out << "ITEM: TIMESTEP\n" << rd.timestep << std::endl;
    out << "ITEM: NUMBER OF ATOMS\n" << rd.num_atoms << std::endl;
    out << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";

    // calculate triclinic box
    Eigen::Matrix3d dumpcell;
    vector_to_dumpcell_converter(
        dumpcell, rd.cellbox_a, rd.cellbox_b, rd.cellbox_c, rd.cellbox_origin);
    out << dumpcell(0, 0) << " " << dumpcell(0, 1) << " " << dumpcell(0, 2) << std::endl;
    out << dumpcell(1, 0) << " " << dumpcell(1, 1) << " " << dumpcell(1, 2) << std::endl;
    out << dumpcell(2, 0) << " " << dumpcell(2, 1) << " " << dumpcell(2, 2) << std::endl;

    out << "ITEM: ATOMS id mol xu yu zu\n";
    int mol = rd.header_map.at("mol");
    for (int i = 0; i < rd.num_atoms; i++){
        out << i + 1 << " " << rd.atoms_all_data(i, mol) << " " << coordinations[i](0) << " "
            << coordinations[i](1) << " " << coordinations[i](2) << std::endl;
    }
}


void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    for (int i = 0; i < 2; i++) std::getline(in, row);
    while(std::getline(in, row)) max_loop++;
}
