#include <random>
#include "../read_dump.h"

std::random_device rd;
std::mt19937 gen(rd());


int main(){
    const std::string dir = "./";
    const std::string fname = "_dump.u.lammpstrj";
    const int N(49), M(512);
    ReadDump::ExtraReadDump rd(dir + fname);

    // test read_all_frames, or read_1frame
    std::cout << "Do you want to load all of the following file first?\n"
        << dir + fname << "\n\n Y/[n]: ";
    char yn = std::cin.get();
    std::cout << std::endl;
    if (yn == 'Y') rd.read_all_frames();

    // test set_want_frames
    int mode;
    std::cout << "Select read_dump mode\n"
        << " 0: all\n 1: timestep\n 2: ratio\n: ";
    std::cin >> mode;
    if (mode == 1){
        std::vector<int> t(3);
        std::cout << "Type 3 timesteps below...\n";
        std::cin >> t[0] >> t[1] >> t[2];
        rd.set_want_frames(t);
    } else if (mode == 2){
        std::vector<double> r(3);
        std::cout << "Type 3 ratios (0 <= ratio <= 1) below...\n";
        std::cin >> r[0] >> r[1] >> r[2];
        rd.set_want_frames(r);
    } else if (mode != 0){
        std::cout << "Invalid read_dump mode\n";
        std::exit(EXIT_FAILURE);
    }

    // reading sequence
    bool f1rst_loop = true;
    while (rd.read_1frame()){
        if (f1rst_loop){
            // test header_validation
            std::string col1, col2, col3;
            std::cout << "Type 3 header names of dump ATOMS:\n";
            std::cin >> col1 >> col2 >> col3;
            rd.header_validation(col1, col2, col3);
            f1rst_loop = false;
        }

        // test check_if_wanted_frame
        if (!rd.check_if_wanted_frame()){
            std::cout << "\n!! skip to read timestep " << rd.timestep << std::endl;
            continue;
        }

        // test add_mol
        std::uniform_int_distribution<> dist(0, N*M-1);
        int id = dist(gen);
        rd.add_column_if_not_exist("mol", N, M);
        int mol = rd.header_map->at("mol");
        int mol_origin = rd.header_map->at("mol_origin");
        
        std::cout << "id: " << id << std::endl;
        std::cout << "added mol: " << rd.atoms_all_data->coeff(id, mol) << std::endl;
        std::cout << "original mol: " << rd.atoms_all_data->coeff(id, mol_origin) << std::endl;

        // test reference to member variables
        std::cout << "timestep: " << rd.timestep << std::endl;
        std::cout << "num_atoms: " << rd.num_atoms << std::endl;
        std::cout << "num_frames: " << rd.num_frames << std::endl;
        std::cout << "header_map:\n";
        for (const auto& [key, val] : *(rd.header_map)){
            std::cout << "  " << key << ": " << val << std::endl;
        }
        std::cout << "cellbox_a: \n" << rd.cellbox_a << std::endl;
        std::cout << "cellbox_b: \n" << rd.cellbox_b << std::endl;
        std::cout << "cellbox_c: \n" << rd.cellbox_c << std::endl;
        std::cout << "cellbox_origin: \n" << rd.cellbox_origin << std::endl;

        // test join_3columns, and direct reference to atoms_all_data
        std::vector<Eigen::Vector3d> coordinate;
        rd.join_3columns(coordinate, "xu", "yu", "zu");
        std::cout << "joinした座標:\n";
        std::cout << coordinate[14920-1] << std::endl;
        std::cout << "直接参照した座標:\n";
        int xu = rd.header_map->at("xu");
        int yu = rd.header_map->at("yu");
        int zu = rd.header_map->at("zu");
        std::cout << rd.atoms_all_data->coeff(14920-1, xu) << std::endl;
        std::cout << rd.atoms_all_data->coeff(14920-1, yu) << std::endl;
        std::cout << rd.atoms_all_data->coeff(14920-1, zu) << std::endl;
        std::cout << std::endl;

        // test max_of_col
        int maxmol = rd.max_of_col("mol");
        int minmol = rd.min_of_col("mol");
        double maxxu_d = rd.max_of_col("xu");
        int maxxu_i = rd.max_of_col("xu");
        std::cout << "max mol: " << maxmol << std::endl;
        std::cout << "min mol: " << minmol << std::endl;
        std::cout << "max xu (double): " << maxxu_d << std::endl;
        std::cout << "max xu (int): " << maxxu_i << std::endl;

        // test append_column
        std::cout << ">>> test append_column >>>\n";
        Eigen::MatrixXd apcol(rd.num_atoms, 1);
        for (size_t i = 0; i < apcol.rows(); i++)
            apcol(i) = i;
        rd.append_column(apcol, "apcol");
        int apcolidx = rd.header_map->at("apcol");
        for (const auto& [key, val] : *(rd.header_map))
            std::cout << "  " << key << ": " << val << std::endl;
        std::cout << rd.atoms_all_data->coeff(id, apcolidx) << std::endl;
        // replace=false
        for (size_t i = 0; i < apcol.rows(); i++)
            apcol(i) = i*2.;
        rd.append_column(apcol, "apcol", false);
        for (const auto& [key, val] : *(rd.header_map))
            std::cout << "  " << key << ": " << val << std::endl;
        std::cout << rd.atoms_all_data->coeff(id, apcolidx) << std::endl;
        // replace=true
        double tmp = rd.atoms_all_data->coeff(id, apcolidx);
        rd.append_column(apcol, "apcol", true);
        for (const auto& [key, val] : *(rd.header_map))
            std::cout << "  " << key << ": " << val << std::endl;
        std::cout << rd.atoms_all_data->coeff(id, apcolidx) << std::endl;
        std::cout << "(check if " << rd.atoms_all_data->coeff(id, apcolidx)
            << " = 2 * " << tmp  << ")" << std::endl;
        std::cout << " <<< test append_column over <<<\n";

        // test append_columns
        double abef = rd.atoms_all_data->coeff(id, apcolidx);
        std::cout << abef << "*" << 3 << "=" << abef*3 << std::endl;
        Eigen::MatrixXd apcols(apcol.rows(), 3);
        apcols << apcol*3, apcol, apcol;
        rd.append_columns(apcols, true, "apcol", "apcol2", "apcol3");
        std::cout << rd.atoms_all_data->coeff(id, apcolidx) << std::endl;
        int a2 = rd.header_map->at("apcol2");
        int a3 = rd.header_map->at("apcol3");
        std::cout << rd.atoms_all_data->coeff(id, a2) << std::endl;
        std::cout << rd.atoms_all_data->coeff(id, a3) << std::endl;
    }
    std::cout << "num_frames: " << rd.num_frames << std::endl;
}
