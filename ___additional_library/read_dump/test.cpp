#include "read_dump.h"


int main(){
    const std::string dir = "./";
    const std::string fname = "dump.u.lammpstrj";
    const std::vector<std::string> headers = {"id", "xu", "yu", "zu"};
    ReadDump::ReadDump rd(dir + fname);

    std::cout << "Do you want to load all of the following file first?\n"
        << dir + fname << "\n\n Y/[n]: ";
    char yn = std::cin.get();
    std::cout << std::endl;
    if (yn == 'Y') rd.read_all_frames();
    while (rd.read_1frame()){
        //rd.header_validation(headers);
        //rd.header_validation("id", "xu", "yu", "zu");
        rd.header_validation("id");
        std::vector<Eigen::Vector3d> coordinate;
        rd.join_3columns(coordinate, "xu", "yu", "zu");

        std::cout << "timestep: " << rd.timestep << std::endl;
        std::cout << "num_atoms: " << rd.num_atoms << std::endl;
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
    }
}
