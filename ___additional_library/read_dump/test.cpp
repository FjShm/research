#include "read_dump.h"


int main(){
    std::string dir = "./";
    std::string fname = "dump.u.lammpstrj";
    ReadDump::ReadDump rd(dir + fname);
    std::vector<std::string> headers = {"id", "xu", "yu", "zz"};

    std::cout << "Do you want to load all of the following file first?\n"
        << dir + fname << "\n\n Y/[n]: ";
    std::string yn;
    std::cin >> yn;
    std::cout << std::endl;
    if (yn == "Y") rd.read_all_frames();
    while (rd.read_1frame()){
        rd.header_validation(headers);
        std::vector<Eigen::Vector3d> coordinate;
        rd.join_3columns(coordinate, "xu", "yu", "zu");

        std::cout << rd.timestep << std::endl;
        //std::cout << rd.num_atoms << std::endl;
        //std::cout << coordinate[14920-1] << std::endl;
        //std::cout << coordinate[18475-1] << std::endl;
    }
}
