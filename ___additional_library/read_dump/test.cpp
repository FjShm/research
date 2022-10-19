#include "read_dump.h"


int main(){
    ReadDump::ReadDump rd("dump.u.lammpstrj");
    while (rd.read_1frame()){
        std::vector<Eigen::Vector3d> coordinate;
        rd.join_3columns(coordinate, "xu", "yu", "zu");

        std::cout << rd.timestep << std::endl;
        std::cout << rd.num_atoms << std::endl;
        std::cout << coordinate[14920-1] << std::endl;
        std::cout << coordinate[18475-1] << std::endl;
    }
}
