#include <random>
#include <read_dump/read_dump.h>
#include "../write_dump.h"

std::random_device rd;
std::mt19937 gen(rd());


int main(){
    const std::string dir = "./";
    const std::string fname = "dump.u.lammpstrj";
    const int N(49), M(512);
    ReadDump::ExtraReadDump rd(dir + fname);
    WriteDump::WriteDump wd(dir + "w_" + fname);

    // test read_all_frames, or read_1frame
    std::cout << "Do you want to load all of the following file first?\n"
        << dir + fname << "\n\n Y/[n]: ";
    char yn = std::cin.get();
    std::cout << std::endl;
    if (yn == 'Y') rd.read_all_frames();


    // reading sequence
    while (rd.read_1frame()){

        // test check_if_wanted_frame
        if (!rd.check_if_wanted_frame()){
            std::cout << "\n!! skip to read timestep " << rd.timestep << std::endl;
            continue;
        }

        //wd.timestep = &rd.timestep;
        //wd.num_atoms = &rd.num_atoms;
        //wd.cellbox_a = &rd.cellbox_a;
        //wd.cellbox_b = &rd.cellbox_b;
        //wd.cellbox_c = &rd.cellbox_c;
        //wd.cellbox_origin = &rd.cellbox_origin;
        //wd.atoms_all_data = rd.atoms_all_data;
        //wd.header_map = rd.header_map;

        wd.set_by_ReadDump(rd);
        wd.set_headers("id", "xu", "yu", "zu");

        wd.write_1frame();
    }
}
