#include "../mytools.h"


int main(){
    // test split
    std::string line_s = "test row  wa.";
    std::string line_d = "1 0  -5.4        ";
    std::vector<std::string> splited_s = std::split(line_s, ' ');
    std::vector<std::string> splited_d = std::split(line_d, ' ');
    for (std::string s : splited_s)
        std::cout << s << std::endl;
    for (std::string d : splited_d)
        std::cout << std::stod(d) << std::endl;

    // test linspace
    std::vector<double> lin = std::linspace(-10, 10, 51);
    for (double d : lin)
        std::cout << d << std::endl;

    // test count_row
    std::string fname = "dump.u.lammpstrj";
    int counter = std::count_rows(fname);
    std::cout << "# of all rows: " << counter << std::endl;
    counter = std::count_rows(fname, "ITEM: TIMESTEP");
    std::cout << "# of 'ITEM: TIMESTEP': " << counter << std::endl;
    counter = std::count_rows(fname, "ITEM:");
    std::cout << "# of 'ITEM:': " << counter << std::endl;
}
