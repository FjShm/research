#include <iostream>
#include <string>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <algorithm>
#include <map>
#include <boost/progress.hpp>
#include <yaml-cpp/yaml.h>
#include <read_dump/read_dump.h>


std::vector<std::string> split(const std::string&, char);
void rotationtxt2rotmatrix(std::string&, Eigen::Matrix3d&, int&);
void write_to_newdump(
        std::ofstream&, ReadDump::ReadDump&, std::vector<Eigen::Vector3d>&);
void count_number_of_rows(const std::string&, int&);

