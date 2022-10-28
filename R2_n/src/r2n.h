#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <boost/progress.hpp>
#include <yaml-cpp/yaml.h>
#include <read_dump/read_dump.h>

void count_number_of_rows(const std::string&, int&);
void compute_R2_n(ReadDump::ReadDump&, int, int, Eigen::VectorXd&);
