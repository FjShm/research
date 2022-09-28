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

void count_number_of_rows(const std::string&, int&);
void compute_R2_n(std::ifstream&, int, int, int, Eigen::VectorXd&, std::string);
