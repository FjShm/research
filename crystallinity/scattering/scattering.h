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
#include <boost/progress.hpp>
#include <yaml-cpp/yaml.h>
#include <read_dump/read_dump.h>

#include <omp.h>


std::vector<std::string> split(const std::string&, char);
void count_number_of_rows(const std::string&, int&);
std::vector<double> linspace(const double&, const double&, int);
void rotationtxt2rotmatrix(std::string&, Eigen::Matrix3d&, int&);

