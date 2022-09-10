#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <algorithm>
#include <boost/progress.hpp>
#include <yaml-cpp/yaml.h>


std::vector<std::string> split(const std::string&, char);
void read_one_timestep_of_dump(
        std::ifstream&, std::vector<Eigen::Vector3d>&, std::vector<int>&, int&, int&, int&);
void count_number_of_rows(const std::string&, int&);
void linspace(const int& num, std::vector<double>& x_list, double&);
double count_match_value(double&, double&, std::vector<double>&, bool, bool);

