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


std::vector<std::string> split(const std::string&, char);
void rotationtxt2rotmatrix(std::string&, Eigen::Matrix3d&, int&);
void read_one_timestep_of_dump(
        std::ifstream&, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&,
        Eigen::Vector3d&, std::vector<Eigen::Vector3d>&, std::vector<int>&, int&, int&);
void write_to_newdump(
        std::ofstream&, int&, int&, Eigen::Vector3d&, Eigen::Vector3d&,
        Eigen::Vector3d&, Eigen::Vector3d&, std::vector<Eigen::Vector3d>&, std::vector<int>&,
        std::vector<int>&, std::vector<int>&, std::vector<int>&);
void count_number_of_rows(const std::string&, int&);
void sqrdDistance_btw_linesegment_linesegment(
        Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&, double&);

