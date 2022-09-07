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


std::vector<std::string> split(const std::string&, char);
void dumpcell_to_vector_converter(
        double&, double&, double&, double&, double&, double&, double&, double&, double&,
        Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);
void vector_to_dumpcell_converter(
        Eigen::Matrix3d&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&);
void rotationtxt2rotmatrix(std::string&, Eigen::Matrix3d&, int&);
void read_one_timestep_of_dump(
        std::ifstream&, Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::MatrixXd&,
        Eigen::MatrixXd&, std::vector<Eigen::MatrixXd>&, std::vector<int>&, int&, int&);
void write_to_newdump(
        std::ofstream&, int&, int&, Eigen::MatrixXd&, Eigen::MatrixXd&,
        Eigen::MatrixXd&, Eigen::MatrixXd&, std::vector<Eigen::MatrixXd>&, std::vector<int>&);
void count_number_of_rows(std::string&, int&);

