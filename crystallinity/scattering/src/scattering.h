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
#include <kd-tree/kdtree.h>
#include <read_dump/read_dump.h>
//#include <write_dump/write_dump.h>
#include <general/mytools.h>


void rotationtxt2rotmatrix(std::string&, Eigen::Matrix3d&, int&);

