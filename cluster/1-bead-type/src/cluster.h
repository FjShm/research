#include <iostream>
#include <string>
//#include <cmath>
//#include <cfloat>
#include <fstream>
//#include <sstream>
#include <vector>
#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/LU>
//#include <algorithm>
#include <boost/progress.hpp>
#include <yaml-cpp/yaml.h>
#include <kd-tree/kdtree.h>
#include <read_dump/read_dump.h>



class Point : public std::array<double, 3>
{
public:
	static const int DIM = 3;

	Point(){}
	Point(Eigen::Vector3d _xyz){ 
		for (int i = 0; i < 3; i++)
            (*this)[i] = _xyz(i);
	}
};


int measure_cluster_size(
    int,
    std::vector<bool>&,
    kdt::KDTree<Point>&,
    const double&,
    std::vector<Eigen::Vector3d>&,
    const int&,
    int cluster_size = 0
);
