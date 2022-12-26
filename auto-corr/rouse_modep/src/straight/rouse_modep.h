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
#include <general/mytools.h>

double calc_auto_corr(
    ReadDump::ExtraReadDump&,
    std::vector<std::vector<Eigen::Vector3d>*>&,
    const int&,
    const double&,
    const int&,
    const int&,
    const int&,
    const int&
);

void calc_Xp(
    ReadDump::ExtraReadDump&,
    std::vector<std::vector<Eigen::Vector3d>*>&,
    const int&,
    const int&,
    const int&
);
