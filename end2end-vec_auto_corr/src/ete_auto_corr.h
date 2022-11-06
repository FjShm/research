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

double calc_auto_corr(ReadDump::ExtraReadDump &, int, int, int, int, int);
