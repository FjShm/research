#ifndef __MYTOOLS_H__
#define __MYTOOLS_H__


#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>



namespace std
{
    // split
    std::vector<std::string> split(const std::string &line, const char &delim){
        std::vector<std::string> elems;
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, delim))
            if (!item.empty())
                elems.push_back(item);
        return elems;
    }

    // linspace
    std::vector<double> linspace(const double &min, const double &max, const int &num){
        std::vector<double> ret(num);
        double dx = (max - min) / ((double)num - 1.);
        for (size_t i = 0; i < num; i++) ret[i] = min + dx*(double)i;
        return ret;
    }

    // count_rows
    int count_rows(const std::string &fname, const std::string &contains=""){
        std::ifstream f(fname);
        std::string buf;
        int counter = 0;
        while (std::getline(f, buf))
            if (contains == "" || buf.find(contains) != std::string::npos)
                counter++;
        return counter;
    }
}


#endif
