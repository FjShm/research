#ifndef __WRITE_DUMP_H__
#define __WRITE_DUMP_H__

#include <cstdlib>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ios>
#include <eigen3/Eigen/Dense>
#include <general/mytools.h>


namespace WriteDump
{
    class WriteDump {
        public:
            int *timestep, *num_atoms;
            Eigen::Vector3d *cellbox_a, *cellbox_b, *cellbox_c, *cellbox_origin;
            Eigen::MatrixXd *atoms_all_data;
            std::map<std::string, int> *header_map;

            WriteDump() {init();}
            WriteDump(const std::string &arg_opath) : opath(arg_opath){
                open(opath);
                init();
            }
            ~WriteDump() {clear();}

            void open(const std::string &arg_opath){
                dump.open(arg_opath, std::ios::out | std::ios::trunc);
            }

            void clear(){
                dump.close();
                init();
                delete dumpcell;
            }

            void write_1frame(){
                if (!header_set){
                    std::vector<std::string> headers
                        = keys_of_map_sorted_by_value(header_map);
                    set_headers(headers);
                }
                _write_1frame();
            }

            void set_headers(const std::vector<std::string> &headers){
                if (header_map == nullptr){
                    std::cout << "Set header_map first.\n";
                    std::exit(EXIT_FAILURE);
                }
                for (std::string s : headers){
                    if (header_validation(s)){
                        write_col_name.push_back(s);
                    } else {
                        std::cout << "'" << s << "' is ignored because it is not exist "
                            << "in dump file header_map.\n";
                    }
                }
                header_set = true;
            }

            template<class... T> void set_headers(T... args){
                if (header_map == nullptr){
                    std::cout << "Set header_map first.\n";
                    std::exit(EXIT_FAILURE);
                }
                for (std::string s : std::initializer_list<std::string>{args...}){
                    if (header_validation(s)){
                        write_col_name.push_back(s);
                    } else {
                        std::cout << "'" << s << "' is ignored because it is not exist "
                            << "in dump file header_map.\n";
                    }
                }
                header_set = true;
            }


        private:
            std::string opath;
            std::ofstream dump;
            Eigen::Matrix3d *dumpcell = new Eigen::Matrix3d;
            std::vector<std::string> write_col_name;
            bool header_set;

            void init(){
                timestep = nullptr;
                num_atoms = nullptr;
                cellbox_a = nullptr;
                cellbox_b = nullptr;
                cellbox_c = nullptr;
                cellbox_origin = nullptr;
                atoms_all_data = nullptr;
                header_map = nullptr;
                write_col_name.clear();
                header_set = false;
            }

            void write_timestep(){
                if (timestep == nullptr){
                    std::cout << "'timestep' is not set.\n";
                    std::exit(EXIT_FAILURE);
                }
                dump << "ITEM: TIMESTEP\n"
                    << *timestep << std::endl;
            }

            void write_number_of_atoms(){
                if (num_atoms == nullptr){
                    std::cout << "'num_atoms' is not set.\n";
                    std::exit(EXIT_FAILURE);
                }
                dump << "ITEM: NUMBER OF ATOMS\n"
                    << *num_atoms << std::endl;
            }

            void write_box(){
                if (cellbox_a == nullptr ||
                    cellbox_b == nullptr ||
                    cellbox_c == nullptr ||
                    cellbox_origin == nullptr)
                {
                    std::cout << "'cellbox_*' is not set.\n";
                    std::exit(EXIT_FAILURE);
                }
                vector_to_dumpcell_converter();
                dump << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";
                std::ios::fmtflags flagsSaved = dump.flags();
                dump << std::scientific;
                for (int i = 0; i < 3; i++){
                    for (int j = 0; j < 3; j++)
                        dump << dumpcell->coeff(i, j) << " ";
                    dump << std::endl;
                }
                dump.flags(flagsSaved);

            }

            void write_atoms(){
                dump << "ITEM: ATOMS";
                for (std::string header : write_col_name)
                    dump << " " + header;
                dump << std::endl;
                for (int rowidx = 0; rowidx < *num_atoms; rowidx++){
                    for (std::string col : write_col_name){
                        dump << atoms_all_data->coeff(rowidx, header_map->at(col)) << " ";
                    }
                    dump << std::endl;
                }
            }

            void vector_to_dumpcell_converter(){
                double xlo,xhi,xy, ylo,yhi,xz, zlo,zhi,yz;
                xlo = cellbox_origin->coeff(0);
                ylo = cellbox_origin->coeff(1);
                zlo = cellbox_origin->coeff(2);
                xy = cellbox_b->coeff(0);
                xz = cellbox_c->coeff(0);
                yz = cellbox_c->coeff(1);
                xhi = cellbox_a->coeff(0) + xlo;
                yhi = cellbox_b->coeff(1) + ylo;
                zhi = cellbox_c->coeff(2) + zlo;

                double xlo_b,xhi_b, ylo_b,yhi_b, zlo_b,zhi_b;
                xlo_b = xlo + std::min({0., xy, xz, xy + xz});
                xhi_b = xhi + std::max({0., xy, xz, xy + xz});
                ylo_b = ylo + std::min({0., yz});
                yhi_b = yhi + std::max({0., yz});
                zlo_b = zlo;
                zhi_b = zhi;

                *dumpcell << xlo_b, xhi_b, xy,
                         ylo_b, yhi_b, xz,
                         zlo_b, zhi_b, yz;
            }

            void _write_1frame(){
                write_timestep();
                write_number_of_atoms();
                write_box();
                write_atoms();
                init();
            }

            template<class... T> bool header_validation(T... headers){
                bool flg = true;
                for (std::string header : std::initializer_list<std::string>{headers...}){
                    if (header_map->count(header) == 0){
                        flg = false;
                    }
                }
                return flg;
            }

            std::vector<std::string> keys_of_map_sorted_by_value(std::map<std::string, int> *map){
                std::vector<std::string> keys(map->size(), "null");
                int counter = 0;
                for (auto const &[key, val]: *map){
                    if (keys[val] == "null"){
                        keys[val] = key;
                        counter++;
                    } else {
                        std::cout << "'" << key << "' and '" << keys[val] << "' have "
                            << "same value at header_map.\n"
                            << "'" << keys[val] << "' is not exported to dumpfile.\n";
                    }
                }
                if (map->size() != counter){
                    keys.resize(counter);
                }
                return keys;
            }

    }; // class WriteDump
}

#endif
