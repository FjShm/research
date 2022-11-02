#ifndef __READ_DUMP_H__
#define __READ_DUMP_H__

#include <cstdlib>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
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
            }

            void write_1frame(){
                if (timestep == nullptr ||
                    num_atoms == nullptr ||
                    cellbox_a == nullptr ||
                    cellbox_b == nullptr ||
                    cellbox_c == nullptr ||
                    cellbox_origin == nullptr ||
                    atoms_all_data == nullptr ||
                    header_map == nullptr)
                {
                    std::cout << "The information required to output"
                        << "the lammpstrj file is missing.\n";
                    std::exit(EXIT_FAILURE);
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
                        write_col_idx.push_back(header_map->at(s));
                        write_col_name.push_back(s);
                    }
                }
            }

            template<class... T> void set_headers(T... args){
                if (header_map == nullptr){
                    std::cout << "Set header_map first.\n";
                    std::exit(EXIT_FAILURE);
                }
                for (std::string s : std::initializer_list<std::string>{args...}){
                    if (header_validation(s))
                        write_col_idx.push_back(header_map->at(s));
                        write_col_name.push_back(s);
                }
            }


        private:
            std::string opath;
            std::ofstream dump;
            Eigen::Matrix3d *dumpcell;
            std::vector<int> write_col_idx;
            std::vector<std::string> write_col_name;

            void init(){
                timestep = nullptr;
                num_atoms = nullptr;
                cellbox_a = nullptr;
                cellbox_b = nullptr;
                cellbox_c = nullptr;
                cellbox_origin = nullptr;
                atoms_all_data = nullptr;
                header_map = nullptr;
                write_col_idx.clear();
                write_col_name.clear();
            }

            void write_timestep(){
                dump << "ITEM: TIMESTEP\n";
            }

            void write_number_of_atoms(){
                dump << "ITEM: NUMBER OF ATOMS\n";
            }

            void write_box(){
                dump << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";
            }

            void write_atoms(){
                dump << "ITEM: ATOMS";
                for (std::string header : write_col_name)
                    dump << " " + header;
                dump << std::endl;
                for (int rowidx = 0; rowidx < *num_atoms; rowidx++){
                    for (int colidx : write_col_idx){
                        dump << atoms_all_data->coeff(rowidx, colidx) << " ";
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
                        std::cout << "'" << header << "' is ignored because it is not exist "
                            << "in dump file header_map.\n";
                        flg = false;
                    }
                }
                return flg;
            }

    }; // class WriteDump
}

#endif
