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


namespace ReadDump
{
    class ReadDump {
        public:
            int timestep, num_atoms, num_frames;
            Eigen::Vector3d cellbox_a, cellbox_b, cellbox_c, cellbox_origin;
            Eigen::MatrixXd *atoms_all_data;
            std::map<std::string, int> *header_map;

            ReadDump() {init();}
            ReadDump(const std::string &arg_ipath) : ipath(arg_ipath){
                fileopen(ipath);
                init();
            }
            ~ReadDump() {clear();}

            void fileopen(const std::string &arg_ipath){
                dump.open(arg_ipath);
            }

            void clear(){
                dump.close();
                header_map->clear();
            }

            void join_3columns(
                std::vector<Eigen::Vector3d> &joined,
                std::string c1,
                std::string c2,
                std::string c3
            ){
                if (joined.size() != num_atoms) joined.resize(num_atoms);
                for (int i = 0; i < num_atoms; i++){
                    joined[i](0) = atoms_all_data->coeff(i, header_map->at(c1));
                    joined[i](1) = atoms_all_data->coeff(i, header_map->at(c2));
                    joined[i](2) = atoms_all_data->coeff(i, header_map->at(c3));
                }
            }

            bool read_1frame(){
                if (!all_frames_loaded){
                    return _read_1frame();
                } else {
                    return change_now_frame(1);
                }
            }

            void read_all_frames(){
                std::cout << ipath << " : now loading...\n";
                while(_read_1frame()){
                    std::cout << "\r>>> timestep: " + std::to_string(timestep);
                    timestep_v.push_back(timestep);
                    num_atoms_v.push_back(num_atoms);
                    ca_v.push_back(cellbox_a);
                    cb_v.push_back(cellbox_b);
                    cc_v.push_back(cellbox_c);
                    co_v.push_back(cellbox_origin);
                    atoms_all_data_v.push_back(atoms_all_data);
                    header_map_v.push_back(header_map);
                }
                std::cout << std::endl << "done" << std::endl;
                all_frames_loaded = true;
            }

            void header_validation(const std::vector<std::string> &headers){
                bool abort = false;
                for (std::string header : headers){
                    if (header_map->count(header) == 0){
                        std::cout << "'" << header << "' is not exist in dump file ATOMS.\n";
                        abort = true;
                    }
                }
                if (abort) std::exit(EXIT_FAILURE);
            }

            template<class... T> void header_validation(T... headers){
                bool abort = false;
                for (std::string header : std::initializer_list<std::string>{headers...}){
                    if (header_map->count(header) == 0){
                        std::cout << "'" << header << "' is not exist in dump file ATOMS.\n";
                        abort = true;
                    }
                }
                if (abort) std::exit(EXIT_FAILURE);
            }

            void set_want_frames(const std::vector<double> &read_ratio){
                for (double rr : read_ratio)
                    want_timesteps.push_back(search_nearest_timestep(rr));
            }

            void set_want_frames(const std::vector<int> &read_timestep){
                for (int ts : read_timestep)
                    want_timesteps.push_back(ts);
            }

            int search_nearest_timestep(const double &ratio){
                if (ratio < 0. || 1. < ratio){
                    std::cout << "Invalid ratio was given.\n" << ratio
                        << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                if (!all_frames_loaded){
                    std::cout << "Call read_all_frames before calling the "
                        << "search_nearest_timestep (read_dump.h)\n";
                    std::exit(EXIT_FAILURE);
                }
                int timestep_max = std::vec_max(timestep_v);
                std::vector<double> diffs(timestep_v.size());
                for (size_t i = 0; i < diffs.size(); i++)
                    diffs[i] = std::abs((double)timestep_v[i]/(double)timestep_max - ratio);
                return timestep_v[std::vec_minid(diffs)];
            }

            bool check_if_wanted_frame(){
                if (!all_frames_loaded) return true;
                if (want_timesteps.size() == 0) return true;
                std::vector<int>::iterator itr
                    = std::find(want_timesteps.begin(), want_timesteps.end(), timestep);
                if (itr == want_timesteps.end()){
                    return false;
                } else {
                    return true;
                }
            }



        private:
            int line_number = 0;
            int now_frame = -1;
            std::string ipath, tmp;
            std::ifstream dump;
            bool all_frames_loaded = false;

            // 全frameのデータを格納するvector
            std::vector<int> timestep_v, num_atoms_v, want_timesteps;
            std::vector<Eigen::Vector3d> ca_v, cb_v, cc_v, co_v;
            std::vector<Eigen::MatrixXd*> atoms_all_data_v;
            std::vector< std::map<std::string, int>* > header_map_v;

            void init(){
                num_frames = 0;
            }

            void read_timestep(){
                dump >> timestep;
                std::getline(dump, tmp);
                line_number++;
            }

            void read_number_of_atoms(){
                dump >> num_atoms;
                std::getline(dump, tmp);
                line_number++;
            }

            void read_box(std::string boxtype){
                if (boxtype == "triclinic"){
                    double xlo_b,xhi_b,xy, ylo_b,yhi_b,xz, zlo_b,zhi_b,yz;
                    dump >> xlo_b >> xhi_b >> xy
                         >> ylo_b >> yhi_b >> xz
                         >> zlo_b >> zhi_b >> yz;
                    dumpcell_to_vector_converter(
                        xlo_b, xhi_b, xy, ylo_b, yhi_b, xz, zlo_b, zhi_b, yz);
                    line_number += 3;
                } else if (boxtype == "orthogonal"){
                    double xlo,xhi, ylo,yhi, zlo,zhi;
                    dump >> xlo >> xhi
                         >> ylo >> yhi
                         >> zlo >> zhi;
                    cellbox_origin << xlo, ylo, zlo;
                    cellbox_a << xhi-xlo, 0., 0.;
                    cellbox_b << 0., yhi-ylo, 0.;
                    cellbox_c << 0., 0., zhi-zlo;
                    line_number += 3;
                } else {
                    std::cout << "Error: undefined box type\n"
                        << "line " << line_number << std::endl;
                }
                std::getline(dump, tmp);
            }

            void read_atoms(int num_column){
                if (header_map->count("id") == 0){
                    for (int r = 0; r < num_atoms; r++){
                        for (int c = 0; c < num_column; c++){
                            dump >> (*atoms_all_data)(r, c);
                        }
                    }
                } else {
                    std::vector<double> buf(num_column);
                    for (int r = 0; r < num_atoms; r++){
                        for (int c = 0; c < num_column; c++){
                            dump >> buf[c];
                        }
                        for (int c = 0; c < num_column; c++){
                            (*atoms_all_data)((int)buf[header_map->at("id")]-1, c) = buf[c];
                        }
                    }
                }
                std::getline(dump, tmp);
                line_number += num_atoms;
            }

            void dumpcell_to_vector_converter(
                    double &xlo_b, double &xhi_b, double &xy,
                    double &ylo_b, double &yhi_b, double &xz,
                    double &zlo_b, double &zhi_b, double &yz
                    ){
                double xlo,xhi, ylo,yhi, zlo,zhi;
                xlo = xlo_b - std::min({0., xy, xz, xy + xz});
                xhi = xhi_b - std::max({0., xy, xz, xy + xz});
                ylo = ylo_b - std::min({0., yz});
                yhi = yhi_b - std::max({0., yz});
                zlo = zlo_b;
                zhi = zhi_b;
                cellbox_origin << xlo, ylo, zlo;
                cellbox_a << xhi - xlo, 0., 0.;
                cellbox_b << xy, yhi - ylo, 0.;
                cellbox_c << xz, yz, zhi - zlo;
            }

            bool _read_1frame(){
                std::string row;
                atoms_all_data = new Eigen::MatrixXd;
                header_map = new std::map<std::string, int>;
                while(std::getline(dump, row)){
                    line_number++;
                    std::vector<std::string> splited_row = std::split(row, ':');
                    if (splited_row[0] != "ITEM"){
                        std::cout << "Error: Invalid dumpfile format.\n"
                            << "line " << line_number << ": " << row << std::endl;
                        std::exit(EXIT_FAILURE);
                    }

                    if (splited_row[1] == " TIMESTEP"){
                        read_timestep();
                        num_frames++;
                    } else if (splited_row[1] == " NUMBER OF ATOMS"){
                        read_number_of_atoms();
                    } else if (splited_row[1] == " BOX BOUNDS xy xz yz pp pp pp"){
                        read_box("triclinic");
                    } else if (splited_row[1] == " BOX BOUNDS xx yy zz pp pp pp"){
                        read_box("orthogonal");
                    } else if (splited_row[1] == " BOX BOUNDS pp pp pp"){
                        read_box("orthogonal");
                    } else {
                        std::vector<std::string> atoms_header = split(splited_row[1], ' ');
                        if (atoms_header[0] != "ATOMS"){
                            std::cout << "Error: Invalid dumpfile format.\n"
                                << "ITEM: ATOMS must be come\n"
                                << "line " << line_number << ": " << row << std::endl;
                            std::exit(EXIT_FAILURE);
                        }

                        // header index をマッピング
                        for (size_t i = 1; i < atoms_header.size(); i++)
                            header_map->insert(std::make_pair(atoms_header[i], i-1));

                        // ATOMS 読み込み
                        int num_column = atoms_header.size() - 1;
                        atoms_all_data->conservativeResize(num_atoms, num_column);
                        read_atoms(num_column);

                        return true; // ATOMSの最後が1frameの最後
                    }
                }
                return false; // ファイル末尾の場合はここ
            }

            bool change_now_frame(int frame, bool absolute = false){
                if (absolute){
                    now_frame = frame;
                } else {
                    now_frame += frame;
                }
                if (now_frame < 0 || num_frames <= now_frame) return false;
                timestep = timestep_v[now_frame];
                num_atoms = num_atoms_v[now_frame];
                cellbox_a = ca_v[now_frame];
                cellbox_b = cb_v[now_frame];
                cellbox_c = cc_v[now_frame];
                cellbox_origin = co_v[now_frame];
                atoms_all_data = atoms_all_data_v[now_frame];
                header_map = header_map_v[now_frame];
                return true;
            }

    }; // class ReadDump


    class ExtraReadDump : public ReadDump
    {
        using ReadDump::ReadDump;

        public:
            template<class... T> void add_column_if_not_exist(std::string colname, T... args){
                if (colname == "mol"){
                    add_mol(args...);
                } else {
                    std::cout << "Invalid column name: " << colname << std::endl
                        << "This message can be ignored but may cause an error.\n"
                        << "(add_column_if_not_exist, read_dump.h)\n";
                }
            }

            double max_of_col(std::string colname){
                int col = header_map->at(colname);
                double val = atoms_all_data->col(col).array().maxCoeff();
                return val;
            }

            double min_of_col(std::string colname){
                int col = header_map->at(colname);
                double val = atoms_all_data->col(col).array().minCoeff();
                return val;
            }

        private:
            void add_mol(int N, int M){
                if (header_map->count("mol") == 0){
                    int id = header_map->at("id");
                    int mol = atoms_all_data->cols();
                    header_map->insert(std::make_pair("mol", mol));
                    atoms_all_data->conservativeResize(num_atoms, mol+1);
                    if (N != -1){
                        for (int i = 0; i < num_atoms; i++)
                            (*atoms_all_data)(i, mol) = ((int)atoms_all_data->coeff(i, id) - 1) / N + 1;
                    } else if (M != -1){
                        N = num_atoms / M;
                        for (int i = 0; i < num_atoms; i++)
                            (*atoms_all_data)(i, mol) = ((int)atoms_all_data->coeff(i, id) - 1) / N + 1;
                    } else {
                        std::cout << "Since there is no 'mol' in ATOMS in the dump file,"
                            << " it is necessary to write N or M in the input file.\n";
                        std::exit(EXIT_FAILURE);
                    }
                }
            }
    }; // class ExtraReadDump
}

#endif
