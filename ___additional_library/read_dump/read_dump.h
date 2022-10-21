#include <cstdlib>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <eigen3/Eigen/Dense>


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
                if (timestep_v.size() == 0){
                    // read_all_framesが実行されていない場合
                    return _read_1frame();
                } else {
                    if (now_frame >= num_frames) return false;
                    timestep = timestep_v[now_frame];
                    num_atoms = num_atoms_v[now_frame];
                    cellbox_a = ca_v[now_frame];
                    cellbox_b = cb_v[now_frame];
                    cellbox_c = cc_v[now_frame];
                    cellbox_origin = co_v[now_frame];
                    atoms_all_data = &atoms_all_data_v[now_frame];
                    header_map = &header_map_v[now_frame];
                    now_frame++;
                    return true;
                }
            }

            void read_all_frames(){
                std::cout << ipath << " : now loading...\n";
                while(_read_1frame()){
                    num_frames++;
                    std::cout << "\rtimestep: " << timestep;
                    timestep_v.push_back(timestep);
                    num_atoms_v.push_back(num_atoms);
                    ca_v.push_back(cellbox_a);
                    cb_v.push_back(cellbox_b);
                    cc_v.push_back(cellbox_c);
                    co_v.push_back(cellbox_origin);
                    atoms_all_data_v.push_back(_atoms_all_data);
                    header_map_v.push_back(_header_map);
                }
                std::cout << std::endl << "done" << std::endl;
            }

            void header_validation(std::vector<std::string> &headers){
                bool abort = false;
                for (std::string header : headers){
                    if (header_map->count(header) == 0){
                        std::cout << "'" << header << "' is not exist in dump file ATOMS.\n";
                        abort = true;
                    }
                }
                if (abort) std::exit(EXIT_FAILURE);
            }

            void header_validation(std::string header){
                if (header_map->count(header) == 0){
                    std::cout << "'" << header << "' is not exist in dump file ATOMS.\n";
                    std::exit(EXIT_FAILURE);
                }
            }


        private:
            int line_number = 0;
            int now_frame = 0;
            std::string ipath, tmp;
            std::ifstream dump;
            Eigen::MatrixXd _atoms_all_data;
            std::map<std::string, int> _header_map;

            // 全frameのデータを格納するvector
            std::vector<int> timestep_v, num_atoms_v;
            std::vector<Eigen::Vector3d> ca_v, cb_v, cc_v, co_v;
            std::vector<Eigen::MatrixXd> atoms_all_data_v;
            std::vector< std::map<std::string, int> > header_map_v;

            std::vector<std::string> split(const std::string &s, char delim){
                std::vector<std::string> elems;
                std::stringstream ss(s);
                std::string item;
                while (std::getline(ss, item, delim)) {
                    if (!item.empty()) {
                        elems.push_back(item);
                    }
                }
                return elems;
            }

            void init(){
                atoms_all_data = &_atoms_all_data;
                header_map = &_header_map;
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
                    std::getline(dump, tmp);
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
                while(std::getline(dump, row)){
                    line_number++;
                    std::vector<std::string> splited_row = split(row, ':');
                    if (splited_row[0] != "ITEM"){
                        std::cout << "Error: Invalid dumpfile format.\n"
                            << "line " << line_number << ": " << row << std::endl;
                        std::exit(EXIT_FAILURE);
                    }

                    if (splited_row[1] == " TIMESTEP"){
                        read_timestep();
                    } else if (splited_row[1] == " NUMBER OF ATOMS"){
                        read_number_of_atoms();
                    } else if (splited_row[1] == " BOX BOUNDS xy xz yz pp pp pp"){
                        read_box("triclinic");
                    } else if (splited_row[1] == " BOX BOUNDS xx yy zz pp pp pp"){
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
    };
}

