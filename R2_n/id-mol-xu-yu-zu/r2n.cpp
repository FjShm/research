#include "r2n.h"


int main(int argc, char* argv[]){
    if (argc < 2){
        std::cout << "Error: Execute with the input yaml file path as an argument.\n"
            << "ex.) $ ./a.out param.yaml" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>();
    const std::string col2 = param["col2"].as<std::string>();
    int N = param["beads_per_chain"].as<int>();
    int M = param["num_chain"].as<int>();

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(ipath, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------

    ReadDump::ReadDump rd(ipath);
    int NM, count(0), timestep, beads_total;
    NM = N * M;
    Eigen::VectorXd R2_n(N+1);
    R2_n *= 0;
    std::ifstream in{ipath};

    //while(std::getline(in, line)){
    //    if (line == "ITEM: TIMESTEP"){
    //        in >> timestep;
    //    } else if (line == "ITEM: NUMBER OF ATOMS"){
    //        in >> beads_total;
    //        if (beads_total != NM){
    //            std::cout << "N, M is wrong. N*M = " << NM
    //            << ", total number of beads = " << beads_total << std::endl;
    //            std::exit(EXIT_FAILURE);
    //        }
    //    } else if (line == "ITEM: ATOMS id mol xu yu zu" ||
    //            line == "ITEM: ATOMS id type xu yu zu"){
    //        if (col2 == "mol"){
    //            compute_R2_n(in, N, M, NM, R2_n, "mol");
    //        } else if (col2 == "type"){
    //            compute_R2_n(in, N, M, NM, R2_n, "type");
    //        }
    //        count++;

    //        // update progress bar
    //        ++show_progress;
    //    } else continue;
    //}
    while(rd.read_1frame()){
        rd.header_validation("id", "xu", "yu", "zu");
        int id, mol;
        if (rd.header_map->count("mol") == 0){
            id = rd.header_map->at("id");
            mol = rd.atoms_all_data->cols();
            rd.header_map->insert(std::make_pair("mol", mol));
            rd.atoms_all_data->conservativeResize(rd.num_atoms, mol+1);
            if (N != -1){
                for (int i = 0; i < rd.num_atoms; i++)
                    (*rd.atoms_all_data)(i, mol) = ((int)rd.atoms_all_data->coeff(i, id) - 1) / N + 1;
            } else if (M != -1){
                N = rd.num_atoms / M;
                for (int i = 0; i < rd.num_atoms; i++)
                    (*rd.atoms_all_data)(i, mol) = ((int)rd.atoms_all_data->coeff(i, id) - 1) / N + 1;
            } else {
                std::cout << "Since there is no 'mol' in ATOMS in the dump file,"
                    << " it is necessary to write N or M in the "
                    << argv[1] << ".\n";
                std::exit(EXIT_FAILURE);
            }
        } else {
            mol = rd.header_map->at("mol");
        }
        compute_R2_n(rd, N, M, NM, R2_n);
        ++count;
        ++show_progress;
    }
    R2_n /= (double) count;

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (int n = 2; n <= N; n++){
        out << n-1 << ' ' << R2_n(n) << std::endl;
    }
}


void compute_R2_n(ReadDump::ReadDump &rd, int N, int M, int NM, Eigen::VectorXd& R2_n){
    Eigen::VectorXd R2_n_tmp(N+1);
    //int id, mol, type;
    //for (int i = 0; i < NM; i++){
    //    in >> id;
    //    if (col2 == "mol"){
    //        in >> mol;
    //    } else if (col2 == "type"){
    //        in >> type;
    //        mol = (id - 1)/N + 1;
    //    } else {
    //        std::cout << "warning: The lammpstrj format corresponds "
    //            << "to id-mol-xu-yu-zu or id-type-xu-yu-zu.\n";
    //        return;
    //    }
    //    int idx = (id - 1) % N;
    //    in >> posx[mol-1][idx] >> posy[mol-1][idx] >> posz[mol-1][idx];
    //}

    std::vector<Eigen::Vector3d> coordinations;
    rd.join_3columns(coordinations, "xu", "yu", "zu");
    Eigen::Vector3d dpos;
    for (int n = 2; n <= N; n++){
        for (int m = 0; m < M; m++){
            for (int start_id = 0; start_id <= N-n; start_id++){
                int end_id = start_id + n - 1;
                dpos = coordinations[m*N + start_id] - coordinations[m*N + end_id];
                R2_n_tmp(n) += dpos.squaredNorm();
            }
        }
        R2_n_tmp(n) /= (M * (N - n + 1)) * (n-1);
        // M: chains
        // N-n+1: # of pairs for each chains
        // n-1: # of bonds
    }
    R2_n += R2_n_tmp;
}

void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    while(std::getline(in, row))
        if (row == "ITEM: TIMESTEP")
            max_loop++;
}
