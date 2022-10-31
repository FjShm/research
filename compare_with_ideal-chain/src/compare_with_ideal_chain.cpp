#include "compare_with_ideal_chain.h"



int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>();
    const int N = param["N"].as<int>();
    const int M = param["M"].as<int>();

    const int NM = N * M;

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(ipath, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    std::string line;
    int count(0), timestep, beads_total;
    double R2(0), R4(0), Rgcm(0);
    std::ifstream in{ipath};

    while(std::getline(in, line)){
        if (line == "ITEM: TIMESTEP"){
            in >> timestep;
            continue;
        } else if (line == "ITEM: NUMBER OF ATOMS"){
            in >> beads_total;
            if (beads_total != NM){
                std::cout << "N, M is wrong. N*M = " << NM
                << ", total number of beads = " << beads_total << std::endl;
                return -1;
            }
            continue;
        } else if (line == "ITEM: BOX BOUNDS xy xz yz pp pp pp"){
            continue;
        } else if (line == "ITEM: ATOMS id xu yu zu"){
            compute_R2_n(in, N, M, NM, R2, R4, Rgcm);
            count++;
        } else continue;
    
        // update progress bar
        ++show_progress;
    }
    R2 /= (double)count;
    R4 /= (double)count;
    Rgcm /= (double)count;

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    out << "<R4/(R2)^2> = " << R4/(R2*R2) << std::endl;
    out << "<R2/RG2>    = " << R2/Rgcm << std::endl;
}


void compute_R2_n(std::ifstream& in, int N, int M, int NM, double& R2, double& R4, double& Rgcm){
    std::vector<Eigen::Vector3d> pos(NM);
    double R2_tmp(0), R4_tmp(0), Rgcm_tmp(0);
    int index;
    for (int i = 0; i < NM; i++){
        in >> index;
        in >> pos[index-1](0) >> pos[index-1](1) >> pos[index-1](2);
    }

    for (int i = 0; i < M; i++){
        int idx_head = i*N;
        int idx_tail = idx_head + N - 1;
        double tmp = (pos[idx_tail] - pos[idx_head]).squaredNorm();
        R2_tmp += tmp;
        R4_tmp += tmp * tmp;

        Eigen::Vector3d Rcmj;
        for (int j = idx_head; j <= idx_tail; j++){
            Rcmj += pos[j];
        }
        Rcmj /= (double)N;

        double RG_cm_tmp = 0;
        for (int j = idx_head; j <= idx_tail; j++){
            RG_cm_tmp += (pos[j] - Rcmj).squaredNorm();
        }
        Rgcm_tmp += RG_cm_tmp / (double)N;
    }

    R2 += R2_tmp / (double)M;
    R4 += R4_tmp / (double)M;
    Rgcm += Rgcm_tmp / (double)M;
}

void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    while(std::getline(in, row))
        if (row == "ITEM: TIMESTEP")
            max_loop++;
}
