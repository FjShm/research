#include "compare_with_ideal_chain.h"



int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>();
    const int N = param["N"].as<int>(-1);
    const int M = param["M"].as<int>(-1);

    // -------------------------------
    // max loop
    int max_loop = std::count_rows(ipath, "TIMESTEP");
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    double R2(0), R4(0), Rgcm(0);
    ReadDump::ExtraReadDump rd(ipath);

    while(rd.read_1frame()){
        rd.header_validation("id", "xu", "yu", "zu");
        rd.add_column_if_not_exist("mol", N, M);
        compute(rd, R2, R4, Rgcm);
    
        // update progress bar
        ++show_progress;
    }
    R2 /= (double)rd.num_frames;
    R4 /= (double)rd.num_frames;
    Rgcm /= (double)rd.num_frames;

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    out << "<R4>/<(R2)>^2 = " << R4/(R2*R2) << std::endl
        << "<R2>/<RG2>    = " << R2/Rgcm << std::endl;
}


void compute(ReadDump::ExtraReadDump& rd, double& R2, double& R4, double& Rgcm){
    std::vector<Eigen::Vector3d> pos;
    rd.join_3columns(pos, "xu", "yu", "zu");
    double R2_tmp(0), R4_tmp(0), Rgcm_tmp(0);

    int M = rd.max_of_col("mol");
    int N = rd.num_atoms / M;
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

