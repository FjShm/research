#include "r2n.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("LOG.out");
    const int dump_frames = param["dump_frames"].as<int>();
    int N = param["beads_per_chain"].as<int>();
    int M = param["num_chain"].as<int>();

    // -------------------------------
    // max loop
    boost::progress_display show_progress(dump_frames);

    // -------------------------------

    ReadDump::ReadDump rd(ipath);
    Eigen::VectorXd R2_n(N-1);

    int count = 0;
    while(rd.read_1frame()){
        rd.header_validation("id", "xu", "yu", "zu");
        compute_R2_n(rd, N, M, R2_n);
        ++count;
        ++show_progress;
    }
    R2_n /= (double)count;

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (int bond = 1; bond < N; bond++)
        out << bond << ' ' << R2_n(bond-1) << std::endl;
}


void compute_R2_n(ReadDump::ReadDump &rd, int N, int M, Eigen::VectorXd& R2_n){
    Eigen::VectorXd R2_n_tmp(N-1);
    std::vector<Eigen::Vector3d> coordinations;
    rd.join_3columns(coordinations, "xu", "yu", "zu");
    // N-n+1: # of pairs for each chains
    // b: # of bonds
    for (int b = 1; b < N; b++){
        for (int m = 0; m < M; m++){
            for (int start_id = 0; start_id < N-b; start_id++){
                int end_id = start_id + b;
                R2_n_tmp(b-1) += (coordinations[m*N + start_id] - coordinations[m*N + end_id]).squaredNorm();
            }
        }
        R2_n_tmp(b-1) /= (M * (N - b)) * b;
    }
    R2_n += R2_n_tmp;
}

