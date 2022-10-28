#include "r2n.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>();
    int N = param["beads_per_chain"].as<int>();
    int M = param["num_chain"].as<int>();

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(ipath, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------

    ReadDump::ReadDump rd(ipath);
    Eigen::VectorXd R2_n(N+1);

    int count = 0;
    while(rd.read_1frame()){
        rd.header_validation("id", "xu", "yu", "zu");
        compute_R2_n(rd, N, M, R2_n);
        ++count;
        ++show_progress;
    }
    R2_n /= (double) count;

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (int n = 2; n <= N; n++)
        out << n-1 << ' ' << R2_n(n) << std::endl;
}


void compute_R2_n(ReadDump::ReadDump &rd, int N, int M, Eigen::VectorXd& R2_n){
    Eigen::VectorXd R2_n_tmp(N+1);
    std::vector<Eigen::Vector3d> coordinations;
    rd.join_3columns(coordinations, "xu", "yu", "zu");
    Eigen::Vector3d dpos;
    // N-n+1: # of pairs for each chains
    // n-1: # of bonds
    for (int n = 2; n <= N; n++){
        for (int m = 0; m < M; m++){
            for (int start_id = 0; start_id <= N-n; start_id++){
                int end_id = start_id + n - 1;
                dpos = coordinations[m*N + start_id] - coordinations[m*N + end_id];
                R2_n_tmp(n) += dpos.squaredNorm();
            }
        }
        R2_n_tmp(n) /= (M * (N - n + 1)) * (n-1);
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
