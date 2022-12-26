#include "rouse_modep.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("ete_auto_corr.log");
    const double dt = param["dt_fs"].as<double>();
    const int ave_frames = param["max_frames_to_average"].as<int>(100);
    const int total_frames_dump = param["total_frames_dump"].as<int>(-1);
    const int p = param["mode-p"].as<int>(1);
    const int N = param["beads_per_chain"].as<int>();
    const int M = param["num_chain"].as<int>();
    constexpr double fs2ns = 1.e-6;

    ReadDump::ExtraReadDump rd(ipath);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(total_frames_dump);
    // -------------------------------

    while (rd.read_1frame()){
        //rd.header_validation("id", "xu", "yu", "zu");
        Eigen::Vector3d xp(0., 0., 0.);
        calc_Xp(xp, rd, p, N, M);
        ++show_progress;
    }

    // output
    //std::ofstream out{opath, std::ios::out | std::ios::trunc};
    //for (int i = 0; i < rd.num_frames; i++)
    //    out << time[i] << " " << auto_corr[i] << std::endl;
}


void _calc_Xp(
    Eigen::Vector3d &xp_tmp,
    std::vector<Eigen::Vector3d> &coods,
    const int &p,
    const int &N,
    const int &start_id
){
    const double sqrt_2_N = std::sqrt(2./(double)N);

    for (int n = start_id; n < start_id + N; n++)
        xp_tmp += std::cos(((double)n + 0.5)*(double)p*M_PI/(double)N) * coods[n];
    xp_tmp *= sqrt_2_N;
}

void calc_Xp(
    Eigen::Vector3d &xp,
    ReadDump::ExtraReadDump &rd,
    const int &p,
    const int &N,
    const int &M
){
    std::vector<Eigen::Vector3d> coods;
    rd.join_3columns(coods, "xu", "yu", "zu");
    for (int m; m < M; m++){
        Eigen::Vector3d xp_tmp(0., 0., 0.);
        _calc_Xp(xp_tmp, coods, p, N, m*N);
    }
    xp /= (double)M;
}
