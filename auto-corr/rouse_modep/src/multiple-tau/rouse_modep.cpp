#include "rouse_modep.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("ete_auto_corr.log");
    const double dt = param["dt_fs"].as<double>();
    const int timesteps_1frame = param["timesteps_1frame"].as<int>();
    const int total_frames_dump = param["total_frames_dump"].as<int>(-1);
    const int p = param["mode-p"].as<int>(1);
    const int N = param["beads_per_chain"].as<int>();
    const int M = param["num_chain"].as<int>();
    const bool normalize_autocorr = param["normalize_autocorr"].as<bool>(false);
    constexpr double fs2ns = 1.e-6;

    ReadDump::ExtraReadDump rd(ipath);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(total_frames_dump);
    // -------------------------------

    multipletau::correlator corr;
    while (rd.read_1frame()){
        //rd.header_validation("id", "xu", "yu", "zu");
        Eigen::Vector3d xp(0., 0., 0.);
        calc_Xp(xp, rd, p, N, M);
        corr(xp);
        ++show_progress;
    }

    std::vector<double> t = corr.get_time_vec();
    std::vector<double> f = corr.get_corr_vec();

    // output
    double t_scale = dt * (double)timesteps_1frame * fs2ns;
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (size_t i = 0; i < t.size(); i++){
        double dat = normalize_autocorr ?
            f[i]/f[0] : f[i];
        out << t[i]*t_scale << " " << dat << std::endl;
    }
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
    for (int m = 0; m < M; m++){
        Eigen::Vector3d xp_tmp(0., 0., 0.);
        _calc_Xp(xp_tmp, coods, p, N, m*N);
        xp += xp_tmp;
    }
    xp /= (double)M;
}
