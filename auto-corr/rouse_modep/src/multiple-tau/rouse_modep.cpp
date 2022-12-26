#include "rouse_modep.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("Xp_auto_corr.log");
    const double dt = param["dt_fs"].as<double>();
    const int timesteps_1frame = param["timesteps_1frame"].as<int>();
    const int total_frames_dump = param["total_frames_dump"].as<int>(-1);
    const int p = param["mode-p"].as<int>(1);
    const int N = param["N"].as<int>();
    const int M = param["M"].as<int>();
    const bool normalize_autocorr = param["normalize_autocorr"].as<bool>(false);
    constexpr double fs2ns = 1.e-6;

    ReadDump::ExtraReadDump rd(ipath);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(total_frames_dump);
    // -------------------------------

    multipletau::correlator corr[M];
    while (rd.read_1frame()){
        //rd.header_validation("id", "xu", "yu", "zu");
        Eigen::Vector3d xp;
        std::vector<Eigen::Vector3d> coods;
        rd.join_3columns(coods, "xu", "yu", "zu");
        for (int m = 0; m < M; m++){
            _calc_Xp(xp, coods, p, N, m*N);
            corr[m](xp);
        }
        ++show_progress;
    }

    // average <Xp>
    double t_scale = dt * (double)timesteps_1frame * fs2ns;
    std::vector<double> t = corr[0].get_time_vec();
    std::vector<double> f = corr[0].get_corr_vec();
    for (int m = 1; m < M; m++){
        std::vector<double> f_ = corr[m].get_corr_vec();
        for (size_t i = 0; i < f_.size(); i++)
            f[i] += f_[i];
    }
    for (size_t i = 0; i < f.size(); i++)
        f[i] /= (double)M;


    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (size_t i = 0; i < t.size(); i++){
        double dat = normalize_autocorr ? f[i]/f[0] : f[i];
        out << t[i]*t_scale << " " << dat << std::endl;
    }
}


void _calc_Xp(
    Eigen::Vector3d &xp,
    std::vector<Eigen::Vector3d> &coods,
    const int &p,
    const int &N,
    const int &start_id
){
    const double sqrt_2_N = std::sqrt(2./(double)N);
    xp = Eigen::Vector3d::Zero();

    for (int n = start_id; n < start_id + N; n++)
        xp += std::cos(((double)n + 0.5)*(double)p*M_PI/(double)N) * coods[n];
    xp *= sqrt_2_N;
}

