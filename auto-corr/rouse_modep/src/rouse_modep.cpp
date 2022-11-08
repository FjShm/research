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
    rd.read_all_frames(total_frames_dump);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(rd.num_frames);
    // -------------------------------

    std::vector<double> time(rd.num_frames), auto_corr(rd.num_frames);
    std::vector<std::vector<Eigen::Vector3d>*> Xp(rd.num_frames, nullptr);
    for (size_t fp = 0; fp < rd.num_frames; fp++){
        // calc time
        rd.jump_frames(0, true);
        double timestep0 = (double)rd.timestep;
        rd.jump_frames(fp);
        time[fp] = ((double)rd.timestep - timestep0)*dt*fs2ns;

        // calc auto corr
        rd.header_validation("id", "xu", "yu", "zu");
        int f0_max = (rd.num_frames-1) - fp;
        int ave_frames_ = std::min(ave_frames, f0_max+1);
        double df = ave_frames_ == 1 ? 0. : (double)f0_max / ((double)ave_frames_-1.);
        auto_corr[fp] =
            calc_auto_corr(rd, Xp, fp, df, p, ave_frames_, N, M);
        ++show_progress;
    }
    for (int i = 0; i < rd.num_frames; i++)
        delete Xp[i];

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (int i = 0; i < rd.num_frames; i++)
        out << time[i] << " " << auto_corr[i] << std::endl;
}


double calc_auto_corr(
    ReadDump::ExtraReadDump &rd,
    std::vector<std::vector<Eigen::Vector3d>*> &Xp,
    const int &fp,
    const double &df,
    const int &num_use_frames,
    const int &p,
    const int &N,
    const int &M
){
    double auto_corr = 0.;
    int counter = 0;
    for (int i = 0; i < num_use_frames; i++){
        int f0 = std::round(i*df);
        rd.jump_frames(f0, true);
        calc_Xp(rd, Xp, p, N, M);
        rd.jump_frames(fp);
        calc_Xp(rd, Xp, p, N, M);
        for (int m = 0; m < M; m++)
            auto_corr += (Xp[f0+fp]->at(m)).dot(Xp[f0]->at(m)) / (Xp[f0]->at(m)).squaredNorm();
    }
    return auto_corr / ((double)num_use_frames * (double)M);
}


void calc_Xp(
    ReadDump::ExtraReadDump &rd,
    std::vector<std::vector<Eigen::Vector3d>*> &Xp,
    const int &p,
    const int &N,
    const int &M
){
    const double sqrt_2_N = std::sqrt(2./N);
    const int now_frame = std::stoi(rd.ref_private_vars("now_frame"));
    if (Xp[now_frame] != nullptr) return;

    std::vector<Eigen::Vector3d> coods;
    rd.join_3columns(coods, "xu", "yu", "zu");

    Xp[now_frame] = new std::vector<Eigen::Vector3d>(M, Eigen::Vector3d::Zero());
    for (int m = 0; m < M; m++){
        int start_id = m*N;
        for (int n = 0; n < N; n++){
            (*Xp[now_frame])[m] += std::cos(((double)n + 0.5)*(double)p*M_PI/(double)N) * coods[start_id+n];
        }
        (*Xp[now_frame])[m] *= sqrt_2_N;
    }
}

