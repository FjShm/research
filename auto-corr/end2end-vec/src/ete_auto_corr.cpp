#include "ete_auto_corr.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("ete_auto_corr.log");
    const double dt = param["dt_fs"].as<double>();
    const int ave_frames = param["max_frames_to_average"].as<int>(100);
    const int total_frames_dump = param["total_frames_dump"].as<int>(-1);
    int N = param["beads_per_chain"].as<int>();
    int M = param["num_chain"].as<int>();
    constexpr double fs2ns = 1.e-6;

    ReadDump::ExtraReadDump rd(ipath);
    rd.read_all_frames(total_frames_dump);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(rd.num_frames);
    // -------------------------------

    std::vector<double> time(rd.num_frames), auto_corr(rd.num_frames);
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
            calc_auto_corr(rd, fp, df, ave_frames_, N, M);
        ++show_progress;
    }

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (int i = 0; i < rd.num_frames; i++)
        out << time[i] << " " << auto_corr[i]
            << " " << auto_corr[i]/auto_corr[0] << std::endl;
}


double calc_auto_corr(
    ReadDump::ExtraReadDump &rd, int fp, const double &df, int num_use_frames, int N, int M
){
    double auto_corr = 0.;
    int counter = 0;
    for (int i = 0; i < num_use_frames; i++){
        int now_frame = std::round(i*df);
        std::vector<Eigen::Vector3d> coord_t, coord_0;
        rd.jump_frames(now_frame, true);
        rd.join_3columns(coord_0, "xu", "yu", "zu");
        rd.jump_frames(fp);
        rd.join_3columns(coord_t, "xu", "yu", "zu");
        for (int m = 0; m < M; m++){
            counter++;
            int start_id = m*N;
            int end_id = start_id + N - 1;
            auto_corr +=
                (coord_t[end_id]-coord_t[start_id]).dot(coord_0[end_id]-coord_0[start_id]);
        }
    }
    return auto_corr / ((double)num_use_frames * (double)M);
}

