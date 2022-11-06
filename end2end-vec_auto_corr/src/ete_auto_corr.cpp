#include "ete_auto_corr.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("ete_auto_corr.log");
    const double dt = param["dt_fs"].as<double>();
    const int max_frames = param["max_frames_to_average"].as<int>(100);
    int N = param["beads_per_chain"].as<int>();
    int M = param["num_chain"].as<int>();
    constexpr double fs2ns = 1.e-6;


    ReadDump::ExtraReadDump rd(ipath);
    rd.read_all_frames();

    // -------------------------------
    // max loop
    boost::progress_display show_progress(rd.num_frames);
    // -------------------------------

    int dframes;
    std::vector<double> time(rd.num_frames), auto_corr(rd.num_frames);
    for (size_t f = 0; f < rd.num_frames; f++){
        // calc time
        rd.jump_frames(0, true);
        int timestep0 = rd.timestep;
        rd.jump_frames(f);
        time[f] = ((double)rd.timestep - timestep0)*dt*fs2ns;

        // calc auto corr
        rd.header_validation("id", "xu", "yu", "zu");
        int frames_to_average_max = rd.num_frames - f;
        dframes = frames_to_average_max > max_frames ?
            frames_to_average_max / max_frames : 1;
        auto_corr[f] =
            calc_auto_corr(rd, f, dframes, frames_to_average_max, N, M);
        ++show_progress;
    }

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (int i = 0; i < rd.num_frames; i++)
        out << time[i] << " " << auto_corr[i] << std::endl;
}


double calc_auto_corr(
    ReadDump::ExtraReadDump &rd, int tframes, int dframes, int num_use_frames, int N, int M
){
    double auto_corr = 0.;
    for (int i = 0; i < num_use_frames; i++){
        int now_frame = i*dframes;
        std::vector<Eigen::Vector3d> coord_t, coord_0;
        rd.jump_frames(now_frame, true);
        rd.join_3columns(coord_0, "xu", "yu", "zu");
        rd.jump_frames(tframes);
        rd.join_3columns(coord_t, "xu", "yu", "zu");
        for (int m = 0; m < M; m++){
            int start_id = m*N;
            int end_id = start_id + N - 1;
            auto_corr +=
                (coord_t[end_id]-coord_t[start_id]).dot(coord_0[end_id]-coord_0[start_id]);
        }
    }
    return auto_corr / ((double)num_use_frames * (double)M);
}
