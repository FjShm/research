#include "ete_auto_corr.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("ete_auto_corr.log");
    const double dt = param["dt_fs"].as<double>();
    const int max_frames = param["max_frames_to_average"].as<int>(100);
    int N = param["beads_per_chain"].as<int>(-1);
    int M = param["num_chain"].as<int>(-1);

    // -------------------------------
    // max loop
    int max_loop = std::count_rows(ipath);
    boost::progress_display show_progress(max_loop);

    // -------------------------------

    ReadDump::ExtraReadDump rd(ipath);
    rd.read_all_frames();

    int df;
    for (size_t f = 0; f < rd.num_frames; f++){
        rd.header_validation("id", "xu", "yu", "zu");
        int frames_to_average_max = rd.num_frames - f;
        df = frames_to_average_max > max_frames ?
            frames_to_average_max / max_frames : 1;
        // for (frame_t0=0; frame_t0<frames_to_average_max; frame_t0+=df)
        // R(t0)
        // change_now_frame(now_frame, false);
        // R(t0+t)
        ++show_progress;
    }
    R2_n /= (double) count;

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
}

