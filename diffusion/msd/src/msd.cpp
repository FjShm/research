#include "msd.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("ete_auto_corr.log");
    const double dt = param["dt_fs"].as<double>();
    const int total_frames_dump = param["total_frames_dump"].as<int>(-1);
    int N = param["beads_per_chain"].as<int>(-1);
    int M = param["num_chain"].as<int>(-1);
    constexpr double fs2ns = 1.e-6;

    ReadDump::ExtraReadDump rd(ipath);
    if (total_frames_dump != -1)
        rd.read_all_frames(total_frames_dump);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(rd.num_frames);
    // -------------------------------

    while (rd.read_1frame()){
        std::vector<>
        ++show_progress;
    }

    // output                                                                       :q
}

