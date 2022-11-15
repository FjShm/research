#include "msd.h"


int main(int argc, char* argv[]){
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string ipath = param["input_dump_path"].as<std::string>();
    const std::string opath = param["output_path"].as<std::string>("MSD.txt");
    const double dt = param["dt_fs"].as<double>();
    const int total_frames_dump = param["total_frames_dump"].as<int>(-1);
    int N = param["beads_per_chain"].as<int>(-1);
    int M = param["num_chain"].as<int>(-1);
    constexpr double fs2ns = 1.e-6;

    ReadDump::ExtraReadDump rd(ipath);
    if (total_frames_dump != -1)
        rd.read_all_frames(total_frames_dump);
    
    std::ofstream out{opath, std::ios::out | std::ios::trunc};

    // -------------------------------
    // max loop
    boost::progress_display show_progress(rd.num_frames);
    // -------------------------------

    bool f1rst_loop = true;
    std::vector<Eigen::Vector3d> initial_positions;
    while (rd.read_1frame()){
        rd.header_validation("id", "xu", "yu", "zu");
        rd.add_column_if_not_exist("mol", N, M);
        std::vector<Eigen::Vector3d> coord;
        rd.join_3columns(coord, "xu", "yu", "zu");
        int mol = rd.header_map->at("mol");
        M = (int)rd.max_of_col("mol");

        std::vector<Eigen::Vector3d> positions(M, Eigen::Vector3d::Zero());
        std::vector<int> N_all(M, 0);
        for (size_t id_ = 0; id_ < rd.num_atoms; id_++){
            int m = (int)rd.atoms_all_data->coeff(id_, mol) - 1;
            positions[m] += coord[id_];
            N_all[m]++;
            if (id_+1 >= rd.num_atoms || m != rd.atoms_all_data->coeff(id_+1, mol))
                positions[m] /= (double)N_all[m];
        }

        if (f1rst_loop){
            initial_positions = positions;
            out << "time ";
            for (int m = 1; m <= M; m++)
                out << m << " ";
            out << "average\n";
            f1rst_loop = false;
        }

        // calc MSD
        out << (double)rd.timestep * fs2ns * dt << " ";
        double sum = 0.;
        for (int m = 0; m < M; m++){
            double distance = (positions[m] - initial_positions[m]).squaredNorm();
            sum += distance;
            out << distance << " ";
        }
        out << sum/(double)M << std::endl;

        ++show_progress;
    }
}

