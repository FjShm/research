#include "scattering.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_dir = param["output_dir"].as<std::string>(".");
    const std::string aspect = param["aspect"].as<std::string>("xz");
    const int frames = param["dump_frames"].as<int>(-1);
    const std::vector<double> kx_all = param["kx"].as< std::vector<double> >();
    const std::vector<double> ky_all = param["ky"].as< std::vector<double> >();
    const std::vector<double> ratio = param["ratio"].as< std::vector<double> >();
    int N = param["N"].as<int>(-1);
    int M = param["M"].as<int>(-1);

    int ax0, ax1;
    if ((aspect == "xy") || (aspect == "yx")){
        ax0 = 0;
        ax1 = 1;
    } else if ((aspect == "yz") || (aspect == "zy")){
        ax0 = 1;
        ax1 = 2;
    } else if ((aspect == "zx") || (aspect == "xz")){
        ax0 = 0;
        ax1 = 2;
    } else {
        std::cout << "input parameter error:\n"
            << "aspect: " << aspect << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // -------------------------------
    // max loop
    if (frames == -1){
        std::cout << "It is recommended to put 'dump_frames' in the input file.\n"
            << "This will automatically start counting the number of frames,\n"
            << "but this process puts a heavy load on the I/O of the system.\n";
    }
    int max_loop = frames == -1 ?
        std::count_rows(dump_path, "TIMESTEP"):
        frames;
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    ReadDump::ExtraReadDump rd(dump_path);
    std::ifstream rotxt{rot_path};

    std::string rotxt_row;
    bool f1rst_loop = true;

    // skip header of rotation.txt
    for (int i = 0; i < 2; i++) std::getline(rotxt, rotxt_row);

    int loop_count = 0;
    int ratio_idx = 0;
    while(rd.read_1frame()){
        std::getline(rotxt, rotxt_row);
        // read rotation.txt
        int timestep;
        Eigen::Matrix3d rot;
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);

        // check timestep
        if (timestep != rd.timestep){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep << std::endl
                << "dump file:    " << rd.timestep << std::endl;
            std::exit(EXIT_FAILURE);
        }

        double d_ratio = std::abs((double)loop_count++ / (double)frames - ratio[ratio_idx]);
        double d_ratio_next = std::abs((double)loop_count / (double)frames - ratio[ratio_idx]);
        if (d_ratio_next < d_ratio && loop_count < frames){
            ++show_progress;
            continue;
        }

        // validation
        rd.header_validation("id", "x", "y", "z");
        rd.add_column_if_not_exist("mol", N, M);
        int mol = rd.header_map->at("mol");

        // set wrapped_coordination, N, M
        std::vector<Eigen::Vector3d> w_coordinations;
        rd.join_3columns(w_coordinations, "x", "y", "z");
        M = M != -1 ? M : rd.max_of_col("mol");
        N = N != -1 ? N : rd.num_atoms/M;

        // rotate coordinations
        for (int i = 0; i < rd.num_atoms; i++)
            w_coordinations[i] = w_coordinations[i].transpose() * rot;

        // calculate scattering function
        std::string path = out_dir + "/LOG." + std::to_string(rd.timestep) + ".out";
        std::ofstream out{path, std::ios::out | std::ios::trunc};
        double Sk, kr;
        for (double ky : ky_all){
            for (double kx : kx_all){
                double sum_coskr(0.0), sum_sinkr(0.0);
                for (int n = 0; n < rd.num_atoms; n++){
                    kr = kx*w_coordinations[n](ax0) + ky*w_coordinations[n](ax1);
                    sum_coskr += std::cos(kr);
                    sum_sinkr += std::sin(kr);
                }
                Sk = (sum_coskr*sum_coskr + sum_sinkr*sum_sinkr)/(double)rd.num_atoms;
                out << Sk << " ";
            }
            out << "\n";
        }
        out.close();

        // update progress bar
        ++show_progress;
        ratio_idx++;
    }
}


void rotationtxt2rotmatrix(std::string &row, Eigen::Matrix3d &rot, int &timestep){
    std::vector<std::string> row_split = std::split(row, ' ');
    timestep = std::stoi(row_split[0]);
    rot(0, 0) = std::stod(row_split[1]);
    rot(0, 1) = std::stod(row_split[2]);
    rot(0, 2) = std::stod(row_split[3]);
    rot(1, 0) = std::stod(row_split[4]);
    rot(1, 1) = std::stod(row_split[5]);
    rot(1, 2) = std::stod(row_split[6]);
    rot(2, 0) = std::stod(row_split[7]);
    rot(2, 1) = std::stod(row_split[8]);
    rot(2, 2) = std::stod(row_split[9]);
}

