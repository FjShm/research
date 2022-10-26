#include "scattering.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_scat_path = param["output_scatter_text_path"].as<std::string>();
    const double freq_max = param["freq_max"].as<double>();
    const double freq_min = param["freq_min"].as<double>();
    const int resolution = param["resolution"].as<int>();
    const std::string aspect = param["aspect"].as<std::string>();
    const std::string read_frame = param["frame"]["type"].as<std::string>();
    std::vector<int> read_timestep;
    std::vector<double> read_ratio;
    // read_frame: "all", "timestep", "ratio"
    if (read_frame == "timestep"){
        read_timestep = param["frame"]["content"].as< std::vector<int> >();
    } else if (read_frame == "ratio"){
        read_ratio = param["frame"]["content"].as< std::vector<double> >();
    }

    const std::vector<double> freq_axis = linspace(freq_min, freq_max, resolution);
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
    int max_loop = 0;
    count_number_of_rows(dump_path, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    ReadDump::ReadDump rd(dump_path);
    std::ifstream rotxt{rot_path};
    std::ofstream out_scat{out_scat_path, std::ios::out | std::ios::trunc};
    std::string rotxt_row;
    if (!rotxt){
        std::cout << rot_path << " is not exist." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // skip header of rotation.txt
    for (int i = 0; i < 2; i++) std::getline(rotxt, rotxt_row);

    // pre-read dump
    if (read_frame != "all"){
        rd.read_all_frames();
        rd.set_want_frames(read_ratio, read_timestep);
    }

    int timestep;
    while(rd.read_1frame()){
        std::getline(rotxt, rotxt_row);
        Eigen::Matrix3d rot;
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);

        // check timestep
        if (timestep != rd.timestep){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep << std::endl
                << "dump file:    " << rd.timestep << std::endl;
            std::exit(EXIT_FAILURE);
        }
        out_scat << rd.timestep;

        // rotate coordinations
        std::vector<Eigen::Vector3d> coordinations;
        rd.join_3columns(coordinations, "x", "y", "z");
        std::vector<Eigen::Vector3d> rotated_coordinations(rd.num_atoms);
        for (int i = 0; i < rd.num_atoms; i++){
            rotated_coordinations[i] = coordinations[i].transpose() * rot;
        }

        // move center of gravity of system to origin
        Eigen::Vector3d cog = {0., 0., 0.};
        for (int i = 0; i < rd.num_atoms; i++){
            cog += rotated_coordinations[i];
        }
        cog /= rd.num_atoms;
        for (int i = 0; i < rd.num_atoms; i++){
            rotated_coordinations[i] -= cog;
        }

        // scatter
        double kr;
        Eigen::MatrixXd Sks(resolution, resolution);
        #pragma omp parallel for
        for (int i = 0; i <= int(resolution*0.5); i++){
            for (int j = 0; j < resolution; j++){
                double sum_coskr(0.), sum_sinkr(0.);
                //#pragma omp parallel for reduction(+:sum_coskr, sum_sinkr)
                for (int n = 0; n < rd.num_atoms; n++){
                    kr = freq_axis[j]*rotated_coordinations[n](ax0)
                        + freq_axis[resolution-i-1]*rotated_coordinations[n](ax1);
                    sum_coskr += std::cos(kr);
                    sum_sinkr += std::sin(kr);
                }
                //double Sk = (sum_coskr*sum_coskr + sum_sinkr*sum_sinkr)/(double)rd.num_atoms;
                //out_scat << " " << Sk;
                Sks(i, j) = (sum_coskr*sum_coskr + sum_sinkr*sum_sinkr)/(double)rd.num_atoms;
            }
        }

        // 点対称のため
        for (int i = 0; i < int(resolution*0.5); i++)
            for (int j = 0; j < resolution; j++)
                Sks(resolution-i-1, resolution-j-1) = Sks(i, j);

        for (int i = 0; i < resolution; i++)
            for (int j = 0; j < resolution; j++)
                out_scat << " " << Sks(i, j);
        out_scat << std::endl;

        // update progress bar
        ++show_progress;
    }
}


std::vector<double> linspace(const double &min, const double &max, int num){
    std::vector<double> ret(num);
    double dx = (max - min) / ((double)num - 1.);
    for (size_t i = 0; i < num; i++) ret[i] = min + dx*(double)i;
    return ret;
}

void dumpcell_to_vector_converter(
        double &xlo_b, double &xhi_b, double &xy,
        double &ylo_b, double &yhi_b, double &xz,
        double &zlo_b, double &zhi_b, double &yz,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin
        ){
    double xlo,xhi, ylo,yhi, zlo,zhi;
    xlo = xlo_b - std::min({0., xy, xz, xy + xz});
    xhi = xhi_b - std::max({0., xy, xz, xy + xz});
    ylo = ylo_b - std::min({0., yz});
    yhi = yhi_b - std::max({0., yz});
    zlo = zlo_b;
    zhi = zhi_b;
    cell_origin << xlo, ylo, zlo;
    a << xhi - xlo, 0., 0.;
    b << xy, yhi - ylo, 0.;
    c << xz, yz, zhi - zlo;
}


std::vector<std::string> split(const std::string &s, char delim){
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) {
            elems.push_back(item);
        }
    }
    return elems;
}


void rotationtxt2rotmatrix(std::string &row, Eigen::Matrix3d &rot, int &timestep){
    std::vector<std::string> row_split = split(row, ' ');
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


void vector_to_dumpcell_converter(
        Eigen::Matrix3d &dumpcell,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin
        ){
    double xlo,xhi,xy, ylo,yhi,xz, zlo,zhi,yz;
    xlo = cell_origin(0);
    ylo = cell_origin(1);
    zlo = cell_origin(2);
    xy = b(0);
    xz = c(0);
    yz = c(1);
    xhi = a(0) + xlo;
    yhi = b(1) + ylo;
    zhi = c(2) + zlo;

    double xlo_b,xhi_b, ylo_b,yhi_b, zlo_b,zhi_b;
    xlo_b = xlo + std::min({0., xy, xz, xy + xz});
    xhi_b = xhi + std::max({0., xy, xz, xy + xz});
    ylo_b = ylo + std::min({0., yz});
    yhi_b = yhi + std::max({0., yz});
    zlo_b = zlo;
    zhi_b = zhi;

    dumpcell << xlo_b, xhi_b, xy,
             ylo_b, yhi_b, xz,
             zlo_b, zhi_b, yz;
}


void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    while(std::getline(in, row)){
        if ("ITEM: TIMESTEP" == row)
            max_loop++;
    }
}

