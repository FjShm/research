#include "cluster.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string outpath = param["output_path"].as<std::string>("LOG.out");
    const int smax = param["smax"].as<int>(12);
    const double encutoff = param["encount_cutoff"].as<double>();
    const int total_frames = param["total_frames"].as<int>(1001);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(total_frames);

    // -------------------------------
    ReadDump::ExtraReadDump rd(dump_path, false); // turn off sort by id
    std::ofstream out{outpath, std::ios::out | std::ios::trunc};

    // output header
    out << "TimeStep";
    for (int i = 1; i <= smax; i++) out << " " << i;
    out << "<" << std::endl;
    while(rd.read_1frame()){
        std::vector<Eigen::Vector3d> coordinations;
        rd.join_3columns(coordinations, "x", "y", "z");
        const int N = rd.num_atoms;

        // 周期境界条件も考えてKD-tree作成
        const Eigen::Vector3d zeros = Eigen::Vector3d::Zero();
        const std::vector<Eigen::Vector3d> _a_ = {zeros, -rd.cellbox_a, rd.cellbox_a};
        const std::vector<Eigen::Vector3d> _b_ = {zeros, -rd.cellbox_b, rd.cellbox_b};
        const std::vector<Eigen::Vector3d> _c_ = {zeros, -rd.cellbox_c, rd.cellbox_c};

        std::vector<Point> points(27*N);
        coordinations.resize(N*27);
        int counter = 0;
        int mirror = 0;
        for (int xi = 0; xi < 3; xi++)
            for (int yi = 0; yi < 3; yi++)
                for (int zi = 0; zi < 3; zi++){
                    for (int i = 0; i < N; i++){
                        Eigen::Vector3d tmp = coordinations[i] + _a_[xi] + _b_[yi] + _c_[zi];
                        points[counter++] = Point(tmp);
                        coordinations[i + N*mirror] = tmp;
                    }
                    mirror++;
                }
        kdt::KDTree<Point> kdtree(points);

        // cluster counterrrr
        std::vector<int> n(smax+2, 0);
        std::vector<bool> isCounted(N, false);
        int cluster_size = 0;
        for (int i = 0; i < N; i++){
            int s = measure_cluster_size(i, isCounted, kdtree, encutoff, coordinations, N);
            s = std::min(s, smax+1);
            n[s]++;
        }

        // output
        out << rd.timestep;
        for (int s = 1; s <= smax; s++) out << " " << n[s];
        out << std::endl;

        // check
        int sum = 0;
        for (int s = 0; s <= smax; s++) sum += s*n[s];
        if (sum != N){
            std::cerr << "Warning: The sum of s*n_s is NOT equal to N\n"
                << "'smax' may be too small.\n"
                << "TimeStep: " << rd.timestep << std::endl;
        }

        // progress
        ++show_progress;
    }
}

int measure_cluster_size(
    int i,
    std::vector<bool> &isCounted,
    kdt::KDTree<Point> &kdtree,
    const double &encutoff,
    std::vector<Eigen::Vector3d> &coordinations,
    const int &N,
    int cluster_size
){
    if (isCounted[i%N] == true) return cluster_size;
    isCounted[i%N] = true;
    cluster_size++;

    Point query = Point(coordinations[i]);
    std::vector<int> idxes = kdtree.radiusSearch(query, encutoff);

    for (int j : idxes){
        cluster_size = measure_cluster_size(
            j,
            isCounted,
            kdtree,
            encutoff,
            coordinations,
            N,
            cluster_size
        );
    }
    return cluster_size;
}
