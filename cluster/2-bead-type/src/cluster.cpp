#include "cluster.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dumppath_par = param["dumppath_particles_make_up_cluster"].as<std::string>();
    const std::string dumppath_core = param["dumppath_cluster_core"].as<std::string>();
    const std::string outpath = param["output_path"].as<std::string>("LOG.out");
    const int smax = param["smax"].as<int>(12);
    const double encutoff = param["encount_cutoff"].as<double>();
    const int total_frames = param["total_frames"].as<int>(1001);

    // -------------------------------
    // max loop
    boost::progress_display show_progress(total_frames);

    // -------------------------------
    ReadDump::ExtraReadDump rd_par(dumppath_par, false); // turn off sort by id
    ReadDump::ExtraReadDump rd_core(dumppath_core, false);
    std::ofstream out{outpath, std::ios::out | std::ios::trunc};

    // output header
    out << "TimeStep";
    for (int i = 1; i <= smax; i++) out << " " << i;
    out << "<" << std::endl;
    while(rd_par.read_1frame() && rd_core.read_1frame()){
        if (rd_par.timestep != rd_core.timestep){
            std::cerr << "TimeStep of '" << dumppath_par << "' and '" << dumppath_core
                << "' does not match.\n"
                << dumppath_par << ": " << rd_par.timestep << "\n"
                << dumppath_core << ": " << rd_core.timestep << "\n";
            std::exit(EXIT_FAILURE);
        }

        // set each coordinations
        std::vector<Eigen::Vector3d> coordinations_par, coordinations_core;
        rd_par.join_3columns(coordinations_par, "x", "y", "z");
        rd_core.join_3columns(coordinations_core, "x", "y", "z");
        const int Np = rd_par.num_atoms;
        const int Nc = rd_core.num_atoms;

        // kdtree
        kdt::KDTree<Point> kdtree_par = create_kdtree(coordinations_par, rd_par, Np);
        kdt::KDTree<Point> kdtree_core = create_kdtree(coordinations_core, rd_core, Nc);

        // count cluster
        std::vector<int> n(smax+2, 0);
        std::vector<bool> isCounted_par(Np, false);
        std::vector<bool> isCounted_core(Nc, false);
        for (int i = 0; i < Nc; i++){
            int s = measure_cluster_size(
                i,
                isCounted_par,
                isCounted_core,
                kdtree_par,
                kdtree_core,
                encutoff,
                coordinations_par,
                coordinations_core,
                Np,
                Nc
            );
            s = std::min(s, smax+1);
            n[s]++;
        }

        // output
        out << rd_core.timestep;
        for (int s = 1; s <= smax; s++) out << " " << n[s];
        out << std::endl;

        // progress
        ++show_progress;
    }
}

int measure_cluster_size(
    int i,
    std::vector<bool> &isCounted_par,
    std::vector<bool> &isCounted_core,
    kdt::KDTree<Point> &kdtree_par,
    kdt::KDTree<Point> &kdtree_core,
    const double &encutoff,
    std::vector<Eigen::Vector3d> &coordinations_par,
    std::vector<Eigen::Vector3d> &coordinations_core,
    const int &Np,
    const int &Nc,
    int cluster_size
){
    // 1. core周りのparticlesを探索
    if (isCounted_core[i%Nc] == true) return cluster_size;
    isCounted_core[i%Nc] = true;

    Point query = Point(coordinations_core[i]);
    std::vector<int> idxes_par = kdtree_par.radiusSearch(query, encutoff);

    // 2. 見つけた全てのparticlesについて, 他のcoreとの会合チェック
    for (int j : idxes_par){
        if (isCounted_par[j%Np] == true) continue;
        isCounted_par[j%Np] = true;
        cluster_size++; // 初出のparticleならclusterとしてカウント

        Point _query = Point(coordinations_par[j]);
        std::vector<int> idxes_core = kdtree_core.radiusSearch(_query, encutoff);
        for (int k: idxes_core){
            cluster_size = measure_cluster_size(
                k,
                isCounted_par,
                isCounted_core,
                kdtree_par,
                kdtree_core,
                encutoff,
                coordinations_par,
                coordinations_core,
                Np,
                Nc,
                cluster_size
            );
        }
    }
    return cluster_size;
}


kdt::KDTree<Point> create_kdtree(
    std::vector<Eigen::Vector3d> &coordinations,
    ReadDump::ExtraReadDump & rd,
    const int &N
){
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
    return kdtree;
}
