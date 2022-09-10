#include "crystallinity-C.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string out_crytxt_path = param["output_crytxt_path"].as<std::string>();
    const int k = param["k"].as<int>();
    const int Nbin = param["Nbin"].as<int>();
    const Eigen::Vector3d elongational_vector(0., 0., 1.);

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(dump_path, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    std::ifstream dump{dump_path};
    std::ofstream out_crytxt{out_crytxt_path, std::ios::out | std::ios::trunc};

    // define x_list
    std::vector<double> x_list(Nbin - 1);
    double dx;
    int num = Nbin - 1;
    linspace(num, x_list, dx);

    // output header
    out_crytxt << "TimeStep/x ";
    for (int i = 0; i < Nbin - 1; i++){
        out_crytxt << x_list[i] << " ";
    }
    out_crytxt << std::endl;

    for (int i = 0; i < max_loop; i++){
        // read dump one timestep
        // coordinationsは取り敢えず1要素
        // 後でresize
        std::vector<Eigen::Vector3d> coordinations(1, Eigen::Vector3d(0., 0., 0.));
        std::vector<int> mols;
        int timestep, num_atoms, beads_per_chain, num_chains;
        read_one_timestep_of_dump(
                dump, coordinations, mols, timestep, num_atoms, beads_per_chain);
        num_chains = num_atoms / beads_per_chain;

        // calculate C
        std::vector<double> C((beads_per_chain - 2*k - 1) * num_chains, 100.);
        int looper = 0;
        for (int nc = 0; nc < num_chains; nc++){
            for (int alpha = k; alpha <= beads_per_chain - 2 - k; alpha++){
                int bead_idx = nc*beads_per_chain + alpha;
                Eigen::Vector3d sum(0., 0., 0.);
                Eigen::Vector3d d(0., 0., 0.);
                for (int kk = -k; kk <= k; kk++){
                    d = coordinations[bead_idx + kk + 1] - coordinations[bead_idx + kk];
                    d /= d.norm();
                    sum += d;
                }
                sum /= sum.norm();
                double C_alpha = sum.dot(elongational_vector);
                C[looper] = C_alpha;
                looper++;
            }
        }

        // calculate distribution
        std::vector<int> y_list(Nbin - 1);
        for (int j = 0; j < Nbin - 1; j++){
            double min = -1. + dx * (double)j;
            double max = -1. + dx * ((double)j + 1.);
            if (j == (Nbin - 2)){
                y_list[j] = count_match_value(min, max, C, true, true);
            } else {
                y_list[j] = count_match_value(min, max, C, true, false);
            }
        }

        // output crystallinity
        out_crytxt << timestep << " ";
        for (int j = 0; j < Nbin - 1; j++)
            out_crytxt << y_list[j] << " ";
        out_crytxt << std::endl;

        // update progress bar
        ++show_progress;
    }
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


void read_one_timestep_of_dump(
        std::ifstream &dump,
        std::vector<Eigen::Vector3d> &coordinations,
        std::vector<int> &mols,
        int &timestep,
        int &num_atoms,
        int &beads_per_chain
    ){
    std::string row;

    // TIMESTEP
    std::getline(dump, row);
    std::getline(dump, row);
    timestep = std::stoi(row);

    // NUMBER OF ATOMS
    std::getline(dump, row);
    std::getline(dump, row);
    num_atoms = std::stoi(row);
    coordinations.resize(num_atoms, Eigen::Vector3d(0., 0., 0.));
    mols.resize(num_atoms);

    // BOX BOUNDS xy xz yz pp pp pp
    for (int i = 0; i < 4; i++ )std::getline(dump, row);

    // ATOMS id mol xu yu zu
    std::getline(dump, row);
    int id, mol, mol_max(0);
    double xu, yu, zu;
    for (int i = 0; i < num_atoms; i++){
        dump >> id >> mol >> xu >> yu >> zu;
        mols[id - 1] = mol;
        coordinations[id - 1] << xu, yu, zu;
        if (mol > mol_max) mol_max = mol;
    }
    std::getline(dump, row);
    beads_per_chain = num_atoms / mol_max;
}
 
void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    while(std::getline(in, row))
        if (row == "ITEM: TIMESTEP")
            max_loop++;
}

void linspace(const int& num, std::vector<double>& x_list, double& dx){
    dx = 2. / (double)num;
    x_list[0] = -1. + 0.5 * dx;
    for (int i = 1; i < num; i++)
        x_list[i] = x_list[0] + dx * (double)i;
}

double count_match_value(
        double& min,
        double& max,
        std::vector<double>& C,
        bool include_min=true,
        bool include_max=false
    ){
    int count = 0;
    for (int i = 0; i < C.size(); i++){
        if (include_min && include_max){
            count += (min <= C[i]) && (C[i] <= max);
        } else if (include_min){
            count += (min <= C[i]) && (C[i] <  max);
        } else if (include_max){
            count += (min <  C[i]) && (C[i] <= max);
        } else {
            count += (min <  C[i]) && (C[i] <  max);
        }
    }
    return count;
}
