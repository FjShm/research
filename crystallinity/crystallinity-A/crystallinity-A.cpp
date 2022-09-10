#include "crystallinity-A.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string out_crytxt_path = param["output_crytxt_path"].as<std::string>();
    const int k = param["k"].as<int>();

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(dump_path, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    std::ifstream dump{dump_path};
    std::ofstream out_crytxt{out_crytxt_path, std::ios::out | std::ios::trunc};

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

        // calculate lambda
        std::vector<double> lambda(beads_per_chain - 2*k - 1, 0.);
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
                double norm = sum.norm();
                lambda[alpha - k] += norm / (2.*(double)k + 1.);
            }
        }

        // average lambda
        for (int alpha = k; alpha <= beads_per_chain - 2 - k; alpha++){
            lambda[alpha - k] /= (double)num_chains;
        }

        // output crystallinity
        out_crytxt << timestep << " ";
        for (int alpha = k; alpha < beads_per_chain - 2 - k; alpha++){
            out_crytxt << lambda[alpha - k] << " ";
        }
        out_crytxt << lambda[beads_per_chain - 2 - k - k] << std::endl;


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
