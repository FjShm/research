#include "calc_R2_n_from_dump.h"


void compute_R2_n(std::ifstream&, int, int, int, std::vector<double>&);


int main(){
    /***************************************
    input file
    ----------------------------------------
        (input filename)
        (output filename)
        (N)
        (M)
    ----------------------------------------
    ***************************************/

    // N,Mを入力から得る
    // 1<=n<=49の配列を生成
    // header
    // セル情報取得
    // あるtimestepにおける全座標をunwrapped状態で取得
    //// chainごとに配列で座標を持つ
    // chainの座標配列をforループで回す
    //// chain内で, 2<=n<=49でforループ, R2を算出
    // R2をnで割る
    // 次のtimestepへ.

    std::string line, ipath, opath;
    int N, M, NM, count(0), timestep, beads_total;
    std::cin >> ipath >> opath >> N >> M;
    NM = N * M;
    std::vector<double> R2_n(N+1, 0);
    std::ifstream in{ipath};

    while(std::getline(in, line)){
        if (line == "ITEM: TIMESTEP"){
            in >> timestep;
            std::cout << "Timestep: " << timestep << std::endl;
        } else if (line == "ITEM: NUMBER OF ATOMS"){
            in >> beads_total;
            if (beads_total != NM){
                std::cout << "N, M is wrong. N*M = " << NM
                << ", total number of beads = " << beads_total << std::endl;
                return -1;
            }
        } else if (line == "ITEM: BOX BOUNDS xy xz yz pp pp pp"){
            continue;
        } else if (line == "ITEM: ATOMS id xu yu zu"){
            compute_R2_n(in, N, M, NM, R2_n);
            count++;
        } else continue;
    }
    R2_n /= (double) count;

    // output
    std::ofstream out{opath, std::ios::out | std::ios::trunc};
    for (int n = 2; n <= N; n++){
        out << n-1 << ' ' << R2_n[n] << std::endl;
    }
}


void compute_R2_n(std::ifstream& in, int N, int M, int NM, std::vector<double>& R2_n){
    std::vector<double> R2_n_tmp(N+1, 0);
    std::vector< std::vector<double> > posx(M, std::vector<double> (N));
    std::vector< std::vector<double> > posy(M, std::vector<double> (N));
    std::vector< std::vector<double> > posz(M, std::vector<double> (N));
    int index, id, mol;
    for (int i = 0; i < NM; i++){
        in >> index;
        mol = (index - 1) / N;
        id  = (index - 1) - mol*N;
        in >> posx[mol][id] >> posy[mol][id] >> posz[mol][id];
    }

    for (int n = 2; n <= N; n++){
        for (int m = 0; m < M; m++){
            for (int start_id = 0; start_id <= N-n; start_id++){
                int end_id = start_id + n - 1;
                double dx = (posx[m][start_id] - posx[m][end_id]);
                double dy = (posy[m][start_id] - posy[m][end_id]);
                double dz = (posz[m][start_id] - posz[m][end_id]);
                R2_n_tmp[n] += dx*dx + dy*dy + dz*dz;
            }
        }
        R2_n_tmp[n] /= (M * (N - n + 1)) * (n-1); // M chains, N-n+1 pairs for each chains, n-1 is the number of bonds
    }
    R2_n += R2_n_tmp;
}
