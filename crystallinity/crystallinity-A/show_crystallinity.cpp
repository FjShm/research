#include "show_crystallinity.h"


int main(int argc, char* argv[]){
    // input
    YAML::Node param = YAML::LoadFile(argv[1]);
    const std::string dump_path = param["input_dump_path"].as<std::string>();
    const std::string rot_path = param["input_rotationtxt_path"].as<std::string>();
    const std::string out_dump_path = param["output_cry_dump_path"].as<std::string>();
    const int k = param["k"].as<int>();
    const int beta = param["beta"].as<int>();
    const double lath = param["lambda_threshold"].as<double>();
    const double dith = param["neighbor_stem_distance_threshold"].as<double>();
    const double thth = param["theta_deg_threshold"].as<double>();
    const double sqrdLath = lath*lath;
    const double sqrdDith = dith*dith;
    const double cos_thth = std::cos(thth * M_PI/180.);

    // -------------------------------
    // max loop
    int max_loop = 0;
    count_number_of_rows(rot_path, max_loop);
    boost::progress_display show_progress(max_loop);

    // -------------------------------
    std::ifstream dump{dump_path};
    std::ifstream rotxt{rot_path};
    std::ofstream out_dump{out_dump_path, std::ios::out | std::ios::trunc};

    std::string rotxt_row;;

    // skip header of rotation.txt
    for (int i = 0; i < 2; i++) std::getline(rotxt, rotxt_row);


    while(std::getline(rotxt, rotxt_row)){
        // read rotation.txt
        int timestep;
        Eigen::Matrix3d rot;
        rotationtxt2rotmatrix(rotxt_row, rot, timestep);


        // read dump one timestep
        // coordinations, molsは
        // read_one_timestep_of_dump内でresize
        Eigen::Vector3d a, b, c, cell_origin;
        std::vector<Eigen::Vector3d> coordinations;
        std::vector<int> mols;
        int timestep_dump, num_atoms;
        read_one_timestep_of_dump(
                dump, a, b, c, cell_origin, coordinations, mols, timestep_dump, num_atoms);

        int mol_max = *std::max_element(mols.begin(), mols.end());
        int bead_per_chain = num_atoms / mol_max;

        // check timestep
        if (timestep != timestep_dump){
            std::cout << "The timestep in the dump file and rotation.txt do not match.\n"
                << "rotation.txt: " << timestep
                << "dump file:    " << timestep_dump << std::endl;
            return 1;
        }

        // calc center of gravity of each polymers
        std::vector<Eigen::Vector3d> cogs(mol_max);
        Eigen::Vector3d sum_tmp;
        sum_tmp << 0., 0., 0.;
        for (int i = 0; i < num_atoms; i++){
            sum_tmp += coordinations[i];
            if (i == num_atoms - 1 || mols[i] < mols[i+1]){
                cogs[mols[i]-1] = sum_tmp / (double)bead_per_chain;
                sum_tmp << 0., 0., 0.;
            }
        }

        // move vector
        std::vector<Eigen::Vector3d> movecs(mol_max);
        for (int i = 0; i < mol_max; i++){
            cogs[i] -= cell_origin;
            int xi, yi, zi;
            zi = std::floor(cogs[i](2) / c(2));
            yi = std::floor((cogs[i](1) - (double)zi*c(1)) / b(1));
            xi = std::floor((cogs[i](0) - (double)yi*b(0) - (double)zi*c(0)) / a(0));
            movecs[i] = (double)xi*a + (double)yi*b + (double)zi*c;
        }

        // move polymer into simulation box
        std::vector<Eigen::Vector3d> position_fixed_coordinations(num_atoms);
        for (int i = 0; i < num_atoms; i++){
            position_fixed_coordinations[i] = coordinations[i] - movecs[mols[i]-1];
        }

        // rotate coordinations
        std::vector<Eigen::Vector3d> output_coordinations(num_atoms);
        for (int i = 0; i < num_atoms; i++){
            output_coordinations[i] = position_fixed_coordinations[i].transpose() * rot;
        }

        // store alpha (bead index)
        int num_stems = mol_max*(bead_per_chain - 2*k - 2*beta);
        int idx = 0;
        std::vector<int> alphas(num_stems);
        for (int i = 0; i < mol_max; i++)
            for (int alpha = k+beta; alpha < bead_per_chain-k-beta; alpha++)
                alphas[idx++] = i*bead_per_chain + alpha;

        // abc座標系
        std::vector<Eigen::Vector3d> abc_coords(num_stems);
        for (int i = 0; i < num_stems; i++){
            abc_coords[i] = position_fixed_coordinations[alphas[i]] - cell_origin;
            double xd, yd, zd;
            zd =  abc_coords[i](2) / c(2);
            yd = (abc_coords[i](1) - (double)zd*c(1)) / b(1);
            xd = (abc_coords[i](0) - (double)yd*b(0) - (double)zd*c(0)) / a(0);
            abc_coords[i] << xd, yd, zd;
        }


        // store lamda
        std::vector<double> lambdas(num_stems);
        std::vector<Eigen::Vector3d> stem_vecs(num_stems);
        std::vector<double> stem_vec_norms(num_stems);
        for (int i = 0; i < num_stems; i++){
            Eigen::Vector3d sum_d;
            sum_d << 0., 0., 0.;
            for (int kk = -k; kk <= k; kk++){
                Eigen::Vector3d d;
                d = position_fixed_coordinations[alphas[i]+kk+beta]
                    - position_fixed_coordinations[alphas[i]+kk-beta];
                d /= d.norm();
                sum_d += d;
            }
            double sum_d_norm = sum_d.norm();
            lambdas[i] = sum_d_norm / (2.*(double)k + 1.);
            stem_vecs[i] = sum_d;
            stem_vec_norms[i] = sum_d_norm;
        }

        // calculate cos
        Eigen::MatrixXd all_cos(num_stems, num_stems);
        Eigen::Vector3d stem_i, stem_j;
        for (int i = 0; i < num_stems; i++){
            stem_i = position_fixed_coordinations[alphas[i]];
            for (int j = i+1; j < num_stems; j++){
                stem_j = position_fixed_coordinations[alphas[j]];
                all_cos(i, j)
                    = std::abs(stem_i.dot(stem_j)) / (stem_vec_norms[i] * stem_vec_norms[j]);
                all_cos(j, i) = all_cos(i, j);
            }
        }

        // abc座標系での壁沿い判定
        double wa, wb, wc;
        Eigen::Vector3d ey, ez;
        ey << 0, 1, 0;
        ez << 0, 0, 1;
        wa = dith / a(0);
        wb = dith / b.dot(ey);
        wc = dith / c.dot(ez);

        // 周期境界条件も考える
        std::vector<int> total_neighbor_stems(num_stems, 0);
        std::vector<int> cry_neighbor_stems(num_stems, 0);
        Eigen::Vector3d a1, a2;
        Eigen::Vector3d b1, b2;
        double sqrdDistance;
        for (int i = 0; i < num_stems; i++){
            stem_i = position_fixed_coordinations[alphas[i]];
            a1 = stem_i - stem_vecs[i]*0.5;
            a2 = stem_i + stem_vecs[i]*0.5;
            for (int j = 0; j < num_stems; j++){
                if (i == j) continue;
                // stem-jは周期境界を考慮
                bool done = false;
                int xd_fl, yd_fl, zd_fl;
                int _xi, xi_, _yi, yi_, _zi, zi_;
                xd_fl = std::floor(abc_coords[j](0));
                yd_fl = std::floor(abc_coords[j](1));
                zd_fl = std::floor(abc_coords[j](2));
                _xi = xi_ = xd_fl;
                _yi = yi_ = yd_fl;
                _zi = zi_ = zd_fl;
                if (abc_coords[j](0) - (double)xd_fl < wa) _xi = xd_fl-1;
                if (abc_coords[j](1) - (double)yd_fl < wb) _yi = yd_fl-1;
                if (abc_coords[j](2) - (double)zd_fl < wc) _zi = zd_fl-1;
                if ((double)xd_fl + 1. - abc_coords[j](0) < wa) xi_ = xd_fl+1;
                if ((double)yd_fl + 1. - abc_coords[j](1) < wb) yi_ = yd_fl+1;
                if ((double)zd_fl + 1. - abc_coords[j](2) < wc) zi_ = zd_fl+1;
                for (int xi = _xi; xi <= xi_; xi++){
                    for (int yi = _yi; yi <= yi_; yi++){
                        for (int zi = _zi; zi <= zi_; zi++){
                            stem_j = position_fixed_coordinations[alphas[j]]
                                + (double)xi * a
                                + (double)yi * b
                                + (double)zi * c;
                            Eigen::Vector3d tmp = stem_i - stem_j;
                            if (std::abs(tmp(0)) > dith ||
                                    std::abs(tmp(1)) > dith ||
                                    std::abs(tmp(2)) > dith) continue;
                            b1 = stem_j - stem_vecs[j]*0.5;
                            b2 = stem_j + stem_vecs[j]*0.5;
                            sqrdDistance_btw_linesegment_linesegment(a1, a2, b1, b2, sqrdDistance);
                            if (sqrdDistance <= sqrdDith){
                                done = true;
                                total_neighbor_stems[i]++;
                                if (all_cos(i, j) >= cos_thth){
                                    cry_neighbor_stems[i]++;
                                }
                                break;
                            }
                        }
                        if (done) break;
                    }
                    if (done) break;
                }
            }
        }

        // output new dump
        write_to_newdump(
                out_dump, timestep, num_atoms, a, b, c,
                cell_origin, output_coordinations, mols,
                total_neighbor_stems, cry_neighbor_stems, alphas
                );


        // update progress bar
        ++show_progress;
    }
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

void read_one_timestep_of_dump(
        std::ifstream &dump,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin,
        std::vector<Eigen::Vector3d> &coordinations,
        std::vector<int> &mols,
        int &timestep_dump,
        int &num_atoms
        ){
    std::string row;

    // TIMESTEP
    std::getline(dump, row);
    std::getline(dump, row);
    timestep_dump = std::stoi(row);

    // NUMBER OF ATOMS
    std::getline(dump, row);
    std::getline(dump, row);
    num_atoms = std::stoi(row);
    coordinations.resize(num_atoms);
    mols.resize(num_atoms);

    // BOX BOUNDS xy xz yz pp pp pp
    double xlo_b,xhi_b,xy, ylo_b,yhi_b,xz, zlo_b,zhi_b,yz;
    std::getline(dump, row);
    dump >> xlo_b >> xhi_b >> xy >> ylo_b >> yhi_b >> xz >> zlo_b >> zhi_b >> yz;
    dumpcell_to_vector_converter(
            xlo_b, xhi_b, xy, ylo_b, yhi_b, xz, zlo_b, zhi_b, yz, a, b, c, cell_origin);
    std::getline(dump, row);

    // ATOMS id mol xu yu zu
    std::getline(dump, row);
    int id, mol;
    double xu, yu, zu;
    for (int i = 0; i < num_atoms; i++){
        dump >> id >> mol >> xu >> yu >> zu;
        mols[id - 1] = mol;
        coordinations[id - 1] << xu, yu, zu;
    }
    std::getline(dump, row);
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

void write_to_newdump(
        std::ofstream &out,
        int &timestep,
        int& num_atoms,
        Eigen::Vector3d &a,
        Eigen::Vector3d &b,
        Eigen::Vector3d &c,
        Eigen::Vector3d &cell_origin,
        std::vector<Eigen::Vector3d> &coordinations,
        std::vector<int> &mols,
        std::vector<int> &total_neighbor_stems,
        std::vector<int> &cry_neighbor_stems,
        std::vector<int> &alphas
        ){
    out << "ITEM: TIMESTEP\n" << timestep << std::endl;
    out << "ITEM: NUMBER OF ATOMS\n" << num_atoms << std::endl;
    out << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";

    // calculate triclinic box
    Eigen::Matrix3d dumpcell;
    vector_to_dumpcell_converter(dumpcell, a, b, c, cell_origin);
    out << dumpcell(0, 0) << " " << dumpcell(0, 1) << " " << dumpcell(0, 2) << std::endl;
    out << dumpcell(1, 0) << " " << dumpcell(1, 1) << " " << dumpcell(1, 2) << std::endl;
    out << dumpcell(2, 0) << " " << dumpcell(2, 1) << " " << dumpcell(2, 2) << std::endl;

    out << "ITEM: ATOMS id mol xu yu zu cry_stems total_neighbor_stems\n";
    int al = 0;
    for (int i = 0; i < num_atoms; i++){
        out << i + 1 << " " << mols[i] << " " << coordinations[i](0) << " "
            << coordinations[i](1) << " " << coordinations[i](2) << " ";
        if (i == alphas[al]){
            out << cry_neighbor_stems[al] << " " << total_neighbor_stems[al] << std::endl;
            al++;
        } else {
            out << 0 << " " << 0 << std::endl;
        }
    }
}


void count_number_of_rows(const std::string &path, int &max_loop){
    std::ifstream in{path};
    std::string row;
    for (int i = 0; i < 2; i++) std::getline(in, row);
    while(std::getline(in, row)) max_loop++;
}


// calc squared distance between two line-segments

void sqrdDistance_btw_line_line(
        Eigen::Vector3d &a1,
        Eigen::Vector3d &a2,
        Eigen::Vector3d &b1,
        Eigen::Vector3d &b2,
        double &tA,
        double &tB,
        double &sqrdDistance
        ){
    Eigen::Vector3d a21 = a2 - a1;
    Eigen::Vector3d b21 = b2 - b1;
    double a, b, c, d, e, f;
    a = a21.squaredNorm();
    b = -2. * a21.dot(b21);
    c = b21.squaredNorm();
    d = 2.*(a1.dot(a21) - a21.dot(b1));
    e = 2.*(b1.dot(b21) - a1.dot(b21));
    f = a1.squaredNorm() + b1.squaredNorm() - 2. * a1.dot(b1);

    tB = -(2.*a*e - b*d) / (4.*a*c - b*b);
    tA = -(b*tB + d) / (2.*a);
    sqrdDistance = f - (d*d)/(4.*a) - (2.*a*e - b*d)*(2.*a*e - b*d)/(4.*a*(4.*a*c - b*b));
}


void sqrdDistance_btw_line_point(
        Eigen::Vector3d &a1,
        Eigen::Vector3d &b1,
        Eigen::Vector3d &b2,
        double &t,
        double &sqrdDistance
        ){
    Eigen::Vector3d b21 = b2 - b1;
    t = (a1.dot(b21) - b1.dot(b21)) / (b21.squaredNorm());
    sqrdDistance = (a1 - b1 - t * b21).squaredNorm();
}


void sqrdDistance_btw_point_point(
        Eigen::Vector3d &a1,
        Eigen::Vector3d &b1,
        double &sqrdDistance
        ){
    sqrdDistance = (a1 - b1).squaredNorm();
}


void sqrdDistance_btw_linesegment_linesegment(
        Eigen::Vector3d &a1,
        Eigen::Vector3d &a2,
        Eigen::Vector3d &b1,
        Eigen::Vector3d &b2,
        double &sqrdDistance
        ){
    double tA, tB;

    // line line
    sqrdDistance_btw_line_line(a1, a2, b1, b2, tA, tB, sqrdDistance);
    if (0. <= tA && tA <= 1. && 0. <= tB && tB <= 1.)
        return;

    // line point
    double t, min;
    std::vector<double> sqrdDistance_min(4, DBL_MAX);
    bool ok = false;
    sqrdDistance_btw_line_point(a1, b1, b2, t, min);
    sqrdDistance_min[0] = min;
    if (0. <= t && t <= 1.){
        ok = true;
    }
    sqrdDistance_btw_line_point(a2, b1, b2, t, min);
    sqrdDistance_min[1] = min;
    if (0. <= t && t <= 1.){
        ok = true;
    }
    sqrdDistance_btw_line_point(b1, a1, a2, t, min);
    sqrdDistance_min[2] = min;
    if (0. <= t && t <= 1.){
        ok = true;
    }
    sqrdDistance_btw_line_point(b2, a1, a2, t, min);
    sqrdDistance_min[3] = min;
    if (0. <= t && t <= 1.){
        ok = true;
    }
    if (ok){
        sqrdDistance = *min_element(sqrdDistance_min.begin(), sqrdDistance_min.end());
        return;
    }

    // point point
    for (int i = 0; i < 4; i++) sqrdDistance_min[i] = DBL_MAX;
    sqrdDistance_btw_point_point(a1, b1, min);
    sqrdDistance_min[0] = min;
    sqrdDistance_btw_point_point(a1, b2, min);
    sqrdDistance_min[1] = min;
    sqrdDistance_btw_point_point(a2, b1, min);
    sqrdDistance_min[2] = min;
    sqrdDistance_btw_point_point(a2, b2, min);
    sqrdDistance_min[3] = min;
    sqrdDistance = *min_element(sqrdDistance_min.begin(), sqrdDistance_min.end());
    return;
}
