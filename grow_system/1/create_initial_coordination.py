import numpy as np
import random

"""
AB粗視化ビーズモデルの初期配置作成プログラム
A-B-A-B-A-.......

ABA = 90°
BAB = 180°

   \ 
    A
     \ 
      B
     /
    A
   /
  B
   \ 
    A
     \ 
      B
     /
    A --- start


Equilibrium bond length
A-B : 0.3142 = AB
B-A : 0.2412 = BA

chain end-to-end distance : (N-1)/2 * AA = ete_distance
cell edge length : (2*ete_distance + AB + BA) * M = cell_edge_length

basic bond vector : (AB1, BA1, AB2, BA2)
A-B 1 : (1/sqrt(2), 1/sqrt(2), 0) * AB
B-A 1 : (-1/sqrt(2), 1/sqrt(2), 0) * BA
A-B 2 : (-1/sqrt(2), 1/sqrt(2), 0) * AB
B-A 2 : (1/sqrt(2), 1/sqrt(2), 0) * BA

chain starting position


random rotation matrix : Rn(theta)
Rn(theta) r = rcos(theta) + n(n \dot r)(1-cos(theta)) + (n \cross r)sin(theta)
r : 任意ベクトル
n : 任意の回転軸ベクトル, n = || np.array(random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)) ||
theta : 回転角, theta = random.uniform(-np.pi, np.pi)
"""


def Cells(cell_edge_x, cell_edge_y, cell_edge_z):
    # f.write(f"{-cell_edge_length*0.5} {cell_edge_length*0.5} xlo xhi\n")
    # f.write(f"{-cell_edge_length*0.5} {cell_edge_length*0.5} ylo yhi\n")
    # f.write(f"{-cell_edge_length*0.5} {cell_edge_length*0.5} zlo zhi\n")
    list_f = f.readlines()
    row = 1
    for line in list_f:
        row += 1
        if "dihedral types" in line:
            break
    cell_edge_x, cell_edge_y, cell_edge_z = (
        np.array(cell_edge_x),
        np.array(cell_edge_y),
        np.array(cell_edge_z),
    )
    list_f.insert(row, "\n\n")
    list_f.insert(row, f"{cell_edge_z.min()-1} {cell_edge_z.max()+1} zlo zhi\n")
    list_f.insert(row, f"{cell_edge_y.min()-1} {cell_edge_y.max()+1} ylo yhi\n")
    list_f.insert(row, f"{cell_edge_x.min()-1} {cell_edge_x.max()+1} xlo xhi\n")

    return list_f


def Masses():
    # (atom type) (質量)
    f.write("Masses\n\n")
    f.write(f"1 {mass[0]}\n2 {mass[1]}\n")
    # for m in range(M):
    #    for n in range(1,N+1):
    #        f.write(f"{m*N + n} {mass[n%2]}\n")
    f.write("\n\n")


def Atoms(cell_edge_x, cell_edge_y, cell_edge_z):
    # (通し番号) (polymer番号) (atom type) (x) (y) (z)
    basic_bond_vector = (
        np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0]) * AB,
        np.array([-1 / np.sqrt(2), 1 / np.sqrt(2), 0]) * BA,
        np.array([-1 / np.sqrt(2), 1 / np.sqrt(2), 0]) * AB,
        np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 0]) * BA,
    )
    f.write("Atoms\n\n")
    for m in range(1, M + 1):
        # chainの向きをランダムに決める
        rnd_theta = random.uniform(-np.pi, np.pi)
        n_vec = np.array(
            [random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)]
        )
        n_vec /= np.linalg.norm(n_vec)
        bond_vector = [
            bbv * np.cos(rnd_theta)
            + n_vec * np.inner(n_vec, bbv) * (1 - np.cos(rnd_theta))
            + np.cross(n_vec, bbv) * np.sin(rnd_theta)
            for bbv in basic_bond_vector
        ]
        # chainのスタート位置(A末端の座標)をランダムに決める
        rnd_x = random.uniform(-cell_edge_length * 0.5, cell_edge_length * 0.5)
        rnd_y = random.uniform(-cell_edge_length * 0.5, cell_edge_length * 0.5)
        rnd_z = random.uniform(-cell_edge_length * 0.5, cell_edge_length * 0.5)
        chain_starting_point = np.array([rnd_x, rnd_y, rnd_z])
        pos = chain_starting_point
        for n in range(1, N + 1):
            f.write(
                f"{(m-1)*N + n:d} {m} {(n-1)%2+1:d} {pos[0]:.5f} {pos[1]:.5f} {pos[2]:.5f}\n"
            )
            cell_edge_x.append(pos[0])
            cell_edge_y.append(pos[1])
            cell_edge_z.append(pos[2])
            pos += bond_vector[(n - 1) % 4]
    f.write("\n\n")
    return cell_edge_x, cell_edge_y, cell_edge_z


def cell_edge(xhi, xlo, yhi, ylo, zhi, zlo, pos):
    if pos[0] < xlo:
        xlo = pos[0]
    elif pos[0] > xhi:
        xhi = pos[0]
    if pos[1] < ylo:
        ylo = pos[1]
    elif pos[1] > yhi:
        yhi = pos[1]
    if pos[2] < zlo:
        zlo = pos[2]
    elif pos[2] > zhi:
        zhi = pos[2]
    return xhi, xlo, yhi, ylo, zhi, zlo


def Bonds():
    # (通し番号) (bond type) (I) (J)
    f.write("Bonds\n\n")
    for m in range(1, M + 1):
        for n in range(1, N):
            f.write(
                f"{(m-1)*(N-1) + n:d} {(n-1)%2+1} {(m-1)*N + n:d} {(m-1)*N + n+1:d}\n"
            )  # bond数は N-1 [個/chain]
    f.write("\n\n")


def Angles():
    # (通し番号) (angle type) (I) (J) (K)
    f.write("Angles\n\n")
    for m in range(1, M + 1):
        for n in range(1, N - 1):
            f.write(
                f"{(m-1)*(N-2) + n:d} {(n-1)%2+1} {(m-1)*N + n:d} {(m-1)*N + n+1:d} {(m-1)*N + n+2:d}\n"
            )  # angle数は N-2 [個/chain]
    f.write("\n\n")


def Dihedrals():
    # (通し番号) (dihedral type) (I) (J) (K) (L)
    f.write("Dihedrals\n\n")
    for m in range(1, M + 1):
        for n in range(1, N - 2):
            f.write(
                f"{(m-1)*(N-3) + n:d} {(n-1)%2+1} {(m-1)*N + n:d} {(m-1)*N + n+1:d} {(m-1)*N + n+2:d} {(m-1)*N + n+3:d}\n"
            )  # dihedral数は N-3 [個/chain]
    f.write("\n\n")


if __name__ == "__main__":
    N = 49
    M = 512
    AB = 2.432
    BA = 3.108  # angstroms
    mass = (28.05316, 40.06386)

    bond_num = (N - 1) * M
    angle_num = (N - 2) * M
    dihedral_num = (N - 3) * M
    ete_distance = np.sqrt(AB * AB + BA * BA) * (N - 1) / 2
    cell_edge_length = 10 * ete_distance + AB + BA
    print(f"N = {N}, M = {M}\nAB = {AB} [nm]\nBA = {BA} [nm]")
    print(f"end-to-end distance = {ete_distance:.4f} [nm]")
    print(f"cell edge length    = {cell_edge_length} [nm]")

    cell_edge_x, cell_edge_y, cell_edge_z = [], [], []
    with open("Init.data", mode="w") as f:
        f.write("LAMMPS data file\n\n")
        f.write(f"{N*M} atoms\n")
        f.write(f"{bond_num} bonds\n")
        f.write(f"{angle_num} angles\n")
        f.write(f"{dihedral_num} dihedrals\n\n")

        f.write(f"{2} atom types\n")
        f.write("2 bond types\n")
        f.write("2 angle types\n")
        f.write("2 dihedral types\n\n")
        Masses()
        cell_edge_x, cell_edge_y, cell_edge_z = Atoms(
            cell_edge_x, cell_edge_y, cell_edge_z
        )
        Bonds()
        Angles()
        Dihedrals()
    with open("Init.data", mode="r+") as f:
        list_f = Cells(cell_edge_x, cell_edge_y, cell_edge_z)
    with open("Init.data", mode="w") as f:
        f.writelines(list_f)

    # subprocess.run(["qsub", "laurel.sh"])
