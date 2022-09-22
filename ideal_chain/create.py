import numpy as np
import argparse
import random
from tqdm import tqdm


class Cell:
    xlo = -5.9220768353553296e01
    xhi = 5.9220768353553296e01
    xy = 0.0
    ylo = -5.9220768353553296e01
    yhi = 5.9220768353553296e01
    xz = 0.0
    zlo = -5.9220768353553296e01
    zhi = 5.9220768353553296e01
    yz = 0.0

    def __init__(self):
        self.row1 = (self.xlo, self.xhi, self.xy)
        self.row2 = (self.ylo, self.yhi, self.xz)
        self.row3 = (self.zlo, self.zhi, self.yz)
        self.cell = (self.row1, self.row2, self.row3)
        self.cell_txt = "ITEM: BOX BOUNDS xy xz yz pp pp pp\n"
        for row in self.cell:
            self.cell_txt += f"{row[0]} {row[1]} {row[2]}\n"


def Args():
    parser = argparse.ArgumentParser(
        description="Generate ideal chains as appropriate."
    )
    parser.add_argument("-b", "--bondlength", default=1.0, type=float)
    parser.add_argument("-N", "--num_atoms", default=49, type=int)
    parser.add_argument("-M", "--num_chains", default=512, type=int)
    parser.add_argument("-s", "--steps", default=1, type=int)
    return parser.parse_args()


if __name__ == "__main__":
    args = Args()
    length = args.bondlength
    N = args.num_atoms
    M = args.num_chains
    NM = N * M
    steps = args.steps

    cell = Cell()
    with open("ideal.lammpstrj", mode="w") as f:
        for step in tqdm(range(steps + 1)):
            f.write(f"ITEM: TIMESTEP\n{step}\n")
            f.write(f"ITEM: NUMBER OF ATOMS\n{NM}\n")
            f.write(cell.cell_txt)
            f.write("ITEM: ATOMS id mol xu yu zu\n")
            for m in range(M):
                n_ = m * N + 1
                m_ = m + 1
                posx = random.uniform(cell.xlo, cell.xhi)
                posy = random.uniform(cell.ylo, cell.yhi)
                posz = random.uniform(cell.zlo, cell.zhi)
                pos = np.array([posx, posy, posz])
                f.write(f"{n_} {m_} {pos[0]} {pos[1]} {pos[2]}\n")
                for n in range(N - 1):
                    n_ += 1
                    bondx = random.uniform(-1, 1)
                    bondy = random.uniform(-1, 1)
                    bondz = random.uniform(-1, 1)
                    bond = np.array([bondx, bondy, bondz])
                    bond *= length / np.linalg.norm(bond)
                    pos += bond
                    f.write(f"{n_} {m_} {pos[0]} {pos[1]} {pos[2]}\n")
