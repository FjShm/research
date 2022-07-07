import numpy as np
from typing import Union


class LAMMPSPotentialTableIO:

    x: list
    E: list
    F: list
    path: str
    section_name: str
    table: str

    def __init__(self, path: str, section_name: str):
        self.path = path
        self.section_name = section_name
        self.x = []
        self.E = []
        self.F = []

    def read(self):
        with open(self.path, mode="r") as f:
            ok = False
            for line in f:
                if ok:
                    line_ = [float(val) for val in line.split(" ")]
                    self.x.append(line_[1])
                    self.E.append(line_[2])
                    self.F.append(line_[3])
                if line == "":
                    continue
                elif line[0] == "#":
                    continue
                if line == self.section_name:
                    ok = True

    def create_table(
        self, Min_: Union[float, None] = None, Max_: Union[float, None] = None
    ):
        if len(self.x) == 0:
            print("There is no data to create as LAMMPS potential Table.")
        else:
            x, E, F = self.x, self.E, self.F
            Min = min(self.x)
            Max = max(self.x)
            # linear extrapolation
            if (Min_ is not None) and Min_ < Min:
                E0 = F[0] * (x[0] - Min_) + E[0]
                x, E, F = [Min_] + x, [E0] + E, [F[0]] + F
                Min = Min_
            if (Max_ is not None) and Max_ > Max:
                En = F[-1] * (x[-1] - Max_) + E[-1]
                x, E, F = x + [Max_], E + [En], F + [F[-1]]
                Max = Max_
            N = len(x)
            table = "# LAMMPS POTENTIAL TABLE\n\n"
            table += f"{self.section_name}\n"
            table += f"N {N} R {Min} {Max}\n"
            for i in range(N):
                n = i + 1
                table += f"{n} {x[i]} {E[i]} {F[i]}\n"
            self.table = table

    def write(self):
        with open(self.path, mode="w") as f:
            f.write(self.table)
