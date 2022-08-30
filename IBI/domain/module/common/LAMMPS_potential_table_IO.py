import os
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
                line = line.split()
                if len(line) == 0:
                    continue
                elif line[0] in ("", "#", "\n", "N"):
                    continue
                elif ok:
                    line_ = [float(val) for val in line]
                    if line_[2] > 10000:
                        continue
                    self.x.append(line_[1])
                    self.E.append(line_[2])
                    self.F.append(line_[3])
                elif line[0] == self.section_name:
                    ok = True
                else:
                    print(f"error format of table {self.path}")
                    print(f"{line}")
                    exit()
        return self

    def create_table(self, Min_=None, Max_=None):
        if len(self.x) == 0:
            print("There is no data to create as LAMMPS potential Table.")
        else:
            x, E, F = self.x, self.E, self.F
            Min = min(self.x)
            Max = max(self.x)
            N = len(x)
            table = "# LAMMPS POTENTIAL TABLE\n\n"
            table += f"{self.section_name}\n"
            table += f"N {N} R {Min} {Max}\n\n"
            for i in range(N):
                n = i + 1
                table += f"{n} {x[i]} {E[i]} {F[i]}\n"
            self.table = table

    def write(self):
        with open(self.path, mode="w") as f:
            f.write(self.table)
        param_dir = os.path.dirname(self.path)
        param_path = os.path.join(param_dir, f"table_{self.section_name}.param")
        table_abspath = os.path.abspath(self.path)
        with open(param_path, mode="w") as f:
            f.write(f"{table_abspath} {self.section_name}")
