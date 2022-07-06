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

    def create_table(self):
        if len(self.x) == 0:
            print("There is no data to create as LAMMPS potential Table.")
        else:
            Min = min(self.x)
            Max = max(self.x)
            N = len(self.x)
            table = "# LAMMPS POTENTIAL TABLE\n\n"
            table += f"{self.section_name}\n"
            table += f"N {N} R {Min} {Max}\n"
            for i in range(N):
                n = i + 1
                table += f"{n} {self.x[i]} {self.E[i]} {self.F[i]}\n"
            self.table = table

    def write(self):
        with open(self.path, mode="w") as f:
            f.write(self.table)
