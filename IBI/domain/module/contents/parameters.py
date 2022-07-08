class PRM:
    kB = 1.380649e-23  # [J/K]
    NA = 6.02214076e23  # [/mol]
    T = 360  # [K]
    cal2J = 4.184  # 1 cal = 4.184 J
    J2cal = 1.0 / cal2J
    kBT = (kB * J2cal * 1e-3) * T * NA  # [kcal/mol]
    rcut = 20  # [Ang]
