import numpy as np
from module.common.LAMMPS_potential_table_IO import LAMMPSPotentialTableIO
from module.contents.parameters import PRM


class CalcUip1PolynominalFacade:
    def __call__(self, Uip1_adapter):
        P_target = Uip1_adapter.P_target
        P_CG = Uip1_adapter.P_CG
        Ui = Uip1_adapter.Ui
        if P_CG is None:
            Uip1 = self.__BI(P_target)
        else:
            Uip1 = self.__IBI(Ui, P_target, P_CG)
        Uip1 = self.__std_at_rcut(Uip1_adapter.x_new, Uip1)
        r = np.array(Uip1_adapter.x_new) - PRM.rcut

        # Truncate values above the upper limit
        r_, Uip1_ = self.__truncate_values_above_upper_limit(r, Uip1)

        # fitting
        ok = False
        deg = 30
        while not ok:
            print(f"maximum degree of polynominal function: {deg}")
            z = np.polyfit(r_, Uip1_, deg, rcond=1e-50)
            if (deg % 2 == 0) and (z[0] < 0):
                deg -= 1
                continue
            elif (deg % 2 == 1) and (z[0] > 0):
                deg -= 1
                continue
            else:
                ok = True

        # create LAMMPS table
        fUip1_fitting = np.poly1d(z)
        fdUip1_fitting = np.polyder(fUip1_fitting)
        Uip1_fitting = list(fUip1_fitting(r))
        dUip1_fitting = list(-fdUip1_fitting(r))

        table = LAMMPSPotentialTableIO("", Uip1_adapter.section_name)
        table.x = list(r + PRM.rcut)
        table.E = Uip1_fitting
        table.F = dUip1_fitting
        table.create_table(Min_=Uip1_adapter.Min, Max_=Uip1_adapter.Max)
        return Uip1, Uip1_fitting, z, table.table

    def __IBI(self, Ui, P_target, P_CG) -> list:
        Ui = np.array(Ui)
        P_target = np.array(P_target)
        P_CG = np.array(P_CG)
        Uip1 = Ui + PRM.kBT * np.log(P_CG / P_target)
        return list(Uip1)

    def __BI(self, P) -> list:
        P = np.array(P)
        return list(-PRM.kBT * np.log(P))

    def __std_at_rcut(self, r: list, U: list) -> list:
        idx = np.abs(np.array(r) - PRM.rcut).argmin()
        return list(np.array(U) - U[idx])

    def __truncate_values_above_upper_limit(self, r: list, U: list) -> list:
        r, U = np.array(r), np.array(U)
        flg_array = U < 5.0
        r, U = list(r[flg_array]), list(U[flg_array])
        return r, U
