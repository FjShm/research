import numpy as np
import math
from module.common.LAMMPS_potential_table_IO import LAMMPSPotentialTableIO
from module.contents.parameters import PRM


class Flatten:

    a = PRM.mypoly_a
    x0 = PRM.mypoly_x0

    def tanh(self, x):
        return 0.5 - np.tanh(self.a * (x - self.x0)) * 0.5

    def dtanh(self, x) -> float:
        cosh = np.cosh(self.a * (x - self.x0))
        return -0.5 * self.a / (cosh * cosh)


class CalcUip1MyPolynominalFacade:
    def __call__(self, Uip1_adapter):
        P_target = Uip1_adapter.P_target
        P_CG = Uip1_adapter.P_CG
        Ui = Uip1_adapter.Ui
        ibi_accelerator = Uip1_adapter.ibi_accelerator
        if P_CG is None:
            Uip1 = self.__BI(P_target)
        else:
            Uip1 = self.__IBI(Ui, P_target, P_CG, a=ibi_accelerator)
        Uip1 = self.__std_at_rcut(Uip1_adapter.x_new, Uip1)
        r = np.array(Uip1_adapter.x_new) - PRM.rcut

        # Truncate values above the upper limit
        r_, Uip1_ = self.__truncate_values_above_upper_limit(r, Uip1)

        # fitting
        ok = False
        deg = 30
        Min_ = Uip1_adapter.Min - PRM.rcut
        while not ok:
            print(f"maximum degree of polynominal function: {deg}")
            z = np.polyfit(r_, Uip1_, deg, rcond=1e-50)
            tmp_poly1d = np.poly1d(z)
            tmp_polyder = np.polyder(tmp_poly1d)
            if tmp_poly1d(Min_) < 0 or -tmp_polyder(Min_) < 0:
                print("energy value or force value is incorrect.")
                print(f"energy (at r = {Uip1_adapter.Min}): {tmp_poly1d(Min_)}")
                print(f"force  (at r = {Uip1_adapter.Min}): {tmp_polyder(Min_)}")
                print("fix polynominal dimention as follows.")
                deg -= 1
                continue
            else:
                ok = True

        # r for table
        r_for_table = r
        len_origin = len(r)
        Min_ = Uip1_adapter.Min - PRM.rcut
        Max_ = Uip1_adapter.Max - PRM.rcut
        if Min_ < r_for_table[0]:
            r_for_table = np.array([Min_] + list(r_for_table))
        if r_for_table[-1] < Max_:
            r_for_table = np.array(list(r_for_table) + [Max_])
        if len(r_for_table) != len_origin:
            r_for_table = np.linspace(r_for_table[0], r_for_table[-1], num=len_origin)

        # create LAMMPS table
        fUip1_fitting = np.poly1d(z)
        fdUip1_fitting = np.polyder(fUip1_fitting)
        Uip1_fitting = fUip1_fitting(r)
        Uip1_fitting_for_table = fUip1_fitting(r_for_table)
        dUip1_fitting = fdUip1_fitting(r)
        dUip1_fitting_for_table = fdUip1_fitting(r_for_table)

        # flatten Uip1_fitting, dUip1_fitting
        tanh = Flatten().tanh(r + PRM.rcut)
        dtanh = Flatten().dtanh(r + PRM.rcut)
        dUip1_fitting = -(dUip1_fitting * tanh + Uip1_fitting * dtanh)
        Uip1_fitting = Uip1_fitting * tanh

        tanh = Flatten().tanh(r_for_table + PRM.rcut)
        dtanh = Flatten().dtanh(r_for_table + PRM.rcut)
        dUip1_fitting_for_table = -(
            dUip1_fitting_for_table * tanh + Uip1_fitting_for_table * dtanh
        )
        Uip1_fitting_for_table = Uip1_fitting_for_table * tanh

        # write table
        table = LAMMPSPotentialTableIO("", Uip1_adapter.section_name)
        r_for_table = list(r_for_table + PRM.rcut)
        # 両端の丸め誤差修正
        if math.isclose(Uip1_adapter.Min, r_for_table[0]):
            r_for_table[0] = Uip1_adapter.Min
        if math.isclose(Uip1_adapter.Max, r_for_table[-1]):
            r_for_table[-1] = Uip1_adapter.Max
        table.x = r_for_table
        table.E = Uip1_fitting_for_table
        table.F = dUip1_fitting_for_table
        table.create_table()
        return Uip1, Uip1_fitting, z, table.table

    def __IBI(self, Ui, P_target, P_CG, a=1) -> list:
        Ui = np.array(Ui)
        P_target = np.array(P_target)
        P_CG = np.array(P_CG)
        Uip1 = Ui + PRM.kBT * np.log(P_CG / P_target) * a
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
