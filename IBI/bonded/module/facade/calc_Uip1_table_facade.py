import numpy as np
from module.common.LAMMPS_potential_table_IO import LAMMPSPotentialTableIO
from module.contents.parameters import PRM


class CalcUip1TableFacade:
    def __call__(self, Uip1_adapter):
        P_target = Uip1_adapter.P_target
        P_CG = Uip1_adapter.P_CG
        Ui = Uip1_adapter.Ui
        if P_CG is None:
            Uip1 = self.__BI(P_target)
        else:
            Uip1 = self.__IBI(Ui, P_target, P_CG)

        x = np.array(Uip1_adapter.x_new)

        # create LAMMPS table
        F = self.__deriv(x, Uip1)

        table = LAMMPSPotentialTableIO("", Uip1_adapter.section_name)
        table.x = list(x)
        table.E = Uip1
        table.F = F
        return Uip1, None, None, table.table

    def __IBI(self, Ui, P_target, P_CG):
        Ui = np.array(Ui)
        P_target = np.array(P_target)
        P_CG = np.array(P_CG)
        Uip1 = Ui + PRM.kBT * np.log(P_CG / P_target)
        return list(Uip1)

    def __BI(self, P):
        P = np.array(P)
        return list(-PRM.kBT * np.log(P))

    def __deriv(self, x, y):
        d = []
        for i in range(len(x) - 1):
            if i == 0:
                d_ = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            elif i == len(x) - 2:
                d_ = (y[-1] - y[-2]) / (x[-1] - x[-2])
            else:
                d_ = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
            d.append(d_)
        return d
