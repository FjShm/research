import numpy as np
from module.common.LAMMPS_potential_table_IO import LAMMPSPotentialTableIO
from module.contents.parameters import PRM
from module.contents.functions.b_spline_scipy import b_spline_scipy


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

        # Array length tuning
        len_sp = len_origin = len(x)
        Min = Uip1_adapter.Min
        Max = Uip1_adapter.Max

        # smooth Uip1 by B-spline
        x_sp = x
        Uip1_sp = b_spline_scipy(x, Uip1, num=100)
        F_sp = list(-np.array(self.__deriv(x, Uip1_sp)))

        # extrapolation
        if x[-1] < Max:
            x_sp, Uip1_sp, F_sp = self.__extrapolate(x_sp, Uip1_sp, F_sp, Max=Max)
            print("extrapolate Max")
        if Min < x[0]:
            x_sp, Uip1_sp, F_sp = self.__extrapolate(x_sp, Uip1_sp, F_sp, Min=Min)
            print("extrapolate Min")

        # spaceing
        x_sp, Uip1_sp = b_spline_scipy(x_sp, Uip1_sp, num=500, spacing=True)
        F_sp = list(-np.array(self.__deriv(x_sp, Uip1_sp)))
        F_sp[-1] = 0

        # create LAMMPS table
        table = LAMMPSPotentialTableIO("", Uip1_adapter.section_name)
        table.x = x_sp
        table.E = Uip1_sp
        table.F = F_sp
        table.create_table()
        return Uip1, (x_sp, Uip1_sp), None, table.table

    def __IBI(self, Ui, P_target, P_CG) -> list:
        Ui = np.array(Ui)
        P_target = np.array(P_target)
        P_CG = np.array(P_CG)
        Uip1 = Ui + PRM.kBT * np.log(P_CG / P_target)
        return list(Uip1)

    def __BI(self, P) -> list:
        P = np.array(P)
        return list(-PRM.kBT * np.log(P))

    def __deriv(self, x, y) -> list:
        d = []
        for i in range(len(x)):
            if i == 0:
                d_ = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            elif i == len(x) - 1:
                d_ = (y[-1] - y[-2]) / (x[-1] - x[-2])
            else:
                d_ = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
            d.append(d_)
        return d

    def __extrapolate(
        self, x: list, y: list, d: list, Min: float = None, Max: float = None,
    ) -> tuple:
        if Min is None:
            # d[-1] = y[-1] / (Max - x[-1])
            d[-1] = 0
            x[-1] = Max
            y[-1] = 0
        else:
            # 最大傾きを用いて線形外挿
            # dmax = d[0]
            dmax = 3000
            for d_ in d[1:]:
                if dmax < d_:
                    dmax = d_
            d[0] = dmax
            y[0] = y[0] + dmax * (x[0] - Min)
            x[0] = Min
        return x, y, d
