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
        Min = Uip1_adapter.Min  # xの定義域の最小
        Miny_fixed = Uip1_adapter.Miny_fixed  # (外挿時に使用)x=Minにおけるyの値
        Mind_fixed = Uip1_adapter.Mind_fixed  # (外挿時に使用)x=Minにおけるdの値
        Max = Uip1_adapter.Max  # xの定義域の最大
        Maxy_fixed = Uip1_adapter.Maxy_fixed  # (外挿時に使用)x=Maxにおけるyの値
        Maxd_fixed = Uip1_adapter.Maxd_fixed  # (外挿時に使用)x=Maxにおけるdの値

        # smooth Uip1 by B-spline
        x_sp = x
        Uip1_sp = b_spline_scipy(x, Uip1, num=100)
        F_sp = list(-np.array(self.__deriv(x, Uip1_sp)))

        # extrapolation
        if x[-1] < Max:
            x_sp, Uip1_sp, F_sp = self.__extrapolate(
                x_sp,
                Uip1_sp,
                F_sp,
                Max=Max,
                Maxy_fixed=Maxy_fixed,
                Maxd_fixed=Maxd_fixed,
            )
            print("extrapolate Max")
        if Min < x[0]:
            x_sp, Uip1_sp, F_sp = self.__extrapolate(
                x_sp,
                Uip1_sp,
                F_sp,
                Min=Min,
                Miny_fixed=Miny_fixed,
                Mind_fixed=Mind_fixed,
            )
            print("extrapolate Min")

        # spacing
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
        self,
        x: list,
        y: list,
        d: list,
        Min: float = None,
        Max: float = None,
        Miny_fixed: float = None,
        Maxy_fixed: float = None,
        Mind_fixed: float = None,
        Maxd_fixed: float = None,
        spacing: 
    ) -> tuple:
        # dはtiltと符号が逆であることに注意して外挿
        if Min is None:
            # 右側へ外挿
            if Maxd_fixed is None:
                d_ex = d[-1]
            else:
                d_ex = Maxd_fixed
            d[-1] = d_ex
            if Maxy_fixed is None:
                y[-1] = y[-1] - d_ex * (Max - x[-1])
            else:
                y[-1] = Maxy_fixed
            x[-1] = Max
        else:
            # 左側へ外挿
            if Mind_fixed is None:
                d_ex = d[0]
            else:
                d_ex = Mind_fixed
            d[0] = d_ex
            if Miny_fixed is None:
                y[0] = y[0] + d_ex * (x[0] - Min)
            else:
                y[0] = Miny_fixed
            x[0] = Min
        return x, y, d
