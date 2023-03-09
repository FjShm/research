import numpy as np
from scipy.interpolate import interp1d
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

        # smooth Uip1 by B-spline
        x_sp = x
        Uip1_sp = b_spline_scipy(x, Uip1, num=100)
        F_sp = list(-np.array(self.__deriv(x, Uip1_sp)))

        # extrapolation
        x_sp, Uip1_sp, F_sp = self.__extrapolate(
            x_sp,
            Uip1_sp,
            F_sp,
            Max=Uip1_adapter.Max,  # xの定義域の最大
            Min=Uip1_adapter.Min,  # xの定義域の最小
            extrapolate_type=Uip1_adapter.extrapolate_type,  # 外挿の関数型
            Maxd_coeff=Uip1_adapter.Maxd_coeff,  # (外挿時に使用)x=Maxにおけるdに掛ける値
            Mind_coeff=Uip1_adapter.Mind_coeff,  # (外挿時に使用)x=Minにおけるdに掛ける値
            Max_parabola_axis=Uip1_adapter.Max_parabola_axis,  # (harmonic外挿時に使用)Max外挿時の放物線の軸
            Min_parabola_axis=Uip1_adapter.Min_parabola_axis,  # (harmonic外挿時に使用)Min外挿時の放物線の軸
        )

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
        extrapolate_type: str = "linear",  # linear, harmonic
        Mind_coeff: float = 1.0,
        Maxd_coeff: float = 1.0,
        Min_parabola_axis: float = None,
        Max_parabola_axis: float = None,
    ) -> tuple:
        # dはtiltと符号が逆であることに注意して外挿
        # linear: y = ax + b
        # harmonic: y = p(x - a)^2 + q
        x_rhs, x_lhs = [], []
        y_rhs, y_lhs = [], []
        lenx = len(x)
        if Max is not None and x[-1] < Max:
            # 右側へ外挿
            print("extrapolate Max")
            x_rhs = np.linspace(x[-1], Max, num=len(x) * 10)[1:]
            if extrapolate_type == "linear":
                a = -d[-1] * Maxd_coeff
                b = y[-1] - a * x[-1]
                y_rhs = a * x_rhs + b
            elif extrapolate_type == "harmonic":
                m = -d[-1] * Maxd_coeff
                a = (
                    np.mean([x[0], x[-1]])
                    if Max_parabola_axis is None
                    else Max_parabola_axis
                )
                p = 0.5 * m / (x[-1] - a)
                q = y[-1] - 0.5 * m * (x[-1] - a)
                y_rhs = p * (x_rhs - a) ** 2 + q
            else:
                exit(f"Invalid 'extrapolate_type': '{extrapolate_type}'")
        elif Min is not None and Min < x[0]:
            # 左側へ外挿
            print("extrapolate Min")
            x_lhs = np.linspace(Min, x[0], num=len(x) * 10)[:-1]
            if extrapolate_type == "linear":
                a = -d[0] * Maxd_coeff
                b = y[0] - a * x[0]
                y_lhs = a * x_lhs + b
            elif extrapolate_type == "harmonic":
                m = -d[0] * Maxd_coeff
                a = (
                    np.mean([x[0], x[-1]])
                    if Max_parabola_axis is None
                    else Max_parabola_axis
                )
                p = 0.5 * m / (x[0] - a)
                q = y[0] - 0.5 * m * (x[0] - a)
                y_lhs = p * (x_lhs - a) ** 2 + q
            else:
                exit(f"Invalid 'extrapolate_type': '{extrapolate_type}'")

        x = list(x_lhs) + list(x) + list(x_rhs)
        y = list(y_lhs) + list(y) + list(y_rhs)

        # spacing
        f = interp1d(x, y)
        x = np.linspace(x[0], x[-1], num=lenx)
        y = f(x)
        x = list(x)
        d = list(-np.array(self.__deriv(x, y)))
        return x, y, d
