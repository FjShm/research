import numpy as np
from scipy.optimize import curve_fit
from module.contents.potential_functions import PotentialFunctions
from module.common.LAMMPS_potential_table_IO import LAMMPSPotentialTableIO
from module.contents.parameters import PRM


class CalcUip1FunctionFacade:
    def __call__(self, Uip1_adapter):
        P_target = Uip1_adapter.P_target
        P_CG = Uip1_adapter.P_CG
        Ui = Uip1_adapter.Ui
        if P_CG is None:
            Uip1 = self.__BI(P_target, shift=Uip1_adapter.shift)
        else:
            Uip1 = self.__IBI(Ui, P_target, P_CG)

        # fitting
        x = np.array(Uip1_adapter.x_new)
        func = PotentialFunctions()(Uip1_adapter.function_type)
        popt, pcov = curve_fit(
            func, x, np.array(Uip1), bounds=Uip1_adapter.bounds, maxfev=10000
        )
        Uip1_fitting = list(func(x, *popt))

        # create LAMMPS table
        dfunc = PotentialFunctions()(Uip1_adapter.function_type, d=True)
        F = list(-dfunc(x, *popt))

        table = LAMMPSPotentialTableIO("", Uip1_adapter.section_name)
        table.x = list(x)
        table.E = Uip1_fitting
        table.F = F
        table.create_table(Min_=Uip1_adapter.Min, Max_=Uip1_adapter.Max)
        return Uip1, Uip1_fitting, popt, table.table

    def __IBI(self, Ui, P_target, P_CG) -> list:
        Ui = np.array(Ui)
        P_target = np.array(P_target)
        P_CG = np.array(P_CG)
        Uip1 = Ui + PRM.kBT * np.log(P_CG / P_target)
        return list(Uip1)

    def __BI(self, P, shift=False) -> list:
        P = np.array(P)
        Uip1 = -PRM.kBT * np.log(P)
        if shift is True:
            Uip1 -= Uip1.min()
        return list(Uip1)
