from module.common.LAMMPS_potential_table_IO import LAMMPSPotentialTableIO


class GetUiTableFacade:
    def __call__(self, Ui_adapter):
        previous_params_filepath = Ui_adapter.previous_params_filepath
        with open(previous_params_filepath, mode="r") as f:
            for line in f:
                coeff = line.split(" ")
                break
        table_path = coeff[0]
        section_name = coeff[1]
        table = LAMMPSPotentialTableIO(table_path, section_name).read()
        x = table.x
        U = table.E
        return list(x), list(U)
