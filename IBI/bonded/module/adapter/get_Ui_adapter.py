from typing import Union
from module.facade.get_Ui_table_facade import GetUiTableFacade
from module.facade.get_Ui_function_facade import GetUiFunctionFacade
from module.facade.get_Ui_none_facade import GetUiNoneFacade


class GetUiAdapter:

    function_type: str
    previous_params_filepath: str
    Min: float
    Max: float

    def __init__(
        self, function_type: str, previous_params_filepath: str, Min: float, Max: float
    ) -> None:
        self.function_type = function_type
        self.previous_params_filepath = previous_params_filepath
        self.Min = Min
        self.Max = Max

    def request(self):
        if self.previous_params_filepath is None:
            x, U = GetUiNoneFacade()()
        elif self.function_type == "table":
            x, U = GetUiTableFacade()(self)
        else:
            x, U = GetUiFunctionFacade()(self)
        return x, U
