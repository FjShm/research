from typing import Union
from module.adapter.get_Ui_adapter import GetUiAdapter


class GetUi:

    x: Union[list, None]
    U: Union[list, None]
    __InputData: dict

    def __init__(self, InputData: dict) -> None:
        self.x = []
        self.U = []
        self.__InputData = InputData

    def get_Ui(self):
        self.x, self.U = GetUiAdapter(
            function_type=self.__InputData["function_type"],
            previous_params_filepath=self.__InputData["previous_params_filepath"],
            Min=self.__InputData["min"],
            Max=self.__InputData["max"],
        ).request()
        return self
