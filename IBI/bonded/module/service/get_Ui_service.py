from typing import Union
from module.adapter.get_Ui_adapter import GetUiAdapter


class GetUi:

    x: Union[list, None]
    U: Union[list, None]
    InputData: dict

    def __init__(self, InputData: dict):
        self.x = []
        self.U = []
        self.InputData = InputData

    def get_Ui(self):
        self.x, self.U = GetUiAdapter(
            function_type=self.InputData["function_type"],
            previous_params_filepath=self.InputData["previous_params_filepath"],
            Min=self.InputData["min"],
            Max=self.InputData["max"],
        ).request()
