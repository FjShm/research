from module.contents.functions.evaluate_functions import non_bonded, bonded


class EvaluateFunctions:
    def __call__(self, eval_type: str):
        if eval_type == "non_bonded":
            return non_bonded
        elif eval_type == "bonded":
            return bonded
