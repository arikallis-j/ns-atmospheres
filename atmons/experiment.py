from .function import *
from .classes import *

def test_cli(name='User', *args, **kwargs):
    print(f"Hello, {name}!")

def experiment_format(exp_name, **param):
    """Formatting to origin printing"""
    def decorater(experiment):
        def wrapper(*args, **kwargs):
            print(f"{'#' * (len(exp_name) + 11)}\nEXPERIMENT {exp_name}\n{'#' * (len(exp_name) + 11)}")
            kwargs['experiment'] = exp_name
            for key, value in param.items():
                kwargs[key] = value 
            data = experiment(*args, **kwargs)
            return data
        return wrapper
    return decorater

def experiment_part_format(**param):
    """Formatting to origin printing"""
    def decorater(experiment):
        def wrapper(*args, **kwargs):
            part_name = ""
            data = []
            for key, value in param.items():
                if not isinstance(value, list):
                    part_name += f"\n# {key}: {value}"
                    kwargs[key] = value 
            do_exp = False
            for key, value in param.items():
                if isinstance(value, list):
                    for par in value:
                        do_exp = True
                        current_name = part_name + f"\n# {key}: {par}"
                        kwargs[key] = par 
                        print(current_name)
                        d = experiment(*args, **kwargs)
                        data.append(d)
            
            if not do_exp:
                print(part_name)
                d = experiment(*args, **kwargs)
                data.append(d)

            return data

        return wrapper
    return decorater

def experiment_item_format(**param):
    """Formatting to origin printing"""
    def decorater(experiment):
        def wrapper(*args, **kwargs):
            item_name = ""
            data = []
            for key, value in param.items():
                if not isinstance(value, list):
                    item_name += f"\n{key}: {value}"
                    kwargs[key] = value
                     
            do_exp = False
            for key, value in param.items():
                if isinstance(value, list):
                    for par in value:
                        do_exp = True
                        current_name = item_name + f"\n{key}: {par}"
                        kwargs[key] = par 
                        print(current_name)
                        print(f"{'-' * 64}")
                        d = experiment(*args, **kwargs)
                        data.append(d)
                        print(f"{'-' * 64}")
            
            if not do_exp:
                print(item_name)
                print(f"{'-' * 64}")
                d = experiment(*args, **kwargs)
                data.append(d)
                print(f"{'-' * 64}")
            return data
        return wrapper
    return decorater

def do_exp(exp_func, param, *args, **kwargs):
    for n in range(len(param)):
        part_param = param[n]['part']
        item_param = param[n]['item']
        @experiment_part_format(**part_param)
        @experiment_item_format(**item_param)
        def experiment(*args, **kwargs):
            data = exp_func(*args, **kwargs)
            return data

        data = experiment(*args, **kwargs)
        return data

class Experiment:
    def __init__(self, name, func=test_cli, **param):
        self.name = name
        self.func = func
        self.param = param

    def do_experiment(self, *args, **kwargs):
        data = self.func(*args, **kwargs)
        return data

    def __call__(self):
        descript = self._create_descript()

        @experiment_format(self.name, **self.param)
        def my_exp(descript, **kwargs):
            data = do_exp(self.do_experiment, descript, **kwargs)
            return data
            
        data = my_exp(descript)
        return data

    def __str__(self):
        descript = f"# EXPERIMENT {self.name} #\n"
        for name, stage in self.__dict__.items():
            if isinstance(stage, Stage):
                descript += f"{name}:\n"
                for key, value in stage.part.items():
                    descript += f"\t{key}: {value}\n"
                for key, value in stage.item.items():
                    descript += f"\t\t{key}: {value}\n"
        return descript
    
    def _create_descript(self):
        descript = []
        for name, stage in self.__dict__.items():
            if isinstance(stage, Stage):
                descript.append({
                    'part': stage.part,
                    'item': stage.item,
                })
        return descript
    
class Stage:
    def __init__(self, part = {}, item = {}):
        self.part = part
        self.item = item

class OneStage(Stage):
    def __init__(self, **param):
        super().__init__()
        self.item = param

class NoneModel(OneStage):
    def __init__(self, **param):
        super().__init__()
        self.item['w_func'] = 'none'
        self.item = param

class Const(OneStage):
    def __init__(self, **param):
        super().__init__()
        self.item['w_func'] = 'const'
        self.item = param

class Sqrt(OneStage):
    def __init__(self, **param):
        super().__init__()
        self.item['w_func'] = 'sqrt'
        self.item = param

class Line(OneStage):
    def __init__(self, **param):
        super().__init__()
        self.item['w_func'] = 'line'
        self.item = param

class Power(OneStage):
    def __init__(self, n=1, **param):
        super().__init__()
        self.item['w_func'] = f'power-{n}'
        self.item = param