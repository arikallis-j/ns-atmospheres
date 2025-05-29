from ..const import *
import json, os

class Phenomenon:
    def __init__(self):
        pass

    def __str__(self):
        return str(self.output()) # TODO: more fancy output
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__

    def save(self, name="data"):
        class_name = self.__class__.__name__.lower() 
        data = {}
        for key, parameter in self.__dict__.items():
            if isinstance(parameter, Q):
                if isinstance(parameter.value, np.ndarray):
                    data[key] = (parameter.value.tolist(), str(parameter.unit))
                else:   
                    data[key] = (parameter.value, str(parameter.unit))
            else:
                if isinstance(parameter, np.ndarray):
                    data[key] = (parameter.tolist(), "")
                else:   
                    data[key] = parameter
                
        if not os.path.isdir(f'results'):
            os.mkdir(f'results')

        with open(f"results/{name}_{class_name}.json", 'w') as f:
            json.dump(data, f)

    def load(self, name="data"):
        class_name = self.__class__.__name__.lower() 

        with open(f"results/{name}_{class_name}.json") as f:
            data = json.load(f)

        for key, parameter in data.items():
            if isinstance(parameter, tuple):
                self.__dict__[key] = parameter[0] * u.Unit(parameter[1])
            if isinstance(parameter, list):
                self.__dict__[key] = np.array(parameter[0]) * u.Unit(parameter[1])
            else:
                self.__dict__[key] = parameter
