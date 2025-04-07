"""
Description of specific functions and wrappers
"""

import yaml

from .graph import *
from .math import *

def loader(name=""):
    """Loader of NS-congif"""
    if name=="":
        config = "__config__"
        grid = "__grid__"
    else:
        config = f"{name}_cfg"
        grid = f"{name}_grd"

    with open(f"configs/{config}.yml", 'r') as file:
        cfg = yaml.safe_load(file)
    
    with open(f"configs/{grid}.yml", 'r') as file:
        grd = yaml.safe_load(file)
    
    return cfg, grd
    
def dumper(name, config, grid):
    """Dumper of NS-congif and NS-grid"""
    with open(f"configs/{name}_cfg.yml", 'w') as file:
        yaml.dump(config, file)

    with open(f"configs/{name}_grd.yml", 'w') as file:
        yaml.dump(grid, file)

def input_config(from_name=""):
    """Inuput config and grid"""
    if from_name=="":
        config = "__config__"
        grid = "__grid__"
    else:
        config = f"{from_name}_cfg"
        grid = f"{from_name}_grd"

    with open(f"configs/{config}.yml", 'r') as file:
        parameters = yaml.safe_load(file)
    
    # TODO: inprove guide print
    print("\nPlease, input values in the same pattern:")
    for key, val in parameters.items():
        val_type = str(type(val))
        print(f"{key}: {val_type} = {parameters[key]}")
    
    print("\nYour Config:")
    for key, val in parameters.items():
        val_type = type(val)
        val_input = input(f"{key}: ")
        parameters[key] = val if val_input=="" else val_type(val_input)

    with open(f"configs/{grid}.yml", 'r') as file:
        griding = yaml.safe_load(file)
    
    # TODO: inprove guide print
    print("Please, input values in the same pattern:")
    for key, val in griding.items():
        val_type = str(type(val))
        print(f"{key}: {val_type} = {griding[key]}")
    
    print("\nYour Grid:")
    for key, val in griding.items():
        val_type = type(val)
        val_input = input(f"{key}: ")
        griding[key] = val if val_input=="" else val_type(val_input)

    return parameters, griding