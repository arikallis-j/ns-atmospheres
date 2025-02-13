"""
Description of the CLI commands
"""
import sys
import fire
from time import time

from .experiments import *
from ..ns.classes import *
from ..ns.data_load import *

def _guide() -> None:
    """Print a guide message."""
    print("This is empty guide now.")
  
def CLI() -> None:
    """Main function to run the calculator."""
    if len(sys.argv) == 1:
        _guide()
    else:
        calculator = Calculator()
        fire.Fire(calculator)

# CLI commands
class Calculator:
    def test(self)-> str:
        """This is a field for experiments"""
        return "Test is okey"
    
    def config(self, name, from_name=""):
        """Creation of NS-config"""
        config, grid = input_config(from_name)
        dumper(name, config, grid)
        return "Config and grid is dumped"

    def init(self, name="base"):
        """Inicialization of NS-object"""
        cfg, grid = loader(name)
        ns = NeurtonStar(cfg, grid)
        print(ns)
        return "Init is okey"

    def origin(self, name="base"):
        """Calculating of NS-burster"""
        cfg, grid = loader(name)
        ns = NeurtonStar(cfg, grid)
        @timer
        @origin_format(ns=ns)
        def burster():
            return ns.burst
        burster()

        return "Origin is okey"

    def calc(self, par, name="base", form='line', save=False):
        """Calculating of NS-burster"""
        cfg, grid = loader(name)
        @timer
        @print_format(form=form, save=save)
        @direction_calc(par)
        def burster():
            return cfg, grid
        burster()

        return "Calculation is okey"

# wrappers
def timer(func):
    """Timer for functions """
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        end = time()
        print(f'Timer: {end-start:.4f} seconds')
        return result
    return wrapper

def origin_format(ns):
    """Formatting to origin printing"""
    # TODO: correct format
    def decorater(func):
        def wrapper(*args, **kwargs):
            print("\n  surf/R_eq^2     rpr             R_eq:")
            print(f"  {sci(ns.surface.area / ns.param.R_eq**2, ndig=6)}     {sci(ns.surface.R_pr, ndig=6)}     {sci(ns.param.R_eq, ndig=6)}")
            print("\n         L/L_Edd       T_eff (keV)   f_c           w")
            shots = func(*args, **kwargs)
            num = 0
            for shot in shots:
                print(f"  {' ' if num<9 else ''}{num+1}     {shot.lum:.6f}      {shot.Epsilon_eff:.6f}      {shot.fc:.6f}      {shot.w:.6f}")
                num += 1
        return wrapper
    return decorater

def print_format(form, save):
    """Printing calculated data"""
    def decorater(func):
        def wrapper(*args, **kwargs):
            data = func(*args, **kwargs)
            graph = Graph(data)
            graph.draw(form, save)
        return wrapper
    return decorater

def direction_calc(par):
    """Direction calculated data"""
    def decorater(func):
        def wrapper(*args, **kwargs):
            cfg, grid = func(*args, **kwargs)
            if par=="lfc":
                return lfc(cfg, grid)
            if par=="lw":
                return lw(cfg, grid)
            if par=="spectra":
                return spectra(cfg, grid)
            else:
                print("Unavable command")
                return None
        return wrapper
    return decorater