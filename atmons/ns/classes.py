"""
Description of model classes
"""

import numpy as np

from .metric import * 
from .radiacional import *     
from .rotational import *       
from .model import * 

class NeurtonStar:
    """
    Main class of the Neurton Star
    """ 
    def __init__(self, config, grid):
        self.n_model = N_MODEL
        self.config = config
        self.grid = grid
        self.output = self.output()

    def __str__(self):
        # TODO: more fancy output
        return str(self.output())
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__