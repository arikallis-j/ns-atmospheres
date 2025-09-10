"""
Description of chemical relationships
"""
from ..const import *
from .math import *

# Hydrogen part and scattering opacity
def X_hyd(chem: str) -> Q[DIMLESS]:
    """Hydrogen mass fraction (by observation)"""
    if chem=="he":
        x_hyd = 0.0
    elif chem=="s1" or chem=="s001":
        x_hyd = 0.7374
    else:
        x_hyd = 0.0 
    return x_hyd << DIMLESS
    
def kappa_e(x_hyd: Q[DIMLESS]) -> Q[CM**2 / GRAM]:
    """Coherent Thomson electron scattering opacity (I.eq.8)"""
    kappa = 0.2 * (1.0 + x_hyd) * CM**2 / GRAM
    return kappa << CM**2 / GRAM