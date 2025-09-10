"""
Description of photometrical relationships
"""
from ..const import *
from .math import *

# Flux calculating
def Flux_edd(grv: Q[CM/SEC**2], kappa: Q[CM**2/GRAM]) -> Q[ERG/(SEC * CM**2)]:
    """Bolometric Eddington flux (I.eq.9)"""
    flux_edd = grv * C_LIGHT / kappa
    return flux_edd << ERG / (SEC * CM**2)

def Flux_SB(T: Q[KEL]) -> Q[ERG/(SEC * CM**2)]:
    """Flux of black body (by definition)"""
    flux_SB = SIGMA_SB * T**4
    return flux_SB << ERG / (SEC * CM**2)

# Temperature calculating
def T_SB(Flux: Q[ERG/(SEC * CM**2)]) -> Q[KEL]:
    """Temperature of black body (by definition)"""
    T = (Flux/SIGMA_SB)**(1/4)
    return T << KEL

def T_obs(Flux: Q[ERG/(SEC * CM**2)], zst: Q[DIMLESS]) -> Q[KEL]:
    """Temperature of observed black body (by relativity's definition)"""
    T = (Flux/SIGMA_SB)**(1/4) / zst
    return T << KEL

def Theta(T: Q[KEL]) -> Q[ERG]:
    """Temperature in ergs (by definition)"""
    theta = K_B * T
    return theta << ERG

def Epsilon(T: Q[KEL]) -> Q[KEV]:
    """Temperature in keV (by definition)"""
    epsilon = K_B * T
    return epsilon << KEV

# Luminosity calculating
def Lumen(Flux: Q[ERG/(SEC * CM**2)], R: Q[CM]) -> Q[ERG/SEC]: 
    """Luminosity (by definition)"""
    lumen = 4.0 * PI * R**2 * Flux
    return lumen << ERG / SEC

def Lumen_obs(Flux: Q[ERG/(SEC * CM**2)], R: Q[CM], zst: Q[DIMLESS]) -> Q[ERG/SEC]:
    """Observed luminosity (by relativity's definition)"""
    lumen = 4.0 * PI * R**2 * Flux / zst**2
    return lumen << ERG / SEC