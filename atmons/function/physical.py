"""
Description of phisical relationships
"""
from ..const import *
from .math import *

# Inercia and Momentum
def I_NS(i: Q[DIMLESS], R: Q[CM], M: Q[GRAM]) -> Q[CM**2 * GRAM]:
    """Neutron star moment of inertia (II.eq.A.10-A.11)"""
    I = i * M * R**2 
    return I << CM**2 * GRAM

def J_NS(I:Q[CM**2 * GRAM], omega: Q[SEC**(-1)]) -> Q[CM**2 * GRAM / SEC]: 
    """Neutron star angular momentum (II.eq.A.10-11)"""
    J = I * omega
    return J << CM**2 * GRAM / SEC

# Schwarzschild parameters
def R_sch(M: Q[GRAM]) -> Q[CM]:
    """Schwarzschild radius (by definition)"""
    r_sch = 2.0 * G_GRAV * M / C_LIGHT**2 
    return r_sch << CM

def zsch(R: Q[CM], R_s: Q[CM]) -> Q[DIMLESS]: 
    """Schwarzschild red shift plus one (I.eq.3)"""
    z = 1.0 / sqrt(1.0 - R_s / R)
    return z << DIMLESS

def Surf(R: Q[CM], zst: Q[DIMLESS]) -> Q[CM**2]:
    """Surface area (by relativity's definition)"""
    surf = (R * zst)**2
    return surf << CM**2

# Average gravity 
def g(R: Q[CM], M: Q[GRAM], zst: Q[DIMLESS]) -> Q[CM / SEC**2]:
    """Acceleration of freefall (by relativity's definition)"""
    grv = (G_GRAV * M / R**2) * zst
    return grv << CM / SEC**2
