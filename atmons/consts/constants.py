"""
Description of physical constants
    C_LIGHT: the speed of light
    G_GRAV: the gravitational constant 
    H_PL: Planck's constant
    SIGMA_SB: the Stefan-Boltzmann constant
    K_B: Boltzmann's constant
"""
import astropy.constants as astro_const

C_LIGHT = astro_const.c.cgs
G_GRAV = astro_const.G.cgs
H_PL = astro_const.h.cgs
SIGMA_SB = astro_const.sigma_sb.cgs
K_B = astro_const.k_B.cgs