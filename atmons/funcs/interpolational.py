"""
Description of interpolational relationships
"""
from ..consts import *
from .math import *

def wfc_inter(T_c, w_b, flux, log_g):
    u_T, u_w = T_c.unit, w_b.unit

    T_c = np.tile(T_c[:, :, None, None], (1, 1, *log_g.shape))
    w_b = np.tile(w_b[:, :, None, None], (1, 1, *log_g.shape))

    T_g = map2(FLUX_REL, T_c, flux)
    w_g = map2(FLUX_REL, w_b, flux)

    T_g = np.tile(T_g[:, None, :, :], (1, 1, 1, 1))
    w_g = np.tile(w_g[:, None, :, :], (1, 1, 1, 1))

    T_f = map2(LOG_G, T_g, log_g)
    w_f = map2(LOG_G, w_g, log_g)

    T_f = T_f.reshape(log_g.shape)
    w_f = w_f.reshape(log_g.shape)

    wwf_T = w_f
    tcf_T = T_f
    return wwf_T << u_w, tcf_T << u_T

def B_inter(B_mod, flux, log_g, E):
    u_B = B_mod.unit
    
    B_int = B_mod << 1 / (SEC * CM**2)
    B_int = np.tile(B_mod[:, :, :, None, None], (1, 1, 1, *flux.shape))

    B_lum = map2(FLUX_REL_SHORT, B_int, flux)
    
    B_lg = map2(LOG_G, B_lum, log_g) 

    B_lg = np.tile(B_lg[:, None, :, :, None], (1, 1, 1, 1, E.shape[2]))
    
    B_en = map2(ENERGY, B_lg, E)
    B_en = B_en.reshape(E.shape)

    return B_en << u_B

