from ..funcs import *
from .phenomenon import *

class SpreadLayer(Phenomenon):
    """Description of SpreadLayer parameters"""
    def __init__(self, cfg=None,  body=None):
        if cfg is None or body is None:
            return None
        
        self.w_func = cfg.w_func
        self.w_par = cfg.w_par
        self.th_star = cfg.th_star * DEG << RAD
        self.possible_th = cfg.possible_th

        # different factors
        kep_part = self.w_par
        chi_omega = body.omega_rot / body.omega_kep
        kep_part = max(kep_part, chi_omega)

        self.kep_part = kep_part
        self.omega_kep_local = body.omega_kep * kep_part
        self.rel_omega = self.omega_kep_local/body.omega_rot