from ..function import *
from .phenomenon import *

class Body(Phenomenon):
    """Description of Neutron Star parameters"""
    def __init__(self, cfg=None):
        if cfg is None:
            return None
        
        self.name = cfg.name

        # basic
        self.r_ns = cfg.r_ns * KM
        self.m_ns = cfg.m_ns * M_SUN
        self.v_rot =  cfg.v_rot * HZ
        self.chem = cfg.chem

        if cfg.rel:
            self.r = rel_r_eq(self.r_ns, self.m_ns, self.v_rot)
            self.m = rel_m_cor(self.m_ns, self.r_ns, self.v_rot)
            self.r_eq, self.m_cor = self.r_ns, self.m_ns
        else:
            self.r, self.m = self.r_ns, self.m_ns
            self.r_eq = r_eq(self.r, self.m, self.v_rot)
            self.m_cor = m_cor(self.m, self.r, self.v_rot)
        
        self.R = R_NS(self.r, self.m)
        self.M = M_NS(self.m)

        # phisical
        self.R_sch = R_sch(self.M) # TODO: QUESTION ABOUT NON-CORRECTED SCHWARZSCHILD RADIUS
        self.zsch = zsch(self.R, self.R_sch)
        self.area_0 = Surf(self.R, self.zsch)
        self.g = g(self.R, self.M, self.zsch)
        self.log_g = log(self.g / self.g.unit)

        # chemical
        self.X_hyd = X_hyd(self.chem)
        self.kappa_e = kappa_e(self.X_hyd)

        # photometrical
        self.Flux_edd = Flux_edd(self.g, self.kappa_e)
        self.T_edd = T_obs(self.Flux_edd, self.zsch)
        self.Theta_edd = Theta(self.T_edd)
        self.Epsilon_edd = Epsilon(self.T_edd)
        self.Lum_edd = Lumen(self.Flux_edd, self.R)
        self.Lum_obs = Lumen_obs(self.Flux_edd, self.R, self.zsch)

        # rotatinal
        self.nu_rot = self.v_rot
        self.incl_ang = (cfg.i_ang * DEG).to(RAD)
        self.sin_i = sin(self.incl_ang)
        self.cos_i = cos(self.incl_ang)
        self.omega_rot = omega(self.nu_rot)

        # relativical
        self.v_cr = nu_crit(self.r, self.m)
        self.v_rel = nu_relative(self.nu_rot, self.v_cr)
        self.M_cor = M_NS(self.m_cor)
        self.R_eq = R_NS(self.r_eq, self.m_cor)
        self.R_sch_cor = R_sch(self.M_cor)

        # metrical
        self.chi = chi_metric(self.R_eq, self.M_cor)
        self.Omega = Omega_metric(self.R_eq, self.M_cor, self.omega_rot)
        self.q_c = q_c_metric(self.chi, self.Omega)
        self.b_c = b_c_metric(self.chi, self.Omega)
        self.i_bar = i_bar_metric(self.chi)
        self.I = I_NS(self.i_bar, self.R_eq, self.M_cor)
        self.J = J_NS(self.I, self.omega_rot)
        self.g_0 = g_0_metric(self.R_eq, self.M_cor, self.chi)

        # keplerian
        self.V_rot = V_rot(self.omega_rot, self.R_eq)
        self.V_kep = V_kep(self.g_0, self.R_eq)
        self.omega_kep = omega_rot(self.V_kep, self.R_eq)