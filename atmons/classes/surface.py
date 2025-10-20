from ..funcs import *
from .phenomenon import *

class Surface(Phenomenon):
    """Description of Surface parameters"""
    def __init__(self, cfg=None, grid = None, body = None, sp_layer = None):
        if cfg is None or grid is None or body is None or sp_layer is None:
            return None
        
        phi, theta, R, dR = grid.phi, grid.theta, grid.R, grid.dR
        sin_ph, cos_ph = grid.sin_ph, grid.cos_ph
        sin_th, cos_th = grid.sin_th, grid.cos_th
        th_star = sp_layer.th_star

        # kappa max
        Omega_bs =  Omega_metric(body.R_eq, body.M_cor, body.omega_rot) 
        g_th_0 = g_metric(1, 0, body.chi, Omega_bs, Omega_bs, body.i_bar)
        kappa_teory = np.sqrt(G_GRAV*body.M_cor/(body.R_eq * body.V_kep**2))
        chi_i = (1 + body.chi * (-1 + 2*body.i_bar) + body.chi**2 * (-2 + 4*body.i_bar - 8*body.i_bar**2))
        kappa_teory = kappa_teory * np.sqrt(g_th_0/chi_i + Omega_bs**2)
        self.kappa_max = kappa_teory
        
        # kappa int
        Omega_bs =  Omega_metric(body.R_eq, body.M_cor, body.omega_rot) 
        g_th_null = 10.0**13.7 / body.g_0.value
        g_th_0 = g_metric(1, 0, body.chi, Omega_bs, Omega_bs, body.i_bar)
        kappa_teory = np.sqrt(G_GRAV*body.M_cor/(body.R_eq * body.V_kep**2))
        chi_i = (1 + body.chi * (-1 + 2*body.i_bar) + body.chi**2 * (-2 + 4*body.i_bar - 8*body.i_bar**2))
        kappa_teory = kappa_teory * np.sqrt((g_th_0 - g_th_null)/chi_i + Omega_bs**2)
        self.kappa_int = kappa_teory

        # kappa crit
        self.kappa_cr = body.omega_cr / body.omega_kep

        # keplerian part
        kep_part_local = sp_layer.kep_part
        if cfg.kep_part_fn=='max':
            kep_part_local = min(kep_part_local, self.kappa_max)
        if cfg.kep_part_fn=='int':
            kep_part_local = min(kep_part_local, self.kappa_int)
        if cfg.kep_part_fn=='cr':
            kep_part_local = min(kep_part_local, self.kappa_cr)
        if cfg.kep_part_fn=='none':
            kep_part_local = sp_layer.kep_part

        omega_kep_local = body.omega_kep * kep_part_local
        self.kep_part = kep_part_local

        # spread layer
        self.W_model = W_model(R, body.R_eq, theta, th_star, sp_layer.w_func, sp_layer.w_par, omega_kep_local, body.omega_rot)
        self.W_base = np.ones(self.W_model.shape) * self.W_model.unit
        self.omega_model = body.omega_rot * self.W_model
        self.omega_base = body.omega_rot * self.W_base

        self.psi = abs(90 * DEG - theta) << RAD
        self.spread_layer = self.psi <= th_star
        self.spread_layer_base = self.psi < 0.0
        self.spread_layer_true = self.spread_layer
        self.Omega_model = Omega_metric(body.R_eq, body.M_cor, self.omega_model)
        self.Omega_base = Omega_metric(body.R_eq, body.M_cor, self.omega_base)

        if sp_layer.w_func=='base':
            self.W_model = self.W_base
            self.omega_model = self.omega_base
            self.Omega_model = self.Omega_base
            self.spread_layer = self.spread_layer_base
            
        # metrical
        self.u = u_metric(R, body.R_sch)
        self.r_bar, self.u_bar = r_u_metric(R, cos_th, body.q_c, body.b_c, body.R_sch_cor)
        self.nu, self.B, self.zeta = nu_B_dzeta_metric(cos_th, self.u_bar, body.q_c, body.b_c)
        self.omega_bar = omega_bar_metric(self.r_bar, self.u_bar, body.J)
        self.beta_ph = beta_ph_metric(R, sin_th, self.omega_bar, self.nu)
        self.g_th = g_metric(sin_th, cos_th, body.chi, self.Omega_base, self.Omega_model, body.i_bar)

        g = self.g_th * body.g_0.value
        g = g * hs(g - 1) + 1 * hs(1 - g)
        log_g = log(g)
        log_g_crit = cfg.log_g_cr
        g_th_null = 10.0**log_g_crit / body.g_0.value
        if (log_g <= log_g_crit).any():
            print("Incorrect gravity")
        self.g_th = np.where(log_g > log_g_crit, self.g_th, self.g_th * 0.0 + g_th_null)
        self.g_th_base = g_metric(sin_th, cos_th, body.chi, self.Omega_base, self.Omega_base, body.i_bar)

        self.f_th = f_theta(R, dR, self.nu, self.B, self.zeta)
        self.sin_eta, self.cos_eta = eta_metric(self.f_th)
        self.beta = beta_metric(R, sin_th, self.nu, self.omega_bar, self.omega_model)
        self.gamma = gamma_metric(self.beta)

        # Gravity
        self.grv = grv_metric(theta, self.g_th, body.g_0)
        self.grv_base = grv_metric(theta, self.g_th_base, body.g_0)
        self.log_g = log(self.grv / self.grv.unit)
        self.log_g_base = log(self.grv_base / self.grv_base.unit)

        # radiational
        self.Flux_edd = Flux_edd(self.grv, body.kappa_e)
        self.Flux_edd_base = Flux_edd(self.grv_base, body.kappa_e)

        # rotational
        self.sin_psi, self.cos_psi = psi_rot(sin_th, cos_th, cos_ph, body.sin_i, body.cos_i)
        self.G_yu = G_yu_rot(self.cos_psi, self.u)
        self.D = D_rot(self.cos_psi, self.u)
        self.sin_a, self.cos_a = alpha_rot(self.cos_psi, self.u, self.G_yu)
        self.cos_chi = chi_rot(sin_th, cos_th, self.sin_psi, self.cos_psi, body.cos_i)
        self.cos_sig = sigma_rot(self.sin_eta, self.cos_eta, self.sin_a, self.cos_a, self.cos_chi, cos_th)
        self.cos_xi = xi_rot(self.sin_a, self.sin_psi, sin_ph, body.sin_i)
        self.delta = delta_rot(self.beta, self.gamma, self.cos_xi)
        self.cos_sig_1 = sigma_1_rot(self.cos_sig, self.delta)

        # integration
        self.dS = dS_metric_1(theta, self.cos_eta, R, grid.n_phi, grid.n_theta, grid.ph_range, grid.th_range)
        
        self.dOmega = self.dS
        self.dOmega_obs = dOmega_rot(self.dS, self.cos_sig, self.D)
        self.dOmega_obs_real = np.where(np.logical_not(self.cos_sig < 0.0), self.dOmega_obs, np.zeros(self.dOmega_obs.shape))
        
        self.area = np.sum(self.dOmega)
        self.area_real = np.sum(self.dOmega_obs_real)

        ph_min, ph_max = grid.ph_range
        self.l_phi = (ph_max - ph_min)

        self.R_pr = sqrt(self.area / (self.l_phi / RAD))
