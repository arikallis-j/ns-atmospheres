from ..function import *
from .phenomenon import *

class Surface(Phenomenon):
    """Description of Surface parameters"""
    def __init__(self, grid = None, body = None, sp_layer = None):
        if grid is None or body is None or sp_layer is None:
            return None
        
        phi, theta, R, dR = grid.phi, grid.theta, grid.R, grid.dR
        sin_ph, cos_ph = grid.sin_ph, grid.cos_ph
        sin_th, cos_th = grid.sin_th, grid.cos_th

        # spread layer
        self.W_model = W_model(R, body.R_eq, theta, sp_layer.th_star, sp_layer.w_func, sp_layer.w_par, sp_layer.omega_kep_local, body.omega_rot)
        self.W_base = np.ones(self.W_model.shape) * self.W_model.unit
        self.omega_model = body.omega_rot * self.W_model
        self.omega_base = body.omega_rot * self.W_base

        self.psi = abs(90 * DEG - theta) << RAD
        self.spread_layer = self.psi <= sp_layer.th_star
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
        self.u = u_metric(R, body.R_sch_cor)
        self.r_bar, self.u_bar = r_u_metric(R, cos_th, body.q_c, body.b_c, body.R_sch_cor)
        self.nu, self.B, self.zeta = nu_B_dzeta_metric(cos_th, self.u_bar, body.q_c, body.b_c)
        self.omega_bar = omega_bar_metric(self.r_bar, self.u_bar, body.J)
        self.beta_ph = beta_ph_metric(R, sin_th, self.omega_bar, self.nu)

        self.g_th = g_metric(sin_th, cos_th, body.chi, self.Omega_model)

        g = self.g_th * body.g_0.value
        g = g * hs(g - 1) + 1 * hs(1 - g)
        log_g = log(g)
        g_th_null = 10.0**13.7 / body.g_0.value
        self.g_th = np.where(log_g > 13.7, self.g_th, self.g_th * 0.0 + g_th_null)
        self.g_th_base = g_metric(sin_th, cos_th, body.chi, self.Omega_base)

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
