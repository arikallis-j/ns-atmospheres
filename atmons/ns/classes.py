"""
Description of model classes
"""

import numpy as np

from .metric import * 
from .radiacional import *     
from .rotational import *       
from .model import * 

class NeurtonStarConfig:
    """???"""
    def __init__(self, config):
        self.name = config['name']

        # key parameters
        self.spec_key = config['spec_key']
        self.fc_key = config['fc_key']
        self.flux_key = config['flux_key']
        self.chem = config['chem']
        self.rel = config['rel']
        self.w_func = config['w_func']

        # init parameters
        self.r_ns = config['r_ns']
        self.m_ns = config['m_ns']
        self.v_rot = config['v_rot']
        self.i_ang = config['i_ang']
        self.th_star = config['th_star']
        self.w_par = config['w_par']

    def __str__(self):
        # TODO: more fancy output
        return str(self.output())
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__

class NeurtonStarParameters:
    """???"""
    def __init__(self, cfg):
        # M & R | M_cor & R_eq
        if cfg.rel:
            self.r = inverse(r_eq, (cfg.r_ns * KM), ((cfg.m_ns * M_SUN), (cfg.v_rot * u.Hz)), base=1.0, pw=5)
            self.m = inverse(m_cor, (cfg.m_ns * M_SUN), (cfg.r_ns * KM, cfg.v_rot * u.Hz), base=1.0, pw=5)
            self.r_eq, self.m_cor = cfg.r_ns, cfg.m_ns
        else:
            self.r, self.m = cfg.r_ns, cfg.m_ns
            self.r_eq, self.m_cor = 0.0, 0.0

        self.r <<= KM
        self.m <<= M_SUN 
        self.r_eq <<= KM
        self.m_cor <<= M_SUN
        
        self.R = R_NS(self.r, self.m)
        self.M = M_NS(self.m)

        # model of surface
        self.R_func = R_metric

        # model of grav_eff
        self.w_func, self.w_par = cfg.w_func, cfg.w_par
        self.th_star = cfg.th_star

        # phisical
        # TODO: ВОПРОС ПРО НЕСКОРРЕКТИРОВАННЫЙ РАДИУС ШВРАЦШИЛЬДА
        self.R_sch = R_sch(self.M)
        self.zsch = zsch(self.R, self.R_sch)
        self.g = g(self.R, self.M, self.zsch)
        self.log_g = log(self.g / self.g.unit)
        self.area_0 = Surf(self.R, self.zsch)

        # chemical
        self.X_hyd = X_hyd(cfg.chem)
        self.kappa_e = kappa_e(self.X_hyd)
        self.w_b = w_b(cfg.chem, cfg.fc_key)
        self.T_c = T_c(cfg.chem, cfg.fc_key)

        # photometrical
        self.Flux_edd = Flux_edd(self.g, self.kappa_e)
        self.T_edd = T_obs(self.Flux_edd, self.zsch)
        self.Theta_edd = Theta(self.T_edd)
        self.Epsilon_edd = Epsilon(self.T_edd)
        self.Lum_edd = Lumen(self.Flux_edd, self.R)
        self.Lum_obs = Lumen_obs(self.Flux_edd, self.R, self.zsch)

        # rotatinal
        self.nu_rot = cfg.v_rot * u.Hz
        self.incl_ang = (cfg.i_ang * u.deg).to(u.rad)
        self.sin_i = sin(self.incl_ang)
        self.cos_i = cos(self.incl_ang)
        self.omega_rot = omega(self.nu_rot)   

        # relativical
        self.v_cr = nu_crit(self.r, self.m)
        self.v_rel = nu_relative(self.nu_rot, self.v_cr)
        self.m_cor = m_cor(self.m, self.r, self.nu_rot, define=cfg.rel, m=self.m_cor)
        self.M_cor = M_NS(self.m_cor)
        self.r_eq = r_eq(self.r, self.m, self.nu_rot, define=cfg.rel, r=self.r_eq)
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

        self.V_kep = np.sqrt(self.g_0 * self.R_eq)
        self.V_rot = self.omega_rot * self.R_eq

    def __str__(self):
        # TODO: more fancy output
        return str(self.output())
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__
    
class GridConfig:
    def __init__(self, grid):
        # key parameters
        self.unnull = grid['unnull']
        self.zsch_key = grid['zsch_key']

        # size parameters
        self.n_phi = grid['n_phi']
        self.n_theta = grid['n_theta']
        self.n_nu = grid['n_nu']

        self.size = self.n_phi, self.n_theta, self.n_nu

        # rng parameters
        self.rng_phi = grid['rng_phi']
        self.rng_theta = grid['rng_theta']
        self.rng_erg = grid['rng_erg']
        
        if self.unnull:
            dphi =  (self.rng_phi[1] - self.rng_phi[0]) / (self.n_phi)
            dtheta = (self.rng_theta[1] - self.rng_theta[0]) / (self.n_theta)
            self.rng_phi = (self.rng_phi[0] + dphi/2, self.rng_phi[1] - dphi/2)
            self.rng_theta = (self.rng_theta[0] + dtheta/2, self.rng_theta[1] - dtheta/2)
        
        self.rng = self.rng_phi, self.rng_theta, self.rng_erg

    def __str__(self):
        # TODO: more fancy output
        return str(self.output())
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__

class NeurtonStarSurface:
    def __init__(self, grid, par):
        ph_range, th_range, nu_range = grid.rng
        N_ph, N_th, N_nu = grid.size

        if grid.zsch_key:
            nu_range = nu_range[0]*par.zsch, nu_range[1]*par.zsch

        # phisical
        self.r_func = par.R_func
        self.r_0, self.R_0 = par.r_eq, par.R_eq 
    
        self.phi_init = np.full((N_th, N_ph), np.linspace(*ph_range, N_ph)).T
        self.theta_init = np.full((N_ph, N_th), np.linspace(*th_range, N_th))
        self.theta_init_sym = np.where(self.theta_init < 90.0, self.theta_init, 180.0 - self.theta_init)
        
        # dimensionful   
        self.phi_init = self.phi_init * DEG
        self.theta_init = self.theta_init * DEG
        self.theta_init_sym = self.theta_init_sym * DEG

        self.phi = self.phi_init << RAD
        self.theta = self.theta_init << RAD
        self.theta_sym = self.theta_init_sym << RAD

        self.r_init = self.r_func(self.r_0, self.phi, self.theta_sym , par)
        self.R = self.r_init << CM

        self.ph_range = (ph_range[0]*DEG) << RAD, (ph_range[1] * DEG) << RAD
        self.th_range = (th_range[0]*DEG) << RAD, (th_range[1] * DEG) << RAD
        self.sin_th, self.cos_th = sin(self.theta), cos(self.theta)
        self.sin_ph, self.cos_ph = sin(self.phi), cos(self.phi)

        # metrical
        self.dR = dR_metric(self.R_0, self.sin_th, self.cos_th, par.chi, par.Omega)   
        self.u = u_metric(self.R, par.R_sch_cor)
        self.r_bar, self.u_bar = r_u_metric(self.R, self.cos_th, par.q_c, par.b_c, par.R_sch_cor)
        self.nu, self.B, self.zeta = nu_B_dzeta_metric(self.cos_th, self.u_bar, par.q_c, par.b_c)
        self.omega_bar = omega_bar_metric(self.r_bar, self.u_bar, par.J)
        self.beta_ph = beta_ph_metric(self.R, self.sin_th, self.omega_bar, self.nu)
        self.g_th = g_metric(self.sin_th, self.cos_th, par.chi, par.Omega)
        self.f_th = f_theta(self.R, self.dR, self.nu, self.B, self.zeta)
        self.sin_eta, self.cos_eta = eta_metric(self.f_th)
        self.beta = beta_metric(self.R, self.sin_th, self.nu, self.omega_bar, par.omega_rot)
        self.gamma = gamma_metric(self.beta)

        # Gravity
        self.grv = grv_metric(self.theta, self.g_th, par.g_0)
        self.log_g = log(self.grv / self.grv.unit)

        # rotational
        self.sin_psi, self.cos_psi = psi_rot(self.sin_th, self.cos_th, self.cos_ph, par.sin_i, par.cos_i)
        self.G_yu = G_yu_rot(self.cos_psi, self.u)
        self.D = D_rot(self.cos_psi, self.u)
        self.sin_a, self.cos_a = alpha_rot(self.cos_psi, self.u, self.G_yu)
        self.cos_chi = chi_rot(self.sin_th, self.cos_th, self.sin_psi, self.cos_psi, par.cos_i)
        self.cos_sig = sigma_rot(self.sin_eta, self.cos_eta, self.sin_a, self.cos_a, self.cos_chi, self.cos_th)
        self.cos_xi = xi_rot(self.sin_a, self.sin_psi, self.sin_ph, par.sin_i)
        self.delta = delta_rot(self.beta, self.gamma, self.cos_xi)
        self.cos_sig_1 = sigma_1_rot(self.cos_sig, self.delta)

        # spectra 
        E, dE = E_base(N_nu, nu_range)
        self.E = np.full((N_ph, N_th, N_nu), np.array(E))
        self.dE = np.full((N_ph, N_th, N_nu), np.array(dE))
        
        # dimensionful   
        self.E = self.E * KEV
        self.dE = self.dE * KEV

        self.nu_E = np.full((N_nu, N_th, N_ph), self.nu.T).T
        self.delta_E = np.full((N_nu, N_th, N_ph), self.delta.T).T
        self.beta_ph_E = np.full((N_nu, N_th, N_ph), self.beta_ph.T).T
        self.cos_xi_E = np.full((N_nu, N_th, N_ph), self.cos_xi.T).T
        self.cos_sig_1_E = np.full((N_nu, N_th, N_ph), self.cos_sig_1.T).T
        self.E_real = E_rad(self.E, self.nu_E, self.delta_E, self.beta_ph_E, self.cos_xi_E)
        self.kappa_E = kappa_E_rad(self.nu_E, self.delta_E, self.beta_ph_E, self.cos_xi_E)

        # integration
        self.dS = dS_metric_1(self.theta, self.cos_eta, self.R, N_ph, N_th, self.ph_range, self.th_range)
        self.area = np.sum(self.dS)
        
        ph_min, ph_max = self.ph_range
        l_phi = (ph_max - ph_min)

        self.R_pr = sqrt(self.area / (l_phi / RAD))
        
        self.dOmega_obs = dOmega_rad(self.dS, self.cos_sig, self.D)
        self.dOmega_obs_E = np.full((N_nu, N_th, N_ph), self.dOmega_obs.T).T 
        self.dOmega_obs_E <<= self.dOmega_obs.unit
        
        self.dOmega = self.dS
        self.dOmega_E = np.full((N_nu, N_th, N_ph), self.dOmega.T).T 
        self.dOmega_E <<= self.dOmega.unit

        self.theta_max = par.th_star * DEG << RAD
        self.G_eff = G_eff(self.theta, self.theta_max, par.w_func, par.w_par)
        self.G_eff_E = np.full((N_nu, N_th, N_ph), self.G_eff.T).T 
        self.Doppler = Doppler(self.theta, self.theta_max, self.phi, par.V_kep, par.V_rot, par.w_func, par.w_par)
        self.Doppler_E = np.full((N_nu, N_th, N_ph), self.Doppler.T).T 
        self.E_dop = self.E * self.Doppler_E
        self.dE_dop = self.dE * self.Doppler_E
        self.E_dop_real = self.E_real * self.Doppler_E

    def __str__(self):
        # TODO: more fancy output
        return str(self.output())
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__

class NeurtonStarShot:
    def __init__(self, lum, n_model, cfg, par, grid, surf):
        # radiational
        self.n_model = n_model
        self.flux = lum * u.Unit()
        self.Flux = self.flux * par.Flux_edd
        self.T_eff = T_SB(self.Flux)
        self.Epsilon_eff = Epsilon(T_obs(self.Flux, par.zsch))
        self.Flux_edd_real = Flux_edd(surf.grv, par.kappa_e)
        self.flux, self.T_eff, self.n_model = T_flux_eff(cfg.flux_key, self.flux, self.T_eff, self.Flux_edd_real, self.n_model)
        self.wwf_T, self.tcf_T = wwf_tcf_T(par.T_c, par.w_b, self.flux, surf.log_g)
        
        # spectra 
        self.tcf_T_E = np.full((grid.n_nu, grid.n_theta, grid.n_phi), self.tcf_T.T).T
        self.wwf_T_E = np.full((grid.n_nu, grid.n_theta, grid.n_phi), self.wwf_T.T).T

        self.tcf_T_E <<= self.tcf_T.unit
        self.wwf_T_E <<= self.wwf_T.unit

        if cfg.spec_key == 'wfc':
            self.rho = rho_rad(surf.E_real, self.tcf_T_E, self.wwf_T_E, spectrum="planc")
        elif cfg.spec_key == 'be':
            self.rho = B_inter(self.flux, surf.log_g, surf.E_real, cfg.chem)
        
        self.I_e = I_e_rad(self.rho, surf.cos_sig_1_E)
    
        self.B_Omega = B_Omega_rad(self.I_e, surf.kappa_E)

        self.B_int = self.B_Omega * surf.dOmega_obs_E
        self.B_int_full = self.B_Omega * surf.dOmega_E
        self.flux_full =  4.0 * PI * np.sum(self.B_int_full * surf.dE, axis=2) / par.Lum_obs

        self.dOmega_obs_real = np.where(np.logical_not(surf.cos_sig < 0.0), surf.dOmega_obs, np.zeros(surf.dOmega_obs.shape))
        self.area_real = np.sum(self.dOmega_obs_real)

        self.cos_sig_E = np.full((grid.n_nu, grid.n_theta, grid.n_phi), surf.cos_sig.T).T
        self.B_int_real = np.where(np.logical_not(self.cos_sig_E < 0.0), self.B_int, np.zeros(self.B_int.shape))
        
        self.flux_real = np.sum(self.B_int_real * surf.dE, axis=2)
        self.B_real = np.sum(self.B_int_real, axis=(0,1)) 
        self.Lum = 4.0 * PI * np.sum(self.B_int_real * surf.dE)
        self.lum = self.Lum / par.Lum_obs

        self.E_null = surf.E[0,0,:]

        self.w, self.fc = w_fc_rad(par.area_0, self.Epsilon_eff, self.E_null, self.B_real)
        
    def __str__(self):
        # TODO: more fancy output
        return str(self.output())
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__

class NeurtonStar:
    """
    Main class of the Neurton Star
    """ 
    def __init__(self, config, grid):
        self.n_model = N_MODEL
        self._init_config(config)
        self._init_grid(grid)
        self._init_ns()
        self._init_surface()
        self._init_shot()

        self.burst = self._burst
        self.burster = self._burster
        self.output = self.output()

    def __str__(self):
        # TODO: more fancy output
        return str(self._output())
    
    def _init_config(self, config):
        self.config = NeurtonStarConfig(config)

    def _init_grid(self, grid):
        self.grid = GridConfig(grid)
    
    def _init_ns(self):
        self.param = NeurtonStarParameters(self.config)
    
    def _init_surface(self):
        self.surface = NeurtonStarSurface(self.grid, self.param)
    
    def _init_shot(self):
        lum = FLUX_REL[0]
        self.shot = NeurtonStarShot(lum, 
            self.n_model, 
            self.config, 
            self.param, 
            self.grid, 
            self.surface)

        self.n_model = self.shot.n_model

    def _shot(self, l):
        lum = FLUX_REL[l]
        self.shot = NeurtonStarShot(lum, self.n_model, 
                                    self.config, 
                                    self.param, 
                                    self.grid, 
                                    self.surface)
        self.n_model = self.shot.n_model
        return self.shot
    
    def _burst(self):
        for l in range(self.n_model):
            shot = self._shot(l)
            if self.n_model < N_MODEL:
                break
            yield shot
    
    def _burster(self):
        shots = []
        for l in range(self.n_model):
            lum = FLUX_REL[l]
            shot = self._shot(lum)
            if self.n_model < N_MODEL:
                break
            shots.append(shot)
        return shots
    
    def __str__(self):
        # TODO: more fancy output
        return str(self.output())
    
    def __call__(self):
        return self.output()
    
    def output(self):
        return self.__dict__
    