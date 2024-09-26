import numpy as np

from .metric import * 
from .radiacional import *     
from .rotational import *       
from .model import * 

myint = np.int64
myfloat = np.float64
mychar = np.str_


class NS_SurfaceNew(Surface):
    def __init__(self, NS, size, rng):
        self.NS = NS
        self.r_func = NS.R_func
        self.r_0, self.R_0 = NS.r_eq, NS.R_eq 

        ph_range, th_range, nu_range = rng
        N_ph, N_th, N_nu = size

        # УМЕНЬШИТЬ РАЗМЕРНОСТЬ ПОСЛЕ ИСПОЛЬЗОВАНИЯ
        # ЛУЧШЕ:
        # - все рабочие переменные без self
        # - после использования -- c self
        self.phi_init = np.full((N_th, N_ph), np.linspace(*ph_range, N_ph)).T
        self.theta_init = np.full((N_ph, N_th), np.linspace(*th_range, N_th))

        self.r_init = self.r_func(self.r_0, self.phi_init, self.theta_init, self.NS)
        self.phi = self.phi_init * RAD
        self.theta = self.theta_init * RAD
        self.R = self.r_init * KM

        self.ph_range = ph_range[0] * RAD, ph_range[1] * RAD
        self.th_range = th_range[0] * RAD, th_range[1] * RAD
        self.sin_th, self.cos_th = sin(self.theta), cos(self.theta)
        self.sin_ph, self.cos_ph = sin(self.phi), cos(self.phi)

        self.dR = dR_metric(self.R_0, self.sin_th, self.cos_th, self.NS.chi, self.NS.Omega)   
        self.u = u_metric(self.R, self.NS.R_sch_cor)
        self.r_bar, self.u_bar = r_u_metric(self.R, self.cos_th, self.NS.q_c, self.NS.b_c, self.NS.R_sch_cor)
        self.nu, self.B, self.zeta = nu_B_dzeta_metric(self.cos_th, self.u_bar, self.NS.q_c, self.NS.b_c)
        self.omega_bar = omega_bar_metric(self.r_bar, self.u_bar, self.NS.J)
        self.beta_ph = beta_ph_metric(self.R, self.sin_th, self.omega_bar, self.nu)
        self.g_th = g_metric(self.sin_th, self.cos_th, self.NS.chi, self.NS.Omega)
        self.f_th = f_theta(self.R, self.dR, self.nu, self.B, self.zeta)
        self.sin_eta, self.cos_eta = eta_metric(self.f_th)
        self.beta = beta_metric(self.R, self.sin_th, self.nu, self.omega_bar, self.NS.omega_rot)
        self.gamma = gamma_metric(self.beta)

        # Gravity
        self.grv = grv_metric(self.theta, self.g_th, self.NS.W, self.NS.w_args, self.NS.g_0)
        self.log_g = log(self.grv)

        # rotational
        self.sin_psi, self.cos_psi = psi_rot(self.sin_th, self.cos_th, self.cos_ph, self.NS.sin_i, self.NS.cos_i)
        self.G_yu = G_yu_rot(self.cos_psi, self.u)
        self.D = D_rot(self.cos_psi, self.u)
        self.sin_a, self.cos_a = alpha_rot(self.cos_psi, self.u, self.G_yu)
        self.cos_chi = chi_rot(self.sin_th, self.cos_th, self.sin_psi, self.cos_psi, self.NS.cos_i)
        self.cos_sig = sigma_rot(self.sin_eta, self.cos_eta, self.sin_a, self.cos_a, self.cos_chi, self.cos_th)
        self.cos_xi = xi_rot(self.sin_a, self.sin_psi, self.sin_ph, self.NS.sin_i)
        self.delta = delta_rot(self.beta, self.gamma, self.cos_xi)
        self.cos_sig_1 = sigma_1_rot(self.cos_sig, self.delta)
        
        # # radiational
        self.Flux_edd_real = Flux_edd(self.grv, self.NS.sigma_T)
        self.flux, self.T_eff = T_flux_eff(self.NS.flux_key, self.NS.flux, self.NS.T_eff, self.Flux_edd_real)
        self.wwf_T, self.tcf_T = wwf_tcf_T(self.NS.T_c, self.NS.w_b, self.flux, self.log_g)
        

        # spectra 
        E, dE = E_base(N_nu, nu_range)
        self.E = np.full((N_ph, N_th, N_nu), np.array(E))
        self.dE = np.full((N_ph, N_th, N_nu), np.array(dE))
        
        nu_E = np.full((N_nu, N_th, N_ph), self.nu.T).T
        delta_E = np.full((N_nu, N_th, N_ph), self.delta.T).T
        beta_ph_E = np.full((N_nu, N_th, N_ph), self.beta_ph.T).T
        cos_xi_E = np.full((N_nu, N_th, N_ph), self.cos_xi.T).T
        cos_sig_1_E = np.full((N_nu, N_th, N_ph), self.cos_sig_1.T).T
        tcf_T_E = np.full((N_nu, N_th, N_ph), self.tcf_T.T).T
        wwf_T_E = np.full((N_nu, N_th, N_ph), self.wwf_T.T).T
        self.E_real = E_rad(self.E, nu_E, delta_E, beta_ph_E, cos_xi_E)
        self.rho = rho_rad(self.E_real, tcf_T_E, spectrum="planc")
        self.I_e = I_e_rad(self.rho, wwf_T_E, cos_sig_1_E)
        self.kappa_E = kappa_E_rad(nu_E, delta_E, beta_ph_E, cos_xi_E)
        self.B_Omega = B_Omega_rad(self.I_e, self.kappa_E)

        # integration
        self.surf = 0
        self.surf_real = 0
        self.dS = dS_metric_1(self.theta, self.cos_eta, self.R, N_ph, N_th, self.ph_range, self.th_range)
        self.dOmega = dOmega_rad(self.dS, self.cos_sig, self.D)
        dOmega_E = np.full((N_nu, N_th, N_ph), self.dOmega.T).T
        self.B_int = self.B_Omega * dOmega_E
        self.surf = np.sum(self.dS)

        self.dOmega_real = np.where(np.logical_not(self.cos_sig < 0.0), self.dOmega, np.zeros(self.dOmega.shape))
        self.surf_real = np.sum(self.dOmega_real)

        cos_sig_E = np.full((N_nu, N_th, N_ph), self.cos_sig.T).T
        self.B_int_real = np.where(np.logical_not(cos_sig_E < 0.0), self.B_int, np.zeros(self.B_int.shape))

        self.B_real = np.sum(self.B_int_real, axis=(0,1))

        ph_min, ph_max = ph_range
        l_phi = ph_max - ph_min

        self.surf /= l_phi
        self.R_pr = sqrt(self.surf)

        self.Lum = np.sum(self.B_real*self.dE)

        self.lum = 4.0 * PI * self.Lum / self.NS.Lum_obs

        self.E_null = self.E[0,0,:]

        self.w, self.fc = w_fc_rad(self.NS.surf_0, self.NS.Epsilon_eff, self.E_null, self.B_real)

class NeurtonStarNew:
    def __init__(self, name="J0000+0000", sys="base", chem="s1",rel=False, 
                 w_key="null", fc_key=1, flux_key=1,
                 r_ns=12.0, m_ns=1.5, v_rot=700.0, i_ang=60.0, lum=0.1):
        # key parameters
        self.name = name
        self.sys = sys
        self.chem = chem
        self.rel = rel
        self.w_key = w_key
        self.fc_key = fc_key
        self.flux_key = flux_key
        
        # init parameters
        self.r_ns = r_ns
        self.m_ns = m_ns
        self.v_rot = v_rot
        self.i_ang = i_ang
        self.lum = lum

        # M & R | M_cor & R_eq
        self.isRelative = self.rel
        if self.isRelative:
            self.r = inverse(r_eq, r_ns, (m_ns, v_rot), base=1.0, pw=5)
            self.m = inverse(m_cor, m_ns, (r_ns, v_rot), base=1.0, pw=5)
            self.r_eq, self.m_cor = r_ns, m_ns
        else:
            self.r, self.m = r_ns, m_ns
            self.r_eq, self.m_cor = 0.0, 0.0

        self.R = R_NS(self.r, self.m)
        self.M = M_NS(self.m)

        # model of surface
        self.R_func = R_metric

        # model of grav_eff
        self.W, self.w_args = W_null, ()
        
        # printing
        self.par = {
            "r_ns": self.r,
            "m_ns": self.m,
            "v_rot": self.v_rot,
            "i_ang": self.i_ang,
            "lum": self.lum,
            "chem": self.chem,
        }

        self.NS_Keys =  np.array([
            self.name,
            self.sys,
            self.chem,
            self.rel,
            self.w_key,
            self.fc_key,
            self.flux_key,
        ], dtype = mychar)

        self.NS_Init = np.array([
            self.r_ns,
            self.m_ns,
            self.v_rot,
            self.i_ang,
            self.lum,
            self.r, 
            self.m,
            self.R,
            self.M,
        ], dtype=myfloat)
        

    def init_paramters(self):
        # phisical
        ### ВОПРОС ПРО НЕСКОРРЕКТИРОВАННЫЙ РАДИУС ШВРАЦШИЛЬДА
        self.R_sch = R_sch(self.M)
        self.zsch = zsch(self.R, self.R_sch)
        self.g = g(self.R, self.M, self.zsch)
        self.log_g = log(self.g)
        self.surf_0 = Surf(self.R, self.zsch)

        # chemical
        self.X_hyd = X_hyd(self.chem)
        self.sigma_T = sigma_T(self.X_hyd)
        self.w_b = w_b(self.chem, self.fc_key)
        self.T_c = T_c(self.chem, self.fc_key)

        # photometrical
        self.Flux_edd = Flux_edd(self.g, self.sigma_T)
        self.T_edd = T_obs(self.Flux_edd, self.zsch)
        self.Theta_edd = Theta(self.T_edd)
        self.Epsilon_edd = Epsilon(self.T_edd)
        self.Lum_edd = Lumen(self.Flux_edd, self.R)
        self.Lum_obs = Lumen_obs(self.Lum_edd, self.zsch)

        # rotatinal
        self.nu_rot = self.v_rot 
        self.incl_ang = self.i_ang * RAD
        self.sin_i = sin(self.incl_ang)
        self.cos_i = cos(self.incl_ang)
        self.omega_rot = omega(self.nu_rot)   

        # relativical
        self.v_cr = nu_crit(self.r, self.m)
        self.v_rel = nu_relative(self.nu_rot, self.v_cr)
        self.m_cor = m_cor(self.m, self.r, self.nu_rot, define=self.isRelative, m=self.m_cor)
        self.M_cor = M_NS(self.m_cor)
        self.r_eq = r_eq(self.r, self.m, self.nu_rot, define=self.isRelative, r=self.r_eq)
        self.R_eq = R_NS(self.r_eq, self.m_cor)
        self.R_sch_cor = R_sch(self.M_cor)

        # metrical
        self.chi, self.Omega = chi_Omega_metric(self.R_eq, self.M_cor, self.omega_rot)
        self.q_c, self.b_c = q_b_metric(self.chi, self.Omega)
        self.i_bar = i_bar_metric(self.chi)
        self.I = I_NS(self.i_bar, self.R_eq, self.M_cor)
        self.J = J_NS(self.I, self.omega_rot)
        self.g_0 = g_0_metric(self.R_eq, self.M_cor, self.chi)

        
        # radiational
        self.flux = self.lum
        self.Flux = self.flux * self.Flux_edd
        self.T_eff = T_SB(self.Flux)
        self.Epsilon_eff = Epsilon(T_obs(self.Flux, self.zsch))


    def init_grid(self, size, rng, unnull=True):
        ph_range, th_range, nu_range = rng
        N_ph, N_th, N_nu = size

        if unnull:
            dphi =  (ph_range[1] - ph_range[0]) / (N_ph)
            dtheta = (th_range[1] - th_range[0]) / (N_th)
            ph_range = (ph_range[0] + dphi/2, ph_range[1] - dphi/2)
            th_range = (th_range[0] + dtheta/2, th_range[1] - dtheta/2)
        rng = ph_range, th_range, nu_range

        self.Surface = NS_SurfaceNew(self, size, rng)

        self.E, self.dE, self.B_real = self.Surface.E, self.Surface.dE, self.Surface.B_real
        self.Lum, self.lum = self.Surface.Lum, self.Surface.lum
        self.w, self.fc = self.Surface.w, self.Surface.fc
        self.surf, self.surf_real = self.Surface.surf, self.Surface.surf_real
        self.R_pr = self.Surface.R_pr

    def calc(self, n_phi=10, n_theta=10, n_nu=10, rng_phi=(0, 360), rng_theta=(0, 180), rng_erg=(0.1, 50.0), unnull=True):
        self.init_paramters()
        rng = rng_phi, rng_theta, rng_erg
        size = n_phi, n_theta, n_nu
        self.init_grid(size, rng, unnull=unnull)