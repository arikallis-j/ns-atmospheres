import numpy as np

from .metric import * 
from .radiacional import *     
from .rotational import *       
from .model import * 

class NS_Point(Point):
    def __init__(self, NS, coord, phi, theta, r_0, r_func):
        super().__init__(NS, coord, phi, theta, r_0=r_0, r_func=r_func)

        # Point parameters
        self.NS = NS
        self.Surface = NS.Surface
        self.sin_th, self.cos_th = sin(self.theta), cos(self.theta)
        self.sin_ph, self.cos_ph = sin(self.phi), cos(self.phi)
        self.dS, self.dOmega, self.B_int = 0.0, 0.0, []
        
        # General Relativity
        ### ВОПРОС ПРО ЮЖНОЕ ПОЛУШАРИЕ
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

        # radiational
        self.Flux_edd_real = Flux_edd(self.grv, self.NS.sigma_T)
        self.flux, self.T_eff = T_flux_eff(self.NS.flux_key, self.NS.flux, self.NS.T_eff, self.Flux_edd_real)
        self.wwf_T, self.tcf_T = wwf_tcf_T(self.NS.T_c, self.NS.w_b, self.flux, self.log_g)

        # spectra 
        self.E_real = E_rad(self.Surface.E, self.nu, self.delta, self.beta_ph, self.cos_xi, self.th, self.ph)
        self.rho = rho_rad(self.E_real, self.tcf_T)
        self.I_e = I_e_rad(self.rho, self.wwf_T, self.cos_sig_1)
        self.kappa_E = kappa_E_rad(self.nu, self.delta, self.beta_ph, self.cos_xi)
        self.B_Omega = B_Omega_rad(self.I_e, self.kappa_E) 

class NS_Surface(Surface):
    def __init__(self, NS, r, r_func):
        super().__init__(NS, NS_Point, r=r, r_func=r_func)
        self.NS = NS
        self.E_size = 0
        self.E, self.dE, self.B_real = [], [], []
        self.Lum, self.lum = 0.0, 0.0
        self.w, self.fc = 0.0, 0.0

    def integrate(self):
        n_phi, n_theta = self.size
        n_f = self.E_size
        ph_min, ph_max = self.ranges[0]
        l_phi = ph_max - ph_min

        self.surf = 0
        self.surf_real = 0
        for i in range(n_theta):
            for j in range(n_phi):
                pnt = self.get_point(i, j)
                pnt.dS = dS_metric(pnt.th, pnt.ph, pnt.R, pnt.cos_eta, pnt.Surface)
                pnt.dOmega = dOmega_rad(pnt.dS, pnt.cos_sig, pnt.D)
                pnt.B_int = B_real_rad(n_f, pnt.B_Omega, pnt.dOmega)
                self.surf += pnt.dS
                if not(pnt.cos_sig < 0.0):
                    self.surf_real += pnt.dOmega
                    for k in range(n_f):
                        self.B_real[k] += pnt.B_int[k]
        
        self.surf /= l_phi
        self.R_pr = sqrt(self.surf)

        for k in range(n_f):
            self.Lum += self.B_real[k]*self.dE[k]  

        self.lum = 4.0 * PI * self.Lum / self.NS.Lum_obs

        self.w, self.fc = w_fc_rad(self.NS.surf_0, self.NS.Epsilon_eff, self.E, self.B_real)

class NeurtonStar(Star):
    def __init__(self, name="J0000+0000", sys="base", chem="s1",rel=False, 
                 w_key="null", fc_key=1, flux_key=1,
                 r_ns=12, m_ns=1.5, v_rot=700.0, i_ang=60.0, lum=0.1):
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

    def init_paramters(self):
        # phisical
        self.R = R_NS(self.r_ns, self.m_ns)
        self.M = M_NS(self.m_ns)
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

        self.Surface = NS_Surface(self, self.r_eq, r_func=self.R_func)

    def init_grid(self, n_phi=10, n_theta=10, n_nu=10, rng_phi=(0, 360), rng_theta=(0, 180), rng_erg=(0.1, 50.0), unnull=True):
        self.Surface.E, self.Surface.dE = E_base(n_nu, rng_erg)
        self.Surface.B_real = [0.0]*n_nu
        self.Surface.E_size = n_nu
        self.Size = (n_phi, n_theta, n_nu)

        super().init_grid(n_phi, n_theta, rng_phi, rng_theta, unnull)

    def integrate(self):
        self.Surface.integrate()

    def calc(self, n_phi=10, n_theta=10, n_nu=10, rng_phi=(0, 360), rng_theta=(0, 180), rng_erg=(0.1, 50.0), unnull=True):
        self.init_paramters()
        self.init_grid(n_phi, n_theta, n_nu, rng_phi, rng_theta, rng_erg, unnull)
        self.integrate()
        self.E, self.dE, self.B_real = self.Surface.E, self.Surface.dE, self.Surface.B_real
        self.Lum, self.lum = self.Surface.Lum, self.Surface.lum
        self.w, self.fc = self.Surface.w, self.Surface.fc
        self.surf, self.surf_real = self.Surface.surf, self.Surface.surf_real
        self.R_pr = self.Surface.R_pr

