from ..consts import *
from ..funcs import *

class Body:
    def __init__(self, R_star, M_star, nu_star, i, C, atm=None):
        if atm is None:
            return None
        # config
        self.R_star = R_star
        self.M_star = M_star
        self.nu_star = nu_star
        self.i = i
        self.C = C

        #physical
        self.Omega_star = 2 * PI * nu_star << 1/SEC

        M_1_4 = M_star/(1.4*M_SUN)
        R_10 = R_star/(10*KM)
        self.nu_cr = 1278 * HZ * sqrt(M_1_4/R_10**3)

        self.nu_bar = self.nu_star / self.nu_cr 

        a1 = 0.001 * M_1_4**1.5
        a2 = 10*a1
        a0 = 1 - a1/1.1
        self.M_1 = self.M_star * (a0 - a1/(self.nu_bar - 1.1) + a2*(self.nu_bar**2))

        r1 = 0.025
        r2 = 0.07 * M_1_4**1.5
        r0 = 0.9766
        self.R_eq = self.R_star * (r0 - r1/(self.nu_bar - 1.07) + r2*(self.nu_bar**2))

        self.R_s = 2*G_GRAV*self.M_star/C_LIGHT**2 << KM
        self.R_s_1 = 2*G_GRAV*self.M_1/C_LIGHT**2 << KM
        
        self.zp1 = 1/sqrt(1 - self.R_s/self.R_star)
        self.g_star = G_GRAV * self.M_star / self.R_star**2  * self.zp1 << CM/SEC**2
        
        #photometrical
        self.kappa = 0.2 * CM**2 / GRAM * (1 + X[self.C])
        self.F_edd_star = self.g_star * C_LIGHT / self.kappa << ERG / (SEC * CM**2)
        self.T_edd_star = (self.F_edd_star / SIGMA_SB)**(1/4) / self.zp1 << KEL
        self.Theta_edd_star = K_B * self.T_edd_star << ERG
        self.L_edd_star = 4 * PI * self.R_star**2 * self.F_edd_star << ERG / SEC
        
        #metrical 
        self.chi = G_GRAV * self.M_1 / (self.R_eq * C_LIGHT**2) << DIMLESS
        self.Omega_bar_star = self.Omega_star * (self.R_eq**3/(G_GRAV * self.M_1))**(1/2) << DIMLESS
        self.I_bar = sqrt(self.chi) * (1.136 - 2.53*self.chi + 5.6 * self.chi**2) 
        self.I = self.I_bar * self.M_1 * self.R_eq**2
        self.J = self.I * self.Omega_star
        self.g_0 = G_GRAV * self.M_1 / (self.R_eq**2 * sqrt(1 - 2*self.chi)) << CM/SEC**2
        self.V_kep = (self.g_0 * self.R_eq)**(1/2) << CM/SEC
        self.Omega_kep = self.V_kep/self.R_eq << 1/SEC

class Surface:
    def __init__(self, phi, theta, Omega_norm, atm=None):
        if atm is None:
            return None
        #config
        self.phi = phi 
        self.theta = theta
        self.Omega = atm.Omega_star * Omega_norm

        a0 = -0.18 * atm.Omega_bar_star**2 + 0.23 * atm.chi * atm.Omega_bar_star**2 - 0.05 * atm.Omega_bar_star**4
        a2 = -0.39 * atm.Omega_bar_star**2 + 0.29 * atm.chi * atm.Omega_bar_star**2 - 0.13 * atm.Omega_bar_star**4
        a4 = +0.04 * atm.Omega_bar_star**2 - 0.15 * atm.chi * atm.Omega_bar_star**2 + 0.07 * atm.Omega_bar_star**4
        self.R = atm.R_eq * (1 + a0*P_0(cos(self.theta)) +  a2*P_2(cos(self.theta)) + a4*P_4(cos(self.theta)))
        self.dR = atm.R_eq * (a2 * dP_2(cos(self.theta), sin(self.theta)) + a4 * dP_4(cos(self.theta), sin(self.theta)))

        self.Omega_bar = self.Omega * (atm.R_eq**3/(G_GRAV*atm.M_1))**(1/2) << DIMLESS

        self.q_c = -0.11 * (self.Omega_bar / atm.chi)**2
        self.b_c = 0.4454 * self.Omega_bar**2 * atm.chi
        self.r_bar = self.R * 0
        err=1e-8
        r_bar = self.R
        P2 = P_2(abs(cos(theta))) << DIMLESS
        converged = np.full(self.theta.shape, False)
        while np.logical_not(converged.all()):
            r_mid = r_bar
            u_mid = atm.R_s / (2*r_mid)

            nu_0 = ln((1.0 - u_mid/2.0) / (1.0 + u_mid/2.0))
            B_0 = (1.0 - u_mid/2.0) * (1.0 + u_mid/2.0)

            nu = nu_0 + (self.b_c/3.0 - self.q_c * P2) * u_mid**3
            B = B_0 + self.b_c * u_mid**2

            r_bar = self.R / (exp(-nu) * B)

            errs = abs(r_bar - r_mid) / r_bar
            self.r_bar = np.where(np.logical_and(errs < err, np.logical_not(converged)), r_bar, self.r_bar)
            converged = np.where(errs < err, True, converged)
        
        self.u_bar = atm.R_s / (2 * self.r_bar) 
        nu_0 = ln((1.0 - self.u_bar/2.0) / (1.0 + self.u_bar/2.0))
        B_0 = (1.0 - self.u_bar/2.0) * (1.0 + self.u_bar/2.0)
        zeta_0 = ln(B_0)
        
        self.nu = nu_0 + (self.b_c/3.0 - self.q_c * P2) * self.u_bar**3
        self.B = B_0 + self.b_c * self.u_bar**2
        self.zeta = zeta_0 + self.b_c * (4/3 * P2 - 1/3.0) * self.u_bar**3

        self.varpi = 2 * G_GRAV * atm.J / (C_LIGHT**2 * self.r_bar**3) * (1 - 3 * self.u_bar) << 1/SEC
        self.beta_1 = self.R * self.varpi * exp(-self.nu) / C_LIGHT * sin(self.theta) << DIMLESS
        self.beta =  self.R * (self.Omega - self.varpi) * exp(-self.nu) / C_LIGHT * sin(self.theta) << DIMLESS
        self.gamma = 1 / sqrt(1 - self.beta**2)

        # gravity
        c_e = 0.776 * atm.chi - 0.791
        c_p = 1.138 - 1.431 * atm.chi
        d_e = (2.431 * atm.chi - 1.315) * atm.chi
        d_p = (0.653 - 2.864 * atm.chi) * atm.chi
        d_60 = (13.47 - 27.13 * atm.chi) * atm.chi
        f_e = -1.172 * atm.chi
        f_p = 0.975 * atm.chi
        e_e = c_e * atm.Omega_bar_star**2 + d_e * atm.Omega_bar_star**4 + f_e * atm.Omega_bar_star**6 
        e_p = c_p * atm.Omega_bar_star**2 + (d_e - d_60) * atm.Omega_bar_star**4 + f_p * atm.Omega_bar_star**6 
        e_60 = d_60 * atm.Omega_bar_star**4

        self.g = atm.g_0 * (1 + e_e * sin(theta)**2 + e_p * cos(theta)**2 + e_60 * abs(cos(theta)))
        
        lambda_chi = 1 + atm.chi * (1 - 2 * atm.I_bar) + atm.chi**2 * (-2 + 4 * atm.I_bar - 8 * atm.I_bar**2)
        self.g_eff = self.g - atm.g_0 * (self.Omega_bar - atm.Omega_bar_star) * sin(theta)**2 * lambda_chi
        
        self.log_g_eff = log(self.g_eff / self.g_eff.unit)

        self.F_edd = C_LIGHT * self.g_eff / atm.kappa << ERG / (SEC * CM**2)
        self.F_edd_base = C_LIGHT * self.g / atm.kappa << ERG / (SEC * CM**2)

        # rotational 
        self.cos_psi = cos(atm.i) * cos(self.theta) + sin(atm.i) * sin(self.theta) * cos(self.phi)
        self.sin_psi = sqrt(1.0 - self.cos_psi**2)
        y = 1 - self.cos_psi
        u = atm.R_s / self.R
        self.G_yu = 1.0 + (u * y)**2 / 112.0 - EULER * u * y / 100.0 * (ln(1.0 - 0.5 * y) + 0.5 * y) 
        self.cos_alpha = 1.0 - y * (1 - u) * self.G_yu
        self.sin_alpha = sqrt(1.0 - self.cos_alpha**2)
        self.D = 1.0 + 3.0 * (u * y)**2 / 112.0 - EULER * u * y / 100.0 * (2.0 * ln(1.0 - 0.5 * y) + y * (1.0 - 0.75 * y) / (1.0 - 0.5 * y))
        self.cos_mu = (cos(atm.i) - cos(self.theta) * self.cos_psi) / (sin(self.theta) * self.sin_psi)
        self.cos_xi = -1.0 * self.sin_alpha * sin(atm.i) * sin(self.phi) / self.sin_psi

        self.f = 1/sqrt(1 + atm.R_s_1/self.R) * self.dR/self.R
        self.sin_eta = self.f / sqrt(1 + self.f**2)
        self.cos_eta = 1 / sqrt(1 + self.f**2)
        
        self.cos_sigma = self.cos_eta * self.cos_alpha + self.sin_eta * self.sin_alpha * self.cos_mu 

        self.delta = 1.0 / (self.gamma * (1.0 - self.beta * self.cos_xi))
        self.cos_sigma_1 = self.delta * self.cos_sigma

        dcos_th = np.zeros(theta.shape)
        ph_min, ph_max = atm.ph_range * DEG
        th_min, th_max = atm.th_range * DEG
        dphi = (ph_max - ph_min) / (atm.N_ph - 1) 
        dtheta = (th_max - th_min) / (atm.N_th - 1)
        theta_0 = theta
        theta_p = theta + dtheta
        theta_m = theta - dtheta
        dcos_th = np.where(theta == PI/2*RAD - dtheta/2, (cos(theta_m) + cos(theta_0)) / 2.0 , dcos_th)
        dcos_th = np.where(theta == PI/2*RAD + dtheta/2, (cos(theta_0) + cos(theta_p)) / 2.0, dcos_th)
        dcos_th = np.where(theta == PI/2*RAD, (abs(cos(theta_m)) + abs(cos(theta_p))) / 2.0, dcos_th)
        dcos_th = np.where(np.logical_and(theta > th_min, theta < th_max), (cos(theta_m) - cos(theta_p)) / 2.0, dcos_th)
        dcos_th = np.where(theta == th_min, 1.0 - (cos(theta_0) + cos(theta_p)) / 2.0, dcos_th)
        dcos_th = np.where(theta == th_max, 1.0 - (abs(cos(theta_0) + cos(theta_m))) / 2.0, dcos_th)
        self.dS = 1 / self.cos_eta * self.R**2 * dcos_th * dphi / RAD << KM**2
        self.dOmega_obs = self.dS * self.cos_sigma * self.D

         
class Spectrum:
    def __init__(self, E, dE, f, T_eff, s_key, atm=None):
        if atm is None:
            return None
        # config
        self.s_key = s_key
        self.E = E
        self.dE = dE
        self.f = f
        self.T_eff = T_eff
        
        self.E = np.tile(self.E[None, None, :], (*atm.phi.shape, 1))
        self.dE = np.tile(self.dE[None, None, :], (*atm.phi.shape, 1))

        kappa_e = 1 / (exp(atm.nu) * atm.delta) / (1.0 + atm.beta_1 * atm.cos_xi)
        kappa_e = np.tile(kappa_e[:,:, None], (1, 1, atm.N_nu))
        self.E_1 = self.E * kappa_e

        self.w_b = w_b(atm.C, atm.fc_key)
        self.T_c = T_c(atm.C, atm.fc_key)
        self.B_mod = B_model(atm.C)

        if s_key == 'wfc':
            wwf_T, tcf_T = wfc_inter(self.T_c, self.w_b, self.f, atm.log_g_eff)
            T_c_int = np.tile(tcf_T[:, :, None], (1, 1, *E.shape))
            w_b_int = np.tile(wwf_T[:, :, None], (1, 1, *E.shape))
            self.rho = rho_rad(self.E_1, T_c_int, w_b_int, spectrum="planc")
        elif s_key == 'be':
            self.rho = B_inter(self.B_mod, self.f, atm.log_g_eff, self.E_1)

        cos_sigma_1 = atm.cos_sigma_1
        cos_sigma_1 = np.tile(cos_sigma_1[:,:, None], (1, 1, atm.N_nu))
        self.I = self.rho * (0.4215 + 0.86775 * cos_sigma_1)

        dOmega_obs = atm.dOmega_obs
        dOmega_obs = np.tile(dOmega_obs[:,:, None], (1, 1, atm.N_nu))
        dB_Omega = (self.E/self.E_1)**3 * self.I * dOmega_obs
        int_cond = np.logical_not(atm.cos_sigma < 0.0)
        int_cond = np.tile(int_cond[:,:, None], (1, 1, atm.N_nu))
        dB_Omega = np.where(int_cond, dB_Omega, np.zeros(dB_Omega.shape))
        self.B = np.sum(dB_Omega, axis=(0,1)) << ERG / (KEV * SEC)
        if len(self.B [self.B <0]) != 0:
            print("Incorrect Spectrum!")
            print(f"{len(self.B [self.B <0])} points have B_real < 0")


class Atmosphere:
    def __init__(self, r, m, nu, i, c, lum, s_key, f_key, fc_key=1, n_model=N_MODEL, N_ph=30, N_th=30, N_nu=500, ph_range = (0,360), th_range=(0,180), nu_range=(1,50)):
        self.body = Body(r*KM, m * M_SUN, nu*HZ, i*DEG, c, self)
        for key, value in self.body.__dict__.items():
            self.__dict__[key] = value
        
        self.N_ph, self.N_th, self.ph_range, self.th_range = N_ph, N_th, ph_range, th_range
        dphi =  (self.ph_range[1] - self.ph_range[0]) / (self.N_ph)
        dtheta = (self.th_range[1] - self.th_range[0]) / (self.N_th)
        self.ph_range = (self.ph_range[0] + dphi/2, self.ph_range[1] - dphi/2)
        self.th_range = (self.th_range[0] + dtheta/2, self.th_range[1] - dtheta/2)

        phi = np.linspace(*self.ph_range, self.N_ph) 
        theta = np.linspace(*self.th_range, self.N_th)
        phi, theta = np.meshgrid(phi, theta, indexing='xy') 
        Omega_norm = self.Omega_kep/self.Omega_star * np.ones(phi.shape)
        self.surf = Surface(phi * DEG, theta * DEG, Omega_norm, self)
        for key, value in self.surf.__dict__.items():
            self.__dict__[key] = value

        self.N_nu, self.nu_range = N_nu, nu_range
        self.fc_key, self.n_model = fc_key, n_model
        E, dE = E_base(self.N_nu, self.nu_range)
        E, dE = np.array(E), np.array(dE)
        T_eff = (lum * self.F_edd_star / SIGMA_SB)**(1/4) << KEL
        flux, t_eff, self.n_model = self._flux_T_eff(f_key, lum, T_eff, self.F_edd_base, self.F_edd, self.n_model)
        self.spec = Spectrum(E * KEV, dE * KEV, flux, t_eff, s_key, self)
        for key, value in self.spec.__dict__.items():
            self.__dict__[key] = value

    def _flux_T_eff(self, flux_key, flux_NS, T_eff_NS, Flux_edd, Flux_edd_sl, N_model):
        if flux_key == 'rel':
            T_eff = (flux_NS * Flux_edd / SIGMA_SB)**(1/4) << KEL
            flux = np.ones(Flux_edd.shape) * flux_NS
        elif flux_key == 'abs':
            T_eff = np.ones(Flux_edd.shape) * T_eff_NS 
            Flux =  SIGMA_SB * T_eff**4
            flux = Flux / Flux_edd
            if (flux >= FLUX_REL[N_model-1]).any():
                N_model -= 1
                print("incorrectly flux")
        elif flux_key == 'sl':
            Flux_0 = flux_NS * Flux_edd
            flux = Flux_0 / Flux_edd_sl 
            T_eff =  (Flux_0 / SIGMA_SB)**(1/4)
        else:
            print("Invalid flux_key value")

        return flux << u.Unit(), T_eff, N_model