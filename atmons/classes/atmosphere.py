from ..funcs import *
from .phenomenon import *

class Atmosphere(Phenomenon):
    """Description of Atmosphere parameters"""
    def __init__(self, cfg = None, grid = None, spec = None, body = None, sp_layer = None, surface = None):
        if cfg is None or grid is None or spec is None or body is None or sp_layer is None or surface is None:
            return None
        
        # config
        self.luminosity = cfg.lum
        self.n_model = cfg.n_model

        # grid and spec param
        phi, theta, R, dR = grid.phi, grid.theta, grid.R, grid.dR

        E, dE = spec.E, spec.dE
        w_b, T_c, B_mod = spec.w_b, spec.T_c, spec.B_mod
        spec_key = spec.spec_key

        self.area_0, self.Lum_obs =  body.area_0, body.Lum_obs
        self.Flux_edd, self.spread_layer_true = surface.Flux_edd, surface.spread_layer_true
        
        # surface parameters
        self.E = np.tile(spec.E[None, None, :], (*phi.shape, 1))
        self.dE = np.tile(spec.dE[None, None, :], (*phi.shape, 1))

        self.nu = np.tile(surface.nu[:, :, None], (1, 1, *E.shape))
        self.delta = np.tile(surface.delta[:, :, None], (1, 1, *E.shape))
        self.beta_ph = np.tile(surface.beta_ph[:, :, None], (1, 1, *E.shape))
        self.cos_xi = np.tile(surface.cos_xi[:, :, None], (1, 1, *E.shape))
        self.cos_sig = np.tile(surface.cos_sig[:, :, None], (1, 1, *E.shape))
        self.cos_sig_1 = np.tile(surface.cos_sig_1[:, :, None], (1, 1, *E.shape))

        self.E_real = E_rad(self.E, self.nu, self.delta, self.beta_ph, self.cos_xi)
        self.kappa = kappa_E_rad(self.nu, self.delta, self.beta_ph, self.cos_xi)

        self.grv = np.tile(surface.grv[:, :, None], (1, 1, *E.shape))
        self.spread_layer = np.tile(surface.spread_layer[:, :, None], (1, 1, *E.shape))
        
        self.dOmega_obs = np.tile(surface.dOmega_obs[:, :, None], (1, 1, *E.shape))
        self.dOmega = np.tile(surface.dOmega[:, :, None], (1, 1, *E.shape))
    
        # radiational parameters
        self.flux = self.luminosity * DIMLESS
        self.Flux = self.flux * body.Flux_edd
        self.T_eff = T_SB(self.Flux)
        self.Epsilon_eff = Epsilon(T_obs(self.Flux, body.zsch))
        
        self.flux, self.t_eff, self.n_model = self._flux_T_eff(cfg.flux_key, self.flux, self.T_eff, surface.Flux_edd_base, surface.Flux_edd, self.n_model)
        
        # spectra 
        if spec_key == 'wfc':
            wwf_T, tcf_T = wfc_inter(T_c, w_b, self.flux, surface.log_g)
            T_c_int = np.tile(tcf_T[:, :, None], (1, 1, *E.shape))
            w_b_int = np.tile(wwf_T[:, :, None], (1, 1, *E.shape))
            self.rho = rho_rad(self.E_real, T_c_int, w_b_int, spectrum="planc")
        elif spec_key == 'be':
            self.rho = B_inter(B_mod, self.flux, surface.log_g, self.E_real)
            
        self.I_e = I_e_rad(self.rho, self.cos_sig_1)
        self.B_Omega = B_Omega_rad(self.I_e, self.kappa)
        self.B_int = self.B_Omega * self.dOmega_obs
        self.B_int_real = np.where(np.logical_not(self.cos_sig < 0.0), self.B_int, np.zeros(self.B_int.shape))     

        calc_B = self.calc_B(self.B_int_real, self.E, self.dE, self.Lum_obs)
        self.B_real, self.flux_real, self.Lum, self.lum, self.E_null, self.dE_null = calc_B
        self.w, self.fc = self.calc_wfc(self.B_real, self.E, body.area_0, self.Epsilon_eff)

        #additional full
        self.calc_B_full(body.Lum_obs, surface.Flux_edd)

        # additional sl
        self.calc_B_sl(body.Lum_obs, surface.spread_layer_true)

    def calc(self):
        calc_B = self.calc_B(self.B_int_real, self.E, self.dE, self.Lum_obs)
        self.B_real, self.flux_real, self.Lum, self.lum, self.E_null, self.dE_null = calc_B
        self.w, self.fc = self.calc_wfc(self.B_real, self.E, self.area_0, self.Epsilon_eff)

        #additional full
        #self.calc_B_full(self.Lum_obs, self.Flux_edd)

        # additional sl
        #self.calc_B_sl(self.Lum_obs, self.spread_layer_true)

    def calc_wfc(self, B_real, E, area_0, Epsilon_eff):
        E_null = E[0,0,:]
        w, fc = w_fc_rad(area_0, Epsilon_eff, E_null, B_real)
        return w, fc
    
    def calc_B(self, B_int_real, E, dE, Lum_obs):
        E_null = E[0,0,:]
        dE_null = dE[0,0,:]
        B_real = np.sum(B_int_real, axis=(0,1)) 
        flux_real = 4.0 * PI * np.sum(B_int_real * dE, axis=2) / Lum_obs
        
        if len(B_real[B_real<0]) != 0:
            print("Incorrect Spectrum!")
            print(f"{len(B_real[B_real<0])} points have B_real < 0")

        Lum = 4.0 * PI * np.sum(B_int_real * dE)
        lum = Lum / Lum_obs

        return B_real, flux_real, Lum, lum, E_null, dE_null

    def calc_B_full(self, Lum_obs, Flux_edd):
        self.B_int_full = self.B_Omega * self.dOmega
        self.B_full = np.sum(self.B_int_full, axis=(0,1)) 
        self.flux_full = 4.0 * PI * np.sum(self.B_int_full * self.dE, axis=2) / Lum_obs
        self.Lum_full = 4.0 * PI * np.sum(self.B_int_full * self.dE)
        self.lum_full = self.Lum_full / Lum_obs
        self.Flux_full = Flux_edd * self.flux
        self.Flux_all = np.sum(self.Flux_full)
    
    def calc_B_sl(self, Lum_obs, spread_layer_true):
        self.B_int_real_sl = np.where(self.spread_layer, self.B_int_real, self.B_int_real * 0.0)
        self.flux_real_sl = 4.0 * PI * np.sum(self.B_int_real_sl * self.dE, axis=2) / Lum_obs
        self.Flux_full_sl = np.where(spread_layer_true, self.Flux_full, self.Flux_full * 0.0)
        self.Flux_sl = np.sum(self.Flux_full_sl)
        self.xi_sl = self.Flux_sl / self.Flux_all
        self.B_real_sl = np.sum(self.B_int_real_sl, axis=(0,1)) 
        self.Lum_sl = 4.0 * PI * np.sum(self.B_int_real_sl * self.dE)
        self.lum_sl = self.Lum_sl / Lum_obs
    
    def _flux_T_eff(self, flux_key, flux_NS, T_eff_NS, Flux_edd, Flux_edd_sl, N_model):
        if flux_key == 'rel':
            T_eff = T_SB(flux_NS * Flux_edd)
            flux = np.ones(Flux_edd.shape) * flux_NS
        elif flux_key == 'abs':
            T_eff = np.ones(Flux_edd.shape) * T_eff_NS 
            Flux = Flux_SB(T_eff)
            flux = Flux / Flux_edd
            if (flux >= FLUX_REL[N_model-1]).any():
                N_model -= 1
                print("incorrectly flux")
        elif flux_key == 'sl':
            Flux_0 = flux_NS * Flux_edd
            flux = Flux_0 / Flux_edd_sl 
            T_eff = T_SB(Flux_0)
        else:
            print("Invalid flux_key value")

        return flux << u.Unit(), T_eff, N_model