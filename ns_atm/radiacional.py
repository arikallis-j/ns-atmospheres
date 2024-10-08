from .math import *

def E_rad(E, nu, delta, beta_ph, cos_xi):
    # (II.eq.B.35)
    E_real = E / (exp(nu) * delta) / (1.0 + beta_ph * cos_xi)
    return E_real

def rho_rad(E_real, tcf_T, spectrum="planc"):
    if spectrum=="planc":
        rho = 1/K_PLANC * 2*PI / (C_LIGHT * H_PL)**2 * ERG * E_real**3 / (exp(E_real / tcf_T) - 1.0)
    elif spectrum=="line":
        FWHM = 50
        sig = 1/(FWHM * 2 * sqrt(2 * ln(2)))
        rho = gaussian(E_real, mu = 1.0, sig=sig)
    return rho
                        
def I_e_rad(rho, wwf_T, cos_sig_1):
    # (II.eq.9)
    I_e = wwf_T * rho * (0.4215 + 0.86775 * cos_sig_1)
    return I_e
                        
def kappa_E_rad(nu, delta, beta_ph, cos_xi):
    # (II.eq.B.35)
    kappa_E = exp(nu) * delta * (1.0 + beta_ph * cos_xi)
    return kappa_E

def dOmega_rad(dS, cos_sig, D_rot):
    # (II.eq.B.33)
    dOmega =  dS * cos_sig * D_rot
    return dOmega

def B_Omega_rad(I_e, kappa_E):
    # (II.eq.B.34)
    B_Omega = kappa_E**3 * I_e
    return B_Omega


def w_fc_rad(surf_0, th_eff_0, E, B_real):
    w_min, w_max = 0.01, 1.05
    fc_min, fc_max = 1.0e-1, 3.0
    N_w, N_fc, N_f = 100, 100, len(E)
    fc, w = 0.0, 0.0
    itr = 1

    while not(itr>5):
        dw = (w_max - w_min) / N_w
        dfc = (fc_max - fc_min) / N_fc    

        rho_f = np.zeros(*E.shape)
        E_wfc = np.full(*E.shape, E)
        B_wfc = np.full(*B_real.shape, B_real)

        w_real = np.linspace(w_min, w_max - dw, N_w)
        w_real = np.full((N_fc, N_w), w_real).T
        w_real = np.full((*E.shape, N_w, N_fc), w_real).T
        fc_real = np.linspace(fc_min, fc_max - dfc, N_fc)
        fc_real = np.full((*E.shape, N_w, N_fc), fc_real).T

        itr += 1

        cont = np.zeros((N_w, N_fc))

        condicional_E = np.logical_and(E_wfc>3.0, E_wfc<20.0)
        rho_f = np.where(condicional_E, rho_rad(E_wfc, fc_real*th_eff_0), np.zeros(E_wfc.shape))
        fc_nu = np.where(condicional_E, w_real * PI * rho_f / E_wfc, np.zeros(E.shape))
        cont = np.sum(np.where(condicional_E, (fc_nu - B_wfc / surf_0 / E_wfc)**2, np.zeros(E_wfc.shape)), axis=2)
        fc_real = fc_real[:,:,0]
        w_real = w_real[:,:,0]
        ind = np.unravel_index(np.argmin(cont, axis=None), cont.shape)
        fc = fc_real[ind]
        w = w_real[ind]

        w_min = w - 2.0 * dw
        w_max = w + 2.0 * dw
        fc_min = fc - 2.0 * dfc
        fc_max = fc + 2.0 * dfc

    return w, fc
