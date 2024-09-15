from .math import *

def E_rad(E, nu, delta, beta_ph, cos_xi, th, ph):
    # (II.eq.B.35)
    N = len(E)
    E_real = [0.0]*N
    for k in range(0, N):
        E_real[k] = E[k] / (exp(nu) * delta)    / (1.0 + beta_ph * cos_xi)
        # if th == 0 and ph == 0:
        #     print(E[k], exp(nu), delta, beta_ph, cos_xi)
        #     print(E_real[k])
    return E_real

def rho_rad(E_real, tcf_T):
    N = len(E_real)
    rho = [0.0]*N
    for k in range(0, N):
        rho[k] = 1/K_PLANC * 2*PI / (C_LIGHT * H_PL)**2 * ERG * E_real[k]**3 / (exp(E_real[k] / tcf_T) - 1.0)
    return rho

def rho_one_rad(E_real, tcf_T):
    rho = 1/K_PLANC * 2*PI / (C_LIGHT * H_PL)**2 * ERG * E_real**3 / (exp(E_real / tcf_T) - 1.0)
    return rho
                        
def I_e_rad(rho, wwf_T, cos_sig_1):
    # (II.eq.9)
    N = len(rho)
    I_e = [0.0]*N
    for k in range(0, N):
        I_e[k] = wwf_T * rho[k] * (0.4215 + 0.86775 * cos_sig_1)
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
    N = len(I_e)
    B_Omega = [0.0]*N
    for k in range(0, N):
        B_Omega[k] = kappa_E**3 * I_e[k]
    return B_Omega

def B_real_rad(N, B_Omega, dOmega):
    # (II.eq.B.34)
    B_new = [0.0]*N
    for k in range(0, N):
        B_new[k] = B_Omega[k] * dOmega
    return B_new

def w_fc_rad(surf_0, θ_eff_0, E, B_real):
    w_min, w_max = 0.01, 1.05
    fc_min, fc_max = 1.0e-1, 3.0
    N_w, N_fc, N_f = 100, 100, len(E)
        
    itr = 0
    while not(itr>4):
        contsum = 1.e50
        dw = (w_max - w_min) / N_w
        dfc = (fc_max - fc_min) / N_fc    

        itr += 1

        for i in range(N_w):
            w_real = w_min + dw*i
            for j in range(N_fc):
                fc_real = fc_min + dfc*j
                cont = 0.0
                for nu in range(N_f):
                    if (E[nu]>3.0 and E[nu]<20.0):
                        rho_f = rho_one_rad(E[nu], fc_real*θ_eff_0)
                        fc_nu = w_real * PI * rho_f / E[nu]
                        cont += (fc_nu - B_real[nu] / surf_0 / E[nu])**2
                if (cont<contsum):
                    fc = fc_real
                    w = w_real
                    contsum = cont

        w_min = w - 2.0 * dw
        w_max = w + 2.0 * dw
        fc_min = fc - 2.0 * dfc
        fc_max = fc + 2.0 * dfc

    return w, fc