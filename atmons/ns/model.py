from .math import *

def W_null(theta, args):
    return 0.0

def W_const(theta, args):
    return 0.75

# def W_line(theta, args):
#     theta_max = args
#     th_max = theta_max * RAD
#     W = (1 - theta/th_max) * hs(1 - theta/th_max)
#     return W
def W_none(theta, psi_max, par):
    W = np.zeros(theta.shape)
    return W

def W_const(theta, psi_max, par=0.5):
    C = par
    psi = abs(90 * DEG - theta) << RAD
    W = C * hs(1 - psi/psi_max)
    return W

def W_sqrt(theta, psi_max, par):
    psi = abs(90 * DEG - theta) << RAD
    W = (1 - (psi/psi_max)**(1/2)) * hs(1 - psi/psi_max)
    return W

def W_line(theta, psi_max, par):
    psi = abs(90 * DEG - theta) << RAD
    W = (1 - (psi/psi_max)) * hs(1 - psi/psi_max)
    return W

def W_quadric(theta, psi_max, par):
    psi = abs(90 * DEG - theta) << RAD
    W = (1 - (psi/psi_max)**(2)) * hs(1 - psi/psi_max)
    return W

def W_power_n(theta, psi_max, par):
    n = par
    psi = abs(90 * DEG - theta) << RAD
    W = (1 - (psi/psi_max)**(n)) * hs(1 - psi/psi_max)
    return W

def W_exp(theta, psi_max, par):
    n = par
    psi = abs(90 * DEG - theta) << RAD
    W = (2 - exp(psi/psi_max * ln(2))) * hs(1 - psi/psi_max)
    return W

def get_w_fun():
    W_fun = {
        'none': W_none,
        'const': W_const,
        'sqrt': W_sqrt,
        'line': W_line,
        'quadric': W_quadric,
        'power-n':W_power_n,
        'exp': W_exp,
    }
    return W_fun

def W_0(theta, theta_max, key_w, par_w):
    W_fun = get_w_fun()
    return W_fun[key_w](theta, theta_max, par_w)

def W_model(theta, theta_max, key_w, par_w, omega_kep, omega_rot):
    part = omega_rot/omega_kep
    W = W_0(theta, theta_max, key_w, par_w)
    W_m = W * (1 - part) + part
    return W_m / part

# def G_eff(theta, theta_max, key_w, par_w):
#     psi = abs(90 * DEG - theta) << RAD
#     W = W_0(theta, theta_max, key_w, par_w)
#     l_bl = 1
#     g_eff = l_bl * (1 - W**2) * hs(1 - psi/theta_max) + hs(psi/theta_max - 1)
#     return g_eff

# def Doppler(theta, theta_max, phi, V_k, V_r, key_w, par_w):
#     psi = abs(90 * DEG - theta) << RAD
#     chi = (phi - 180 * DEG) << RAD
#     W = W_0(theta, theta_max, key_w, par_w)
#     beta = V_k/C_LIGHT * W 
#     beta_s = sin(chi) * beta

#     doppler = sqrt((1 - beta_s)/(1 + beta_s))
#     D = doppler * hs(1 - psi/theta_max) + hs(psi/theta_max - 1)
#     return D

