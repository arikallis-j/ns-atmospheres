"""
Description of model relationships
"""
from ..const import *
from .math import *

def W_none(rho, t, xi, par_w):
    W = np.zeros(t.shape)
    return W

def W_const(rho, t, xi, par_w):
    W = hs(1 - t)
    return W

def V_const(rho, t, xi, par_w):
    V = hs(1 - t)
    W = V / rho
    return W

def Vk_const(rho, t, xi, par_w):
    V = hs(1 - t) 
    W = V / np.sqrt(rho)
    return W

def W_line(rho, t, xi, par_w):
    W = (1 - t) * hs(1 - t)
    return W

def V_line(rho, t, xi, par_w):
    V = (1 - t) * hs(1 - t)
    W = V / rho
    return W

def Vk_line(rho, t, xi, par_w):
    V = (1 - t) * hs(1 - t)
    W = V / np.sqrt(rho)
    return W

def W_power_n(rho, t, xi, par_w):
    n = par_w
    W = (1 - t**(n)) * hs(1 - t)
    return W

def V_power_n(rho, t, xi, par_w):
    n = par_w
    V = (1 - t**(n)) * hs(1 - t)
    W = V / rho
    return W

def Vk_power_n(rho, t, xi, par_w):
    n = par_w
    V = (1 - t**(n)) * hs(1 - t)
    W = V / np.square(rho)
    return W

def get_x_fun():
    X_fun = {
        'base': W_none,
        'const': W_const,
        'vconst': V_const,
        'vkconst': Vk_const,
        'line': W_line,
        'vline': V_line,
        'vkline': Vk_line,
        'power-n':W_power_n,
        'vpower-n':V_power_n,
        'vkpower-n':Vk_power_n,
        # 'sqrt': W_sqrt,
        # 'quadric': W_quadric,
        # 'exp': W_exp,
    }
    return X_fun


def x_0(rho, t, xi, key_w, par_w):
    x_fun = get_x_fun()
    if key_w[:6:] == 'power-':
        par_w = float(key_w[-1:5:-1][::-1])
        key_w = 'power-n'
    elif key_w[:7:] == 'vpower-':
        par_w = float(key_w[-1:6:-1][::-1])
        key_w = 'vpower-n'
    elif key_w[:8:] == 'vkpower-':
        par_w = float(key_w[-1:7:-1][::-1])
        key_w = 'vkpower-n'
        
    return x_fun[key_w](rho, t, xi, par_w)


def W_model(r, r_eq, theta, psi_max, key_w, par_w, omega_kep, omega_rot):
    xi = omega_kep/omega_rot
    psi = abs(90 * DEG - theta) << RAD
    t = psi / psi_max
    rho = r/r_eq
    x = x_0(rho, t, xi, key_w, par_w)
    W = x * (xi - 1) + 1
    return W 

# additional functions

# def W_sqrt(r, r_eq, theta, psi_max, par)):
#     psi = abs(90 * DEG - theta) << RAD
#     W = (1 - (psi/psi_max)**(1/2)) * hs(1 - psi/psi_max)
#     return W

# def W_line(r, r_eq, theta, psi_max, par)):
#     psi = abs(90 * DEG - theta) << RAD
#     W = (1 - (psi/psi_max)) * hs(1 - psi/psi_max)
#     return W

# def W_quadric(r, r_eq, theta, psi_max, par)):
#     psi = abs(90 * DEG - theta) << RAD
#     W = (1 - (psi/psi_max)**(2)) * hs(1 - psi/psi_max)
#     return W

# def W_power_n(r, r_eq, theta, psi_max, par)):
#     n = par
#     psi = abs(90 * DEG - theta) << RAD
#     W = (1 - (psi/psi_max)**(n)) * hs(1 - psi/psi_max)
#     return W

# def W_exp(r, r_eq, theta, psi_max, par)):
#     n = par
#     psi = abs(90 * DEG - theta) << RAD
#     W = (2 - exp(psi/psi_max * ln(2))) * hs(1 - psi/psi_max)
#     return W