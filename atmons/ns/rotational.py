from .math import *

def psi_rot(sin_th: float, cos_th: float, cos_ph: float, sin_i: float, cos_i: float) -> tuple[float, float]:
    """Calculate sin and cos of psi angle."""
    # (II.eq.B.3)
    cos_psi = cos_i * cos_th + sin_i * sin_th * cos_ph
    sin_psi = sqrt(1.0 - cos_psi**2)
    return sin_psi, cos_psi

def G_yu_rot(cos_psi: float, u_th: float) -> float:
    """Calculate G_yu based on cos_psi and u_th."""
    y = 1.0 - cos_psi # (II.eq.B.7-8)
    # (II.eq.B.8)
    G_yu = 1.0 + (u_th * y)**2 / 112.0 - EULER * u_th * y / 100.0 * (ln(1.0 - 0.5 * y) + 0.5 * y) 
    return G_yu
        
def D_rot(cos_psi: float, u_th: float) -> float:
    """Calculate D based on cos_psi and u_th."""
    y = 1.0 - cos_psi # (II.eq.B.7-8)
    # dxdy (II.eq.B.9)
    D = 1.0 + 3.0 * (u_th * y)**2 / 112.0 - EULER * u_th * y / 100.0 * (2.0 * ln(1.0 - 0.5 * y) + y * (1.0 - 0.75 * y) / (1.0 - 0.5 * y))
    return D

def alpha_rot(cos_psi: float, u_th: float, G_yu: float) -> tuple[float, float]:
    """Calculate sin and cos of alpha angle."""
    # (II.eq.B.7)
    cos_a = 1.0 - (1.0 - cos_psi) * (1.0 - u_th) * G_yu
    sin_a = sqrt(1.0 - cos_a**2)
    return sin_a, cos_a

        
def chi_rot(sin_th: float, cos_th: float, sin_psi: float, cos_psi: float, cos_i: float) -> float:
    """Calculate cos of chi angle."""
    # (II.eq.B.30)
    cos_chi = (cos_i - cos_th * cos_psi) / (sin_th * sin_psi)
    return cos_chi

def sigma_rot(sin_eta: float, cos_eta: float, sin_a: float, cos_a: float, cos_chi: float, cos_th: float) -> float:
    """Calculate cos of sigma angle."""
    # (II.eq.B.29)
    cos_sig = cos_eta * cos_a + sin_eta * sin_a * cos_chi * cos_th / abs(cos_th)
    return cos_sig

def xi_rot(sin_a: float, sin_psi: float, sin_ph: float, sin_i: float) -> float:
    """Calculate cos of xi angle."""
    # (II.eq.B.17)
    cos_xi = -1.0 * sin_a * sin_i * sin_ph / sin_psi
    return cos_xi

def delta_rot(beta: float, gamma: float, cos_xi: float) -> float:
    """Calculate delta."""
    # (II.eq.B.14)
    delta = 1.0 / (gamma * (1.0 - beta * cos_xi))
    return delta

def sigma_1_rot(cos_sig: float, delta: float) -> float:
    """Calculate cos of sigma_1 angle."""
    # (II.eq.B.32)
    cos_sig_1 = cos_sig * delta
    return cos_sig_1