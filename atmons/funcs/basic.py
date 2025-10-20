"""
Description of basic relationships
"""
from ..consts import *
from .math import *
 
# Mass and Radius in CGS system
def r_min(m_ns: Q[M_SUN]) -> Q[KM]:
    """Minimum neutron star radius (from ???)"""
    r_m = 2.95 * KM * (1.5 * m_ns/M_SUN)
    return r_m << u.km

def R_NS(r_ns: Q[KM], m_ns: Q[M_SUN]) -> Q[CM]:
    """Neutron star radius"""
    r_m = r_min(m_ns)
    if r_ns <= r_m:
        print("Causility error")
        print(f"r_ns < r_min = {r_ns} < {r_m:.2f}")
        r_ns = 0
    return r_ns << CM

def M_NS(m_ns: Q[M_SUN]) -> Q[GRAM]:
    """Neutron star mass"""
    return m_ns << GRAM

# Critical frequency of rotation
def omega(nu: Q[HZ]) -> Q[SEC**(-1)]:
    """Cyclical frequency (by definition)"""
    omg =  2.0 * PI * nu
    return omg << SEC**(-1)

def nu_crit(r_ns: Q[KM], m_ns: Q[M_SUN]) -> Q[HZ]:
    """Maximum possible rotation frequency (II.eq.1)"""
    nu = 1278.0 * HZ * (10.0*KM/r_ns)**1.5 * sqrt(m_ns/(1.4*M_SUN))
    return nu << HZ

def nu_relative(nu_rot: Q[HZ], nu_crit: Q[HZ]) -> Q[DIMLESS]:
    """Relative rotation frequency (by definition)"""
    nu_bar = nu_rot / nu_crit
    return nu_bar << DIMLESS

# Velocity of rotation 
def omega_rot(V: Q[CM / SEC], R: Q[CM]) -> Q[SEC**(-1)]:
    """Angular frequency of rotation (by definition)"""
    omega = V / R
    return omega << SEC**(-1)

def V_rot(omega: Q[SEC**(-1)], R: Q[CM]) -> Q[CM / SEC]:
    """Linear velocity of rotation (by definition)"""
    V = omega * R
    return V << CM / SEC

def V_kep(g: Q[CM/SEC**2], R: Q[CM]) -> Q[CM / SEC]:
    """Linear velocity of keplerian movement"""
    V = sqrt(g * R)
    return V << CM / SEC

# Mass and Radius in relative units
def r_eq(r_ns: Q[KM], m_ns: Q[M_SUN], v_rot: Q[HZ]) -> Q[KM]:
    """Equatorial radius of Neutron Star (II.eq.2)"""
    v_crit = nu_crit(r_ns, m_ns)
    nu_rel = nu_relative(v_rot, v_crit)
    r_1 = 0.025 
    r_2 = 0.07 * (m_ns/(1.4*M_SUN))**1.5

    r = r_ns * (0.9766 - r_1 / (nu_rel - 1.07) + r_2 * nu_rel**2)

    return r << KM
    
def m_cor(m_ns: Q[M_SUN], r_ns: Q[KM], v_rot: Q[HZ]) -> Q[M_SUN]:
    """Corrected mass of Neutron Star (II.eq.3)"""
    v_crit = nu_crit(r_ns, m_ns)
    nu_rel = nu_relative(v_rot, v_crit)

    a_1 = 0.001 * (m_ns/(1.4*M_SUN))**1.5
    a_2 = 10.0 * a_1
    a_0 = 1.0 - a_1/1.1

    m = m_ns * (a_0 - a_1 / (nu_rel - 1.1) + a_2 * nu_rel**2)

    return m << M_SUN

# Inverse Mass and Radius in relative units
def rel_r_eq(r_ns: Q[KM], m_ns: Q[M_SUN], v_rot: Q[HZ]) -> Q[KM]:   
    """Inverse of corrected mass of Neutron Star (II.eq.3)"""
    r = inverse(r_eq, r_ns, (m_ns, v_rot), base=1.0, pw=5)
    return r

def rel_m_cor(m_ns: Q[M_SUN], r_ns: Q[KM], v_rot: Q[HZ]) -> Q[M_SUN]:
    """Inverse of equatorial radius of Neutron Star (II.eq.2)"""
    m = inverse(m_cor, m_ns, (r_ns, v_rot), base=1.0, pw=5)
    return m

def max_latitude(g_mod, g_base, theta, psi_star):
    """Maximum latitude for the model"""
    l_arr = len(g_mod[0, ::]) // 2
    g_base_th = g_base[0,:l_arr-1:]
    g_mod_th = g_mod[0,:l_arr-1:]
    theta_th = theta[0,:l_arr-1:]
    th_max = 90*DEG - psi_star << RAD
    for k in range(1,len(g_mod_th)):
        if (g_mod_th[k] - g_base_th[k])*(g_mod_th[k-1] - g_base_th[k-1]) < 0:
            th_max = (theta_th[k] + theta_th[k-1])/2
    psi_max = 90*DEG - th_max << DEG
    return psi_max