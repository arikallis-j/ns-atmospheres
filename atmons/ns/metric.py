from .phisical import *

# for Neutron Star 

def chi_metric(R: Q[u.cm], M: Q[u.g]) -> Q[1]:
    """Calculate chi for the neutron star metric (II.eq.A.1))"""
    chi = G_GRAV * M / (R * C_LIGHT**2) #dze
    return chi << u.Unit()

def Omega_metric(R: Q[u.cm], M: Q[u.g], omega_rot: Q[u.s**(-1)]) -> Q[1]:
    """Calculate Omega for the neutron star metric (II.eq.A.1)"""
    Omega = omega_rot * sqrt(R**3 / (G_GRAV * M)) #epsm
    return Omega << u.Unit()

def q_c_metric(chi: Q[1], Omega: Q[1]) -> Q[1]:
    """Calculate q_c for the neutron star metric (II.eq.A.8)"""
    q_c = -0.11 * (Omega / chi)**2 #q
    return q_c << u.Unit()

def b_c_metric(chi: Q[1], Omega: Q[1]) -> Q[1]:
    """Calculate b_c for the neutron star metric (II.eq.A.9)"""
    b_c = 0.4454 * Omega**2 * chi #bet
    return b_c << u.Unit()

def i_bar_metric(chi: Q[1]) -> Q[1]:
    """Calculate i_bar for the neutron star metric (II.eq.A.11)"""
    i_m = sqrt(chi) * (1.136 - 2.53 * chi + 5.6 * chi**2)
    return i_m << u.Unit()

def g_0_metric(R: Q[u.cm], M: Q[u.g], chi: Q[1]) -> Q[u.cm / u.s**2]:
    """Calculate g_0 for the neutron star metric (II.eq.A.15)"""
    g0 = (G_GRAV * M / R**2) / sqrt(1.0 - 2.0 * chi)
    return g0 << u.cm / u.s**2


# for Neurton Star Radius

def R_metric(R_eq, phi, theta, NS):
    chi, Omega = NS.chi, NS.Omega 
    Omega2 = Omega**2

    a0 = -0.18 * Omega2 + 0.23 * chi * Omega2 - 0.05 * Omega2**2
    a2 = -0.39 * Omega2 + 0.29 * chi * Omega2 + 0.13 * Omega2**2
    a4 = 0.04 * Omega2 - 0.15 * chi * Omega2 + 0.07 * Omega2**2

    return R_eq * (1.0 + P_0(cos(theta)) * a0 + P_2(cos(theta)) * a2 + P_4(cos(theta)) * a4)

def dR_metric(R_eq, sin_th, cos_th, chi, Omega):
    cos_th = np.where(cos_th < 0.0, -cos_th, cos_th)
    Omega2 = Omega**2

    a2 = -0.39 * Omega2 + 0.29 * chi * Omega2 + 0.13 * Omega2**2
    a4 = 0.04 * Omega2 - 0.15 * chi * Omega2 + 0.07 * Omega2**2

    return R_eq * (a2 * dP_2(cos_th, sin_th) + a4 * dP_4(cos_th, sin_th))

def u_metric(R_theta, R_sch):
    return R_sch/R_theta

# for Neurton Star Point

def r_u_metric(R, cos_th, q_c, b_c, R_sch, err=1e-8):
    R_0 = R_sch / 2.0
    P2 = P_2(cos_th)

    u_mid, r_mid = np.zeros(cos_th.shape), np.zeros(cos_th.shape)
    r_i, r_bar = np.full(cos_th.shape, R)*R.unit, np.zeros(cos_th.shape)

    converged = np.full(cos_th.shape, False)
    while np.logical_not(converged.all()):
        r_mid = r_i
        u_mid = R_0 / r_mid

        nu_0 = ln((1.0 - u_mid/2.0) / (1.0 + u_mid/2.0))
        B_0 = (1.0 - u_mid/2.0) * (1.0 + u_mid/2.0)

        nu = nu_0 + (b_c/3.0 - q_c * P2) * u_mid**3
        B = B_0 + b_c * u_mid**2

        r_i = R / (exp(-nu) * B)

        errs = abs(r_i - r_mid) / r_i
        r_bar = np.where(np.logical_and(errs < err, np.logical_not(converged)), r_i, r_bar)
        converged = np.where(errs < err, True, converged)
    r_bar = r_i
    u_bar =  R_0 / r_bar # (II.eq.A.5)
   
    return r_bar, u_bar

def nu_B_dzeta_metric(cos_th, u_bar, q_c, b_c): 
    P2 = P_2(cos_th)

    # (II.eq.A.4)
    nu_0 = ln((1.0 - u_bar/2.0) / (1.0 + u_bar/2.0))
    B_0 = (1.0 - u_bar/2.0) * (1.0 + u_bar/2.0)
    dzeta_0 = ln(B_0)

    # (II.eq.A.7)
    nu = nu_0 + (b_c/3.0 - q_c * P2) * u_bar**3
    B = B_0 + b_c * u_bar**2
    zeta = dzeta_0 + b_c * u_bar**2 * (4.0 * P2 - 1.0) / 3.0

    return nu, B, zeta

def omega_bar_metric(r_bar, u_bar, J):
    omega_bar = (2.0 * G_GRAV * J) / (C_LIGHT**2 * r_bar**3) * (1.0 - 3.0 * u_bar) # (II.eq.A.10)
    return omega_bar

def beta_metric(R, sin_th, nu, omega_bar, omega_rot):
    beta = R * exp(-nu) * (omega_rot - omega_bar) * sin_th / C_LIGHT # (II.eq.B.37)
    return beta

def gamma_metric(beta):
    gamma = 1.0 / sqrt(1.0 - beta**2) # (II.eq.B.16)
    return gamma

def beta_ph_metric(R, sin_th, omega_bar, nu):
    beta_ph = R * omega_bar / C_LIGHT * exp(-nu) * sin_th # (II.eq.B.36)
    return beta_ph << u.Unit()

def g_metric(sin_th, cos_th, chi, Omega):
    cos_th = np.where(cos_th < 0.0, -cos_th, cos_th)
    Omega2 = Omega**2

    ce = 0.776 * chi - 0.791 # (II.eq.A.16)
    cp = 1.138 - 1.431 * chi # (II.eq.A.17)

    de = (2.431 * chi - 1.315) * chi
    dp = (0.653 - 2.864 * chi) * chi

    fe = -1.172 * chi
    fp = 0.975 * chi

    d60 = (13.47 - 27.13 * chi) * chi
    f60 = 1.69

    # approximation for the gravity from AlGendy & Morsink 2014
    ee = (ce + (de + fe * Omega2) * Omega2) * Omega2
    ep = (cp + (dp - d60 + (fp - f60) * Omega2) * Omega2) * Omega2
    ep = (cp + (dp - d60 + fp * Omega2) * Omega2) * Omega2

    g_th = 1.0 + ee * sin_th**2 + ep * cos_th**2 + ( d60 + f60 * Omega2 ) * Omega2**2 * cos_th

    return g_th

def grv_metric(theta, g_th, W, w_args, g_0):
    g = g_th * g_0 * (1 - W(theta, w_args)**2)
    return g

def f_theta(R, dR, nu, B, zeta):
    f =  B * exp(- zeta - nu) * dR/R # ??? (II.eq.B.25)
    return f

def eta_metric(f_th):
    sin_eta = f_th / sqrt(1 + f_th**2)
    cos_eta = 1 / sqrt(1 + f_th**2)
    return sin_eta, cos_eta

def dS_metric(th, ph, R, cos_eta, Surface):
    n_phi, n_theta = Surface.size
    theta = Surface.get_point(th, ph).theta
    if n_phi==0:
        return 0.0
    ph_min, ph_max = Surface.ranges[0]
    th_min, th_max = Surface.ranges[1]
    dphi = (ph_max - ph_min) / (n_phi - 1)
    dtheta = (th_max - th_min) / (n_theta - 1)
    if theta == PI/2 - dtheta/2:
        theta_0 = Surface.get_point(th - 1, ph).theta
        theta_1 = Surface.get_point(th, ph).theta
        dcos_th = (cos(theta_0) + cos(theta_1)) / 2.0
    elif theta == PI/2 + dtheta/2:
        theta_0 = Surface.get_point(th + 1, ph).theta
        theta_1 = Surface.get_point(th, ph).theta
        dcos_th = (cos(theta_0) + cos(theta_1)) / 2.0
    elif theta == PI/2:
        theta_0 = Surface.get_point(th - 1, ph).theta
        theta_1 = Surface.get_point(th + 1, ph).theta
        dcos_th = (abs(cos(theta_0)) + abs(cos(theta_1))) / 2.0
    elif th > 0 and th < (n_theta - 1):
        theta_0 = Surface.get_point(th - 1, ph).theta
        theta_1 = Surface.get_point(th + 1, ph).theta
        dcos_th = (cos(theta_0) - cos(theta_1)) / 2.0
    elif th == 0:
        theta_0 = Surface.get_point(th, ph).theta
        theta_1 = Surface.get_point(th + 1, ph).theta
        dcos_th = 1.0 - (cos(theta_0) + cos(theta_1)) / 2.0
    elif th == (n_theta - 1):
        theta_0 = Surface.get_point(th - 1, ph).theta
        theta_1 = Surface.get_point(th, ph).theta
        dcos_th = 1.0 - (abs(cos(theta_0) + cos(theta_1))) / 2.0

    dS = 1 / cos_eta * R**2 * abs(dcos_th) * dphi
    return dS


def dS_metric_1(theta, cos_eta, R, N_ph, N_th, ph_range, th_range):
    dcos_th = np.zeros(theta.shape)
    if N_ph==0:
        return dcos_th
    
    ph_min, ph_max = ph_range
    th_min, th_max = th_range
    dphi = (ph_max - ph_min) / (N_ph - 1)
    dtheta = (th_max - th_min) / (N_th - 1)

    theta_0 = theta
    theta_p = theta + dtheta
    theta_m = theta - dtheta

    dcos_th = np.where(theta == PI/2*RAD - dtheta/2, (cos(theta_m) + cos(theta_0)) / 2.0 , dcos_th)
    dcos_th = np.where(theta == PI/2*RAD + dtheta/2, (cos(theta_0) + cos(theta_p)) / 2.0, dcos_th)
    dcos_th = np.where(theta == PI/2*RAD, (abs(cos(theta_m)) + abs(cos(theta_p))) / 2.0, dcos_th)
    dcos_th = np.where(np.logical_and(theta > th_min, theta < th_max), (cos(theta_m) - cos(theta_p)) / 2.0, dcos_th)
    dcos_th = np.where(theta == th_min, 1.0 - (cos(theta_0) + cos(theta_p)) / 2.0, dcos_th)
    dcos_th = np.where(theta == th_max, 1.0 - (abs(cos(theta_0) + cos(theta_m))) / 2.0, dcos_th)

    dS = 1 / cos_eta * R**2 * abs(dcos_th) * dphi
    return dS / RAD
