from .const import *
from .math import *
 
def r_min(m_ns: Q[u.Msun]) -> Q[u.km]:
    """Minimum neutron star radius (from ???)"""
    r_m = 2.95 * KM * (1.5 * m_ns/M_SUN)
    return r_m << u.km

def R_NS(r_ns: Q[u.km], m_ns: Q[u.Msun]) -> Q[u.cm]:
    """Neutron star radius"""
    r_m = r_min(m_ns)
    if r_ns <= r_m:
        print("Causility error")
        print(f"r_ns < r_min = {r_ns} < {r_m:.2f}")
        r_ns = 0
    return r_ns << u.cm

def M_NS(m_ns: Q[u.Msun]) -> Q[u.g]:
    """Neutron star mass"""
    return m_ns << u.g

def I_NS(i: Q[1], R: Q[u.cm], M: Q[u.g]) -> Q[u.cm**2 * u.g]:
    """Neutron star moment of inertia (II.eq.A.10-A.11)"""
    I = i * M * R**2 
    return I << u.cm**2 * u.g

def J_NS(I:Q[u.cm**2 * u.g], omega: Q[u.s**(-1)]) -> Q[u.cm**2 * u.g / u.s]: 
    """Neutron star angular momentum (II.eq.A.10-11)"""
    J = I * omega
    return J << u.cm**2 * u.g / u.s

def R_sch(M: Q[u.g]) -> Q[u.cm]:
    """Schwarzschild radius (by definition)"""
    r_sch = 2.0 * G_GRAV * M / C_LIGHT**2 
    return r_sch << u.cm

def zsch(R: Q[u.cm], R_s: Q[u.cm]) -> Q[1]: 
    """Schwarzschild red shift plus one (I.eq.3)"""
    z = 1.0 / sqrt(1.0 - R_s / R)
    return z << u.Unit()

def g(R: Q[u.cm], M: Q[u.g], zst: Q[1]) -> Q[u.cm / u.s**2]:
    """Acceleration of freefall (by relativity's definition)"""
    grv = (G_GRAV * M / R**2) * zst
    return grv << u.cm / u.s**2

def Surf(R: Q[u.cm], zst: Q[1]) -> Q[u.cm**2]:
    """Surface area (by relativity's definition)"""
    surf = (R * zst)**2
    return surf << u.cm**2

def X_hyd(chem: str) -> Q[1]:
    """Hydrogen mass fraction (by observation)"""
    if chem=="he":
        x_hyd =  0.0
    elif chem=="s1" or chem=="s001":
        x_hyd = 0.7374
    else:
        x_hyd = 0.0 
    return x_hyd << u.Unit()
    
def kappa_e(x_hyd: Q[1]) -> Q[u.cm**2 / u.g]:
    """Coherent Thomson electron scattering opacity (I.eq.8)"""
    kappa = 0.2 * (1.0 + x_hyd) * u.cm**2 / u.g
    return kappa << u.cm**2 / u.g

def Flux_edd(grv: Q[u.cm/u.s**2], kappa: Q[u.cm**2/u.g]) -> Q[u.erg/(u.s*u.cm**2)]:
    """Bolometric Eddington flux (I.eq.9)"""
    flux_edd = grv * C_LIGHT / kappa
    return flux_edd << u.erg / (u.s * u.cm**2)

def Flux_SB(T: Q[u.K]) -> Q[u.erg/(u.s * u.cm**2)]:
    """Flux of black body (by definition)"""
    flux_SB = SIGMA_SB * T**4
    return flux_SB << u.erg / (u.s * u.cm**2)

def T_SB(Flux: Q[u.erg/(u.s * u.cm**2)]) -> Q[u.K]:
    """Temperature of black body (by definition)"""
    T = (Flux/SIGMA_SB)**(1/4)
    return T << u.K

def T_obs(Flux: Q[u.erg/(u.s*u.cm**2)], zst: Q[1]) -> Q[u.K]:
    """Temperature of observed black body (by relativity's definition)"""
    T = (Flux/SIGMA_SB)**(1/4) / zst
    return T << u.K

def Theta(T: Q[u.K]) -> Q[u.erg]:
    """Temperature in ergs (by definition)"""
    theta = K_B * T
    return theta << u.erg 

def Epsilon(T: Q[u.K]) -> Q[u.keV]:
    """Temperature in keV (by definition)"""
    epsilon = K_B * T
    return epsilon << u.keV 

def Lumen(Flux: Q[u.erg/(u.s*u.cm**2)], R: Q[u.cm]) -> Q[u.erg/u.s]: 
    """Luminosity (by definition)"""
    lumen = 4.0 * PI * R**2 * Flux
    return lumen << u.erg / u.s

def Lumen_obs(Flux: Q[u.erg/(u.s*u.cm**2)], R: Q[u.cm], zst: Q[1]) -> Q[u.erg/u.s]:
    """Observed luminosity (by relativity's definition)"""
    lumen = 4.0 * PI * R**2 * Flux / zst**2
    return lumen << u.erg / u.s

def omega(nu: Q[u.Hz]) -> Q[u.s**(-1)]:
    """Cyclical frequency (by definition)"""
    omg =  2.0 * PI * nu
    return omg << u.s**(-1)

def nu_crit(r_ns: Q[u.km], m_ns: Q[u.Msun]) -> Q[u.Hz]:
    """Maximum possible rotation frequency (II.eq.1)"""
    nu = 1278.0 * u.Hz * (10.0*KM/r_ns)**1.5 * sqrt(m_ns/(1.4*M_SUN))
    return nu << u.Hz

def nu_relative(nu_rot: Q[u.Hz], nu_crit: Q[u.Hz]) -> Q[0]:
    """Relative rotation frequency (by definition)"""
    nu_bar = nu_rot / nu_crit
    return nu_bar << u.Unit()

def r_eq(r_ns: Q[u.km], m_ns: Q[u.Msun], v_rot: Q[u.Hz], define=False, r: Q[u.km] = 0*u.km) -> Q[u.km]:
    """Equatorial radius of Neutron Star (II.eq.2)"""
    v_crit = nu_crit(r_ns, m_ns)
    nu_rel = nu_relative(v_rot, v_crit)
    r_1 = 0.025 
    r_2 = 0.07 * (m_ns/(1.4*M_SUN))**1.5

    if not(define):
        r = r_ns * (0.9766 - r_1 / (nu_rel - 1.07) + r_2 * nu_rel**2)

    return r << u.km
    
def m_cor(m_ns: Q[u.Msun], r_ns: Q[u.km], v_rot: Q[u.Hz], define=False, m: Q[u.Msun] = 0*u.Msun) -> Q[u.Msun]:
    """Corrected mass of Neutron Star (II.eq.3)"""
    v_crit = nu_crit(r_ns, m_ns)
    nu_rel = nu_relative(v_rot, v_crit)

    a_1 = 0.001 * (m_ns/(1.4*M_SUN))**1.5
    a_2 = 10.0 * a_1
    a_0 = 1.0 - a_1/1.1

    if not(define):
        m = m_ns * (a_0 - a_1 / (nu_rel - 1.1) + a_2 * nu_rel**2)

    return m << u.Msun

def E_base(N, range_E):
    E  = np.array([0.0]*N)
    dE = np.array([0.0]*N)
    E_min, E_max = range_E
    
    # Заполнение сетки энергий фотонов
    E[0] = E_min #keV
    E[N-1] = E_max #keV
    dlog_E = (log(E[N-1]) - log(E[0])) / (N-1)

    for nu in range(1, N-1):
        log_E = log(E[0]) + dlog_E * nu
        E[nu] = exp10(log_E)

    # Заполнение сетки изменений энергий 
    dE[0] = (E[1] - E[0]) / 2.0
    dE[N-1] = (E[N-1] - E[N-2]) / 2.0

    for nu in range(1, N-1):
        dE[nu] = (E[nu+1] - E[nu-1]) / 2.0

    return E, dE

def E_doppler(E, dE, B_d, E_d, dE_d):
    grid_a, grid_b, e = E.shape
    grid_ad, grid_bd, ed = E_d.shape
    E1 = np.zeros((grid_a, grid_b, e+1)) * E.unit
    eps = np.zeros((grid_a, grid_b, e)) * B_d.unit * dE_d.unit
    B = np.zeros((grid_a, grid_b, e)) * B_d.unit 

    E_d1 = np.zeros((grid_ad, grid_bd, ed+1)) * E_d.unit
    eps_d = np.zeros((grid_ad, grid_bd, ed)) * B_d.unit * dE_d.unit
    E1[:,:,0] = E[:,:,0]
    E_d1[:,:,0] = E_d[:,:,0]

    E1[:,:,1:] = E1[:,:,:-1] + dE
    E_d1[:,:,1:] = E_d1[:,:,:-1]  + dE_d
    eps_d = B_d * dE_d 

    left_border = (E_d1[:,:,:-1] + dE_d) >= E1[:,:,0:1]
    right_border = (E_d1[:,:,:-1] + dE_d) <= E1[:,:,1:2]
    border = left_border & right_border
    deps_mk = np.where(border, eps_d * ((dE_d - (E1[:,:,0:1] - E_d1[:,:,:-1]))/(dE_d)), eps_d * 0)
    
    left_border_in = (E_d1[:,:,:-1]) > E1[:,:,:-1]
    right_border_in = (E_d1[:,:,:-1]) <= E1[:,:,1:]
    border_in = left_border_in & right_border_in

    deps_k = np.where(border_in, eps_d * ((dE)/(dE_d) - (E_d1[:,:,:-1] - E1[:,:,:-1])/(dE_d)), eps_d * 0)
    deps_mk = eps_d - deps_k 
    eps += deps_k + deps_mk

    B = eps/dE

    return B

def w_b(chem: str, fc_key: int) -> list[list[Q[1]]]:
    """Dilution factor (by the previous fitting)"""
    w_b = [[0.0]*9]
    for k in range(N_MODEL - 1):
        y = [[0.0]*9]
        w_b = w_b + y

    if chem=="s001":
        name = "fcol_S001.dat"
    elif chem=="s1":
        name = "fcol_S1.dat"
    elif chem=="he":
        name = "fcol_He.dat"
    else:
        print("There is not this chemical composion")
        return w_b

    with open(f'spectra/{name}', 'r') as f:
        count, j, ig = 0, 0, 0
        for line in f:
            j = count % N_MODEL
            ig = count // N_MODEL 
            if fc_key == 1:
                a1, a2, gr, T_eff, f_c, a4, a5, a6, a7, w_fc, a8, a9 = line.split()
            else:
                a1, a2, gr, T_eff, a3, f_c, a4, a5, a6, a7, w_fc, a8 = line.split()
            w_fc = float(w_fc)
            f_c = float(f_c)
            T_eff = float(T_eff)
            w_b[j][ig] = w_fc * f_c ** (-4)
            count += 1

    return w_b << u.Unit()

def T_c(chem: str, fc_key: int) -> list[list[Q[u.keV]]]:
    """Color temperature (by the previous fitting)"""
    T_c = [[0.0]*9]
    for k in range(N_MODEL - 1):
        y = [[0.0]*9]
        T_c = T_c + y

    if chem=="s001":
        name = "fcol_S001.dat"
    elif chem=="s1":
        name = "fcol_S1.dat"
    elif chem=="he":
        name = "fcol_He.dat"
    else:
        print("There is not this chemical composion")
        return T_c

    with open(f'spectra/{name}', 'r') as f:
        count, j, ig = 0, 0, 0
        for line in f:
            ig = count // (N_MODEL)
            j = count % N_MODEL
            if fc_key == 1:
                a1, a2, gr, T_eff, f_c, a4, a5, a6, a7, w_fc, a8, a9 = line.split()
            else:
                a1, a2, gr, T_eff, a3, f_c, a4, a5, a6, a7, w_fc, a8 = line.split()
            w_fc = float(w_fc)
            f_c = float(f_c)
            T_eff = float(T_eff)
            T_c[j][ig] = T_eff * f_c
            count += 1

    return T_c << u.keV

def T_flux_eff(key_l, flux_NS, T_eff_NS, Flux_edd_NS, N_model):
    if key_l == 1:
        T_eff = T_SB(flux_NS * Flux_edd_NS)
        flux = np.full(Flux_edd_NS.shape, flux_NS) * u.Unit()
    else:
        T_eff = T_eff_NS 
        flux = Flux_SB(T_eff) / Flux_edd_NS
        if (flux >= FLUX_REL[N_model-1]).any():
            N_model -= 1
            print("incorrectly flux")
            # print(f"""incorrectly flux
            # flux_i  = {flux_NS}
            # flux_rel[N_l] = {FLUX_REL[N_model-1]})
            # flux_i ≥ flux_rel[N_l] = {flux_NS}≥{FLUX_REL[N_model-1]}) = true
            # """)
    return flux << u.Unit(), T_eff, N_model

def wwf_tcf_T(T_c, w_b, flux_i, log_g):
    u_T, u_w = T_c.unit, w_b.unit
    T_c = np.array(T_c)
    w_b = np.array(w_b)

    T_c = np.tile(T_c[:, :, None, None], (1, 1, *log_g.shape))
    w_b = np.tile(w_b[:, :, None, None], (1, 1, *log_g.shape))

    T_g = map2(FLUX_REL, T_c, flux_i)
    w_g = map2(FLUX_REL, w_b, flux_i)

    T_g = np.tile(T_g[:, None, :, :], (1, 1, 1, 1))
    w_g = np.tile(w_g[:, None, :, :], (1, 1, 1, 1))

    T_f = map2(LOG_G, T_g, log_g)
    w_f = map2(LOG_G, w_g, log_g)

    T_f = T_f.reshape(log_g.shape)
    w_f = w_f.reshape(log_g.shape)

    wwf_T = w_f
    tcf_T = T_f
    return wwf_T << u_w, tcf_T << u_T

def B_inter(flux, log_g, E, chem):
    N_int = np.array(SPEC[chem]) 
    B_int = np.zeros(N_int.shape)

    for k in range(len(ENERGY)):
        E_f = ENERGY_LONG[k]
        E_l = ENERGY_LONG[k+1]
        E_shark =  exp10((log(E_f) + log(E_l))/2)
        B_int[:,:,k] =  N_int[:,:,k] / (E_l - E_f) * E_shark

    B_int = np.tile(B_int[:, :, :, None, None], (1, 1, 1, *flux.shape))
    # B_int = np.tile(B_int[:, :, :, None, None, None], (1, 1, 1, *E.shape))
    # flux = np.tile(flux[None, :, :], (B_int.shape[2], 1, 1))
    # log_g = np.tile(log_g[None, :, :], (B_int.shape[2], 1, 1))

    #print(len(FLUX_REL_SHORT), B_int.shape, flux.shape)
    
    B_lum = map2(FLUX_REL_SHORT, B_int, flux)

    #print(len(LOG_G), B_lum.shape, log_g.shape)
    
    B_lg = map2(LOG_G, B_lum, log_g)    
    
    B_lg = np.tile(B_lg[:, None, :, :, None], (1, 1, 1, 1, E.shape[2]))
    
    #print(len(ENERGY), B_lg.shape, E.shape)

    B_en = map2(ENERGY, B_lg, E)
    # B_en = B_en.transpose(0, 3, 1, 2)
    # B_en = B_en.reshape(E.shape)

    # B_en = B_lg[::10, :, :]
    # B_en = B_en.transpose(2, 0, 1)
    B_en = B_en.reshape(E.shape)

    B_en = B_en / (1.0 * u.erg / u.keV) << u.Unit()
    B_en = B_en / PI 
    
    return B_en * u.erg / (u.s * u.cm**2 * u.keV)

