from .const import *
from .math import *

def r_min(m_ns):
    return 2.95 * m_ns * 1.5

def R_NS(r_ns, m_ns):
    r_m = r_min(m_ns)
    if r_ns <= r_m:
        print("Causility error")
        print(f"r_ns < r_min = {r_ns} < {round(r_m,2)}")
        return 0
    else:
        return r_ns * KM

def M_NS(m_ns):
    return m_ns * M_SUN

def I_NS(i, R, M):
    return i * M * R**2  #tin (II.eq.A.10-A.11)

def J_NS(I, omega):
    return I * omega # angular momentum (II.eq.A.10-11)

def R_sch(M):
    return 2.0 * G_GRAV * M / C_LIGHT**2

def zsch(R, R_s):
    return 1.0 / sqrt(1.0 - R_s / R)

def g(R, M, zst):
    return (G_GRAV * M / R**2) * zst

def Surf(R, zst):
    return (R * zst)**2

def X_hyd(chem):
    if chem=="he":
        return 0.0
    elif chem=="s1" or chem=="s001":
        return 0.7374
    else:
        return 0.0
    
def sigma_T(x_hyd):
    return 0.2 * (1.0 + x_hyd)

def Flux_edd(grv, sigma):
    return grv * C_LIGHT / sigma

def Flux_SB(T):
    return SIGMA_SB * T**4

def T_SB(Flux):
    return (Flux/SIGMA_SB)**(1/4)

def T_obs(Flux, zst):
    return (Flux/SIGMA_SB)**(1/4) / zst

def Theta(T):
    return K_B * T

def Epsilon(T):
    return (K_B * T) * KEV

def Lumen(Flux, R):
    return 4.0 * PI * R**2 * Flux

def Lumen_obs(Lum, zst):
    return Lum / zst**2

def omega(nu):
    return 2.0 * PI * nu

def cos_i(i):
    return cos(i)

def sin_i(i):
    return sin(i)

def nu_crit(r_ns, m_ns):
    return 1278.0 * (10.0/r_ns)**1.5 * sqrt(m_ns/1.4)

def nu_relative(nu_rot, nu_crit):
    return nu_rot / nu_crit 

def r_eq(r_ns, m_ns, v_rot, define=False, r=0.0):
    if define:
        return r
    v_crit = nu_crit(r_ns, m_ns)
    nu_rel = nu_relative(v_rot, v_crit)

    r_1 = 0.025 
    r_2 = 0.07 * (m_ns/1.4)**1.5
    return r_ns * (0.9766 - r_1 / (nu_rel - 1.07) + r_2 * nu_rel**2)
    
def m_cor(m_ns, r_ns, v_rot, define=False, m=0.0):
    if define:
        return m
    v_crit = nu_crit(r_ns, m_ns)
    nu_rel = nu_relative(v_rot, v_crit)

    a_1 = 0.001 * (m_ns/1.4)**1.5
    a_2 = 10.0 * a_1
    a_0 = 1.0 - a_1/1.1

    return m_ns * (a_0 - a_1 / (nu_rel - 1.1) + a_2 * nu_rel**2)

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

def w_b(chem, fc_key):
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

    with open(f'ns_atm/spectra/{name}', 'r') as f:
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

    return w_b

def T_c(chem, fc_key):
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

    with open(f'ns_atm/spectra/{name}', 'r') as f:
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

    return T_c

def T_flux_eff(key_l, flux_NS, T_eff_NS, Flux_edd_NS):
    if key_l == 1:
        T_eff = T_SB(flux_NS * Flux_edd_NS)
        flux = flux_NS
    else:
        T_eff = T_eff_NS 
        flux = Flux_SB(T_eff) / Flux_edd_NS
        if flux_NS >= FLUX_REL[N_MODEL-1]:
            # N_l = iif - 1
            print(f"""incorrectly flux
            flux_i  = {flux_NS}
            flux_rel[N_l] = {FLUX_REL[N_MODEL-1]})
            flux_i ≥ flux_rel[N_l] = {flux_NS}≥{FLUX_REL[N_MODEL-1]}) = true
            """)
    return flux, T_eff

def wwf_tcf_T(T_c, w_b, flux_i, log_g):
    T_c_ig, w_ig = [], []
    T_g, w_g = [], []

    for ig in range(0,9):
        T_c_ig, w_ig = [], []
        # Внутренний цикл для заполнения T_c_ig и w_ig
        for i in range(0, N_MODEL):
            T_c_ig.append(T_c[i][ig])
            w_ig.append(w_b[i][ig])
        T_c_ig = np.array(T_c_ig)
        w_ig = np.array(w_ig)
        # Вызов функции интерполяции для T_c

        T_i = map1(FLUX_REL, T_c_ig, N_MODEL, flux_i)  # Предполагаем, что map1 возвращает пару значений
        T_g.append(T_i)
            
        # Вызов функции интерполяции для w
        w_i = map1(FLUX_REL, w_ig, N_MODEL, flux_i)  # Аналогичное предположение для вывода функции map1
        w_g.append(w_i)

    T_f = map1(LOG_G,T_g,9,log_g)
    w_f = map1(LOG_G,w_g,9,log_g)  
        
    wwf_T = w_f
    tcf_T = T_f
    return wwf_T, tcf_T
