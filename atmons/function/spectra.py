"""
Description of spectral relationships
"""
from ..const import *
from .math import *
import json, yaml

# Base spectra grid
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

# Deluted spectra coefficients
def w_b(chem: str, fc_key: int) -> list[list[Q[DIMLESS]]]:
    """Dilution factor (by the previous fitting)"""
    with open(f"spectra/fcol_{chem}.json") as f:
        FCOL = np.array(json.load(f))

    if fc_key == 1:
        f_c = FCOL[:,:,4]
        w_fc = FCOL[:,:, 9]
    else:
        f_c = FCOL[:,:,5]
        w_fc = FCOL[:,:, 10]

    w_b = w_fc * f_c ** (-4)

    return w_b << u.Unit()

def T_c(chem: str, fc_key: int) -> list[list[Q[KEV]]]:
    """Color temperature (by the previous fitting)"""
    with open(f"spectra/fcol_{chem}.json") as f:
        FCOL = np.array(json.load(f))

    if fc_key == 1:
        T_eff = FCOL[:,:,3]
        f_c = FCOL[:,:,4]
    else:
        T_eff = FCOL[:,:,3]
        f_c = FCOL[:,:,5]

    T_c = T_eff * f_c

    return T_c << KEV

# Model spectra
def B_model(chem: str):
    with open(f"spectra/spec_{chem}.json") as f:
        SPEC = np.array(json.load(f))

    N_int = SPEC
    B_int = np.zeros(N_int.shape)

    for k in range(len(ENERGY)):
        E_f = ENERGY_LONG[k]
        E_l = ENERGY_LONG[k+1]
        E_shark =  exp10((log(E_f) + log(E_l))/2)
        B_int[:,:,k] =  N_int[:,:,k] / (E_l - E_f) * E_shark
    
    B_int = B_int / (1.0 * ERG / KEV) << u.Unit()
    B_int = B_int / PI 
    B_int = B_int * ERG / (SEC * CM**2 * KEV)

    return B_int << ERG / (SEC * CM**2 * KEV)