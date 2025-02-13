from ..ns.classes import *

def lfc(cfg, grid):
    NS = NeurtonStar(cfg, grid)
    lum, fc = [], []
    for shot in NS.burst:
        lum.append(shot.lum)
        fc.append(shot.fc)
    return lum, fc

def lw(cfg, grid):
    NS = NeurtonStar(cfg, grid)
    lum, w = [], []
    for shot in NS.burst:
        lum.append(shot.lum)
        w.append(shot.w)
    return lum, w

def spectra(config, grid):
    I_arr = [0, 45, 90]
    configs = [config for k in range(I_arr)]
    grids = [grid for k in range(I_arr)]
    return None
