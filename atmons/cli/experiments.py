from ..ns.classes import NeurtonStar

def lfc(cfg, grid):
    NS = NeurtonStar(cfg, grid)
    lum, fc = [], []
    for shot in NS.burst():
        lum.append(shot.lum)
        fc.append(shot.fc)
    return lum, fc

def lw(cfg, grid):
    NS = NeurtonStar(cfg, grid)
    lum, w = [], []
    for shot in NS.burst():
        lum.append(shot.lum)
        w.append(shot.w)
    return lum, w

def spectra(cfg, grid):
    NS = NeurtonStar(cfg, grid)
    spec, erg = [], []
    for shot in NS.burst():
        spec.append(shot.B_real)
        erg.append(shot.E_null)
    return None
