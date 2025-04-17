from ..ns.classes import NeutronStar

def lfc(cfg, grid):
    NS = NeutronStar(cfg, grid)
    lum, fc = [], []
    for shot in NS.burst():
        lum.append(shot.lum)
        fc.append(shot.fc)
    return lum, fc

def lw(cfg, grid):
    NS = NeutronStar(cfg, grid)
    lum, w = [], []
    for shot in NS.burst():
        lum.append(shot.lum)
        w.append(shot.w)
    return lum, w

def spectra(cfg, grid):
    NS = NeutronStar(cfg, grid)
    spec, erg = [], []
    for shot in NS.burst():
        spec.append(shot.B_real)
        erg.append(shot.E_null)
    return None
