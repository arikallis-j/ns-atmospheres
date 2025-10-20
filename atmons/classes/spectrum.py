from ..funcs import *
from .phenomenon import *

class Spectrum(Phenomenon):
    """Description of Grid parameters"""
    def __init__(self, cfg=None, body=None):
        if cfg is None or body is None:
            return None

        self.spec_key = cfg.spec_key
        self.n_nu = cfg.n_nu 
        self.rng_erg = cfg.rng_erg

        if cfg.zsch_key:
            self.rng_erg = self.rng_erg[0]*body.zsch, self.rng_erg[1]*body.zsch

        # phisical 
        E, dE = E_base(self.n_nu, self.rng_erg)
        self.E = np.array(E) * KEV
        self.dE = np.array(dE) * KEV

        if self.spec_key == 'wfc':
            self.w_b = w_b(body.chem, cfg.fc_key)
            self.T_c = T_c(body.chem, cfg.fc_key)
            self.B_mod = None
        elif self.spec_key == 'be':
            self.w_b = None
            self.T_c = None
            self.B_mod = B_model(body.chem)

        # phisical
        self.shape = self.E.shape
        self.size = self.E.size
    
    def get_batch(self, N_batches):
        part_size = self.size / N_batches
        if part_size < 1:
            part_size = 1
        else:
            part_size = int(part_size)

        for k in range(N_batches):
            a, b = k*part_size, (k+1)*part_size
            if k==N_batches-1:
                b = self.size
            window = (a,b)
            batch = SubSpectrum(self, window)
            yield batch

    
class SubSpectrum(Phenomenon):
    def __init__(self, spec=None, window=(0,0)):
        if spec is None:
            return None
        
        a, b = window
        
        self.spec_key = spec.spec_key
        self.n_nu = spec.n_nu 
        self.rng_erg = spec.rng_erg

        self.w_b = spec.w_b
        self.T_c = spec.T_c
        self.B_mod = spec.B_mod
        self.shape = spec.shape
        self.size = spec.size

        self.E, self.dE  = spec.E[a:b], spec.dE[a:b]