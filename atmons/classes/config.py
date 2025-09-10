from ..const import *
from dataclasses import dataclass

@dataclass
class BodyConfig:
    name: str = 'J0000+0000'
    chem: str = 's1'
    rel: bool = False#True

    r_ns: float = 12.0 #15.48 
    m_ns: float = 1.5 #1.519
    v_rot: float = 600.0
    i_ang: float = 60.0

@dataclass
class SpreadLayerConfig:
    w_func: str = 'base'
    th_star: float = 45.0
    w_par: tuple = 1
    possible_th: bool = True

@dataclass
class GridConfig:
    n_phi: int =  100
    n_theta: int = 100

    rng_phi: tuple[float, float] = (0.0, 360.0)
    rng_theta: tuple[float, float] = (0.0, 180.0)

    unnull: bool = True

@dataclass
class SpectrumConfig:
    spec_key: str = 'wfc'
    n_nu: int = 500
    rng_erg: tuple[float, float] = (0.1, 20.0)
    zsch_key: bool = True

    fc_key: str = 1

@dataclass
class AtmosphereConfig:
    flux_key: str = 'rel'
    lum: float = 0.1
    n_model: float = N_MODEL
    

@dataclass
class SystemConfig:
    N_batches: int = 100
    M_batches: int = 50


@dataclass
class Config:
    name: str = 'J0000+0000'
    chem: str = 's1'
    rel: bool = False
    r_ns: float = 12.0
    m_ns: float = 1.5
    v_rot: float = 600.0
    i_ang: float = 60.0

    w_func: str = 'base'
    th_star: float = 45.0
    w_par: tuple = 1

    n_phi: int =  100
    n_theta: int = 100
    rng_phi: tuple[float, float] = (0.0, 360.0)
    rng_theta: tuple[float, float] = (0.0, 180.0)
    unnull: bool = True

    spec_key: str = 'wfc'
    n_nu: int = 500
    rng_erg: tuple[float, float] = (0.1, 20.0)
    zsch_key: bool = False
    fc_key: str = 1

    flux_key: str = 'rel'
    lum: float = 0.1
    n_model: float = N_MODEL

    N_batches: int = 100
    M_batches: int = 50

    def __call__(self):
        body = BodyConfig(self.name, 
                          self.chem,
                          self.rel,
                          self.r_ns,
                          self.m_ns,
                          self.v_rot,
                          self.i_ang)
        
        sp_layer = SpreadLayerConfig(self.w_func, 
                                     self.th_star,
                                     self.w_par)
        
        grid = GridConfig(self.n_phi, 
                          self.n_theta,
                          self.rng_phi,
                          self.rng_theta,
                          self.unnull)
        
        spectrum = SpectrumConfig(self.spec_key, 
                                  self.n_nu,
                                  self.rng_erg,
                                  self.zsch_key,
                                  self.fc_key)
        
        atmosphere = AtmosphereConfig(self.flux_key, 
                                      self.lum,
                                      self.n_model)

        system = SystemConfig(self.N_batches,
                              self.M_batches)

        return (body, 
                sp_layer,
                grid,
                spectrum,
                atmosphere,
                system)