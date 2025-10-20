from ..funcs import *
from .phenomenon import *

class Grid(Phenomenon):
    """Description of Grid parameters"""
    def __init__(self, cfg=None, body=None):
        if cfg is None or body is None:
            return None

        self.r_0 = body.r_eq
        self.R_0 = body.R_eq
        self.n_phi = cfg.n_phi 
        self.n_theta = cfg.n_theta
        self.rng_phi = cfg.rng_phi 
        self.rng_theta = cfg.rng_theta 

        if cfg.unnull:
            dphi =  (self.rng_phi[1] - self.rng_phi[0]) / (self.n_phi)
            dtheta = (self.rng_theta[1] - self.rng_theta[0]) / (self.n_theta)
            self.rng_phi = (self.rng_phi[0] + dphi/2, self.rng_phi[1] - dphi/2)
            self.rng_theta = (self.rng_theta[0] + dtheta/2, self.rng_theta[1] - dtheta/2)

        # phisical
        self.phi_init = np.full((self.n_theta, self.n_phi), np.linspace(*self.rng_phi, self.n_phi)).T * DEG
        self.theta_init = np.full((self.n_phi, self.n_theta), np.linspace(*self.rng_theta, self.n_theta)) * DEG
        self.phi = self.phi_init << RAD
        self.theta = self.theta_init << RAD
        
        self.sin_th, self.cos_th = sin(self.theta), cos(self.theta)
        self.sin_ph, self.cos_ph = sin(self.phi), cos(self.phi)

        self.r_init = R_metric(self.r_0, self.theta , body.chi, body.Omega)
        self.R = self.r_init << CM
        self.dR = dR_metric(self.R_0, self.sin_th, self.cos_th, body.chi, body.Omega)   

        self.ph_range = np.array([self.rng_phi[0], self.rng_phi[1]]) * DEG << RAD
        self.th_range = np.array([self.rng_theta[0], self.rng_theta[1]]) * DEG << RAD
        
        self.shape = self.phi.shape
        self.size = self.phi.size

    def squeeze(self, parameter):
        arr = self.__dict__[parameter]
        arr_squeezed = arr.reshape(-1)
        return arr_squeezed
    
    def get_batch(self, N_batches):
        flat_grid = FlatGrid(self)
        for batch in flat_grid.get_batch(N_batches):
            yield batch

    
class FlatGrid(Phenomenon):
    def __init__(self, grid=None):
        if grid is None:
            return None
        self.n_phi = grid.n_phi 
        self.n_theta = grid.n_theta
        self.ph_range = grid.ph_range 
        self.th_range = grid.th_range 

        self.shape = grid.shape
        self.size = grid.size

        self.phi, self.theta = grid.squeeze('phi'), grid.squeeze('theta')
        self.R, self.dR = grid.squeeze('R'), grid.squeeze('dR')
        self.sin_ph, self.cos_ph = grid.squeeze('sin_ph'), grid.squeeze('cos_ph')
        self.sin_th, self.cos_th = grid.squeeze('sin_th'), grid.squeeze('cos_th')

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
            batch = SubGrid(self, window)
            yield batch

class SubGrid(Phenomenon):
    def __init__(self, flat_grid=None, window=(0,0)):
        if flat_grid is None:
            return None
        a, b = window

        self.n_phi = flat_grid.n_phi 
        self.n_theta = flat_grid.n_theta
        self.ph_range = flat_grid.ph_range 
        self.th_range = flat_grid.th_range 

        self.shape = flat_grid.shape
        self.size = flat_grid.size

        self.phi, self.theta = flat_grid.phi[a:b], flat_grid.theta[a:b] 
        self.R, self.dR = flat_grid.R[a:b], flat_grid.dR[a:b] 
        self.sin_ph, self.cos_ph = flat_grid.sin_ph[a:b], flat_grid.cos_ph[a:b] 
        self.sin_th, self.cos_th = flat_grid.sin_th[a:b], flat_grid.cos_th[a:b] 
        
        self.phi = self.phi.reshape(*self.phi.shape, 1)
        self.theta = self.theta.reshape(*self.theta.shape, 1)
        self.R = self.R.reshape(*self.R.shape, 1)
        self.dR = self.dR.reshape(*self.dR.shape, 1)
        self.sin_ph = self.sin_ph.reshape(*self.sin_ph.shape, 1)
        self.cos_ph = self.cos_ph.reshape(*self.cos_ph.shape, 1)
        self.sin_th = self.sin_th.reshape(*self.sin_th.shape, 1)
        self.cos_th = self.cos_th.reshape(*self.cos_th.shape, 1)