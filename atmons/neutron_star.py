from .function import *
from .classes import *

@dataclass
class NeutronStar:
    body: Body
    sp_layer: SpreadLayer
    grid: Grid
    spectrum: Spectrum
    surface: Surface
    atmosphere: Atmosphere

@dataclass
class NeutronStarSurface:
    body: Body
    sp_layer: SpreadLayer
    grid: Grid
    surface: Surface

def build_Neutron_Star(
        v_rot=600, i_ang=45, chem='s1',
        N_ph=30, N_th=30, N_nu=500,
        w_func='const', th_star=40, w_par=1.0,
        spec_key = 'wfc', flux_key='rel',
        lum=0.1, r_ns=12, m_ns=1.4,
    ):


    bc = BodyConfig(chem=chem, v_rot=v_rot, i_ang=i_ang, r_ns=r_ns, m_ns=m_ns)
    body = Body(bc)

    slc = SpreadLayerConfig(w_func=w_func, th_star=th_star, w_par=w_par)
    sp_layer = SpreadLayer(slc, body)

    gc = GridConfig(n_phi=N_ph, n_theta=N_th)
    grid = Grid(gc, body)
    surface = Surface(grid, body, sp_layer)

    sc = SpectrumConfig(n_nu=N_nu, spec_key=spec_key)
    spectrum = Spectrum(sc, body)

    ac = AtmosphereConfig(lum=lum, flux_key=flux_key)
    atmosphere = Atmosphere(ac, grid, spectrum, body, sp_layer, surface)
    neutron_star = NeutronStar(body, sp_layer, grid, spectrum, surface, atmosphere)

    return neutron_star

def build_Surface(
        v_rot=600, i_ang=45, chem='s1',
        N_ph=30, N_th=30, unnull = True,
        w_func='const', th_star=40, w_par=1.0,
    ):


    bc = BodyConfig(chem=chem, v_rot=v_rot, i_ang=i_ang)
    body = Body(bc)

    slc = SpreadLayerConfig(w_func=w_func, th_star=th_star, w_par=w_par)
    sp_layer = SpreadLayer(slc, body)

    gc = GridConfig(n_phi=N_ph, n_theta=N_th, unnull=unnull)
    grid = Grid(gc, body)
    surface = Surface(grid, body, sp_layer)

    neutron_star_surface = NeutronStarSurface(body, sp_layer, grid, surface)

    return neutron_star_surface


class NeutronStarComplex(Phenomenon):
    """Description of Atmosphere parameters"""
    def __init__(self, body_cfg = None, sp_layer_cfg = None, 
                      grid_cfg = None, spec_cfg = None, 
                      atmos_cfg = None, sys_cfg = None):
        
        if sys_cfg is None:
            sys_cfg = SystemConfig()
        N_batches, M_batches = sys_cfg.N_batches, sys_cfg.M_batches
        
        if body_cfg is None:
            body_cfg = BodyConfig()
        body = Body(body_cfg)

        if sp_layer_cfg is None:
            sp_layer_cfg = SpreadLayerConfig()
        sp_layer = SpreadLayer(sp_layer_cfg, body)

        if grid_cfg is None:
            grid_cfg = GridConfig()
        grid = Grid(grid_cfg, body)

        if spec_cfg is None:
            spec_cfg = SpectrumConfig()
        spectrum = Spectrum(spec_cfg, body)

        if atmos_cfg is None:
            atmos_cfg = AtmosphereConfig()

        surface = Surface(grid, body, sp_layer)

        # whole_atm = Atmosphere(ac, grid, spectrum, body, slayer, surface)

        # for key in whole_atm.__dict__.keys():
        #     whole_atm.__dict__[key] = None

        # whole_atm.save(name='empty')

        atm, atm_g, atm_gs = Atmosphere(), Atmosphere(), Atmosphere()
        atm.load(name="empty")
        atm_g.load(name="empty")
        atm_gs.load(name="empty")

        g = 0
        for subgrid in grid.get_batch(N_batches):
            s = 0
            subsurf = Surface(subgrid, body, sp_layer)
            for subspec in spectrum.get_batch(M_batches):
                print(f"Grid...{g}|{N_batches} Spectrum...{s}|{M_batches}")
                cur_atm = Atmosphere(atmos_cfg, subgrid, subspec, body, sp_layer, subsurf)
                if s==0:
                    for key in atm.__dict__.keys():
                        param = cur_atm.__dict__[key]
                        atm_gs.__dict__[key] = param
                else:
                    for key in atm.__dict__.keys():
                        param = cur_atm.__dict__[key]
                        if isinstance(param, np.ndarray):
                            if len(param.shape) == 3:
                                atm_gs.__dict__[key] = np.concatenate((atm_gs.__dict__[key], param), axis=2)
                        else:
                            atm_gs.__dict__[key] = param
                s += 1

            if g==0:
                for key in atm.__dict__.keys():
                    param = atm_gs.__dict__[key]
                    atm_g.__dict__[key] = param
            else:
                for key in atm.__dict__.keys():
                    param = atm_gs.__dict__[key]
                    if isinstance(param, np.ndarray):
                        if len(param.shape) == 3:
                            atm_g.__dict__[key] = np.concatenate((atm_g.__dict__[key], param), axis=0)
                        if len(param.shape) == 2:
                            atm_g.__dict__[key] = np.concatenate((atm_g.__dict__[key], param), axis=0)
                    else:
                        atm_g.__dict__[key] = param
            g += 1

        for key in atm.__dict__.keys():
            param = atm_g.__dict__[key]
            if isinstance(param, np.ndarray):
                if len(param.shape) == 3:
                    atm.__dict__[key] = param.reshape((*grid.shape,*spectrum.shape))
                if len(param.shape) == 2:
                    atm.__dict__[key] = param.reshape((*grid.shape,))
                if len(param.shape) == 0:
                    atm.__dict__[key] = param
            else:
                atm.__dict__[key] = param

        atm.calc()

        self.system = sys_cfg
        self.body = body
        self.sp_layer = sp_layer
        self.grid = grid
        self.spectrum = spectrum
        self.surface = surface
        self.atmosphere = atm