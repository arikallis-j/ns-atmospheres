from atmons import *

N_batches, M_batches = 5, 10
N, M = 20, 50

bc = BodyConfig()
body = Body(bc)

slc = SpreadLayerConfig(w_func='const', th_star=40)
slayer = SpreadLayer(slc, body)

gc = GridConfig(n_phi=N, n_theta=N)
grid = Grid(gc, body)

sc = SpectrumConfig(n_nu=M, spec_key='be')
spectrum = Spectrum(sc, body)

ac = AtmosphereConfig(lum=0.1, flux_key='abs')
surface = Surface(grid, body, slayer)

whole_atm = Atmosphere(ac, grid, spectrum, body, slayer, surface)

for key in whole_atm.__dict__.keys():
    whole_atm.__dict__[key] = None

whole_atm.save(name='empty')

atm, atm_g, atm_gs = Atmosphere(), Atmosphere(), Atmosphere()
atm.load(name="empty")
atm_g.load(name="empty")
atm_gs.load(name="empty")

g = 0
for subgrid in grid.get_batch(N_batches):
    s = 0
    subsurf = Surface(subgrid, body, slayer)
    for subspec in spectrum.get_batch(M_batches):
        print(f"Grid...{g}|{N_batches} Spectrum...{s}|{M_batches}")
        cur_atm = Atmosphere(ac, subgrid, subspec, body, slayer, subsurf)
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


cur_atm = Atmosphere(ac, grid, spectrum, body, slayer, surface)
print(cur_atm.B_real - atm.B_real) 


