# ns-atmospheres

Simple python library to work with neutron star's atmospheres.

## Set-up library
1. Install the library
```bash
pip install git+https://github.com/arikallis-j/ns-atmospheres.git#egg=atmons
```

2. Copy in your directory model spectra and configs
```bash
$spectra = "https://raw.githubusercontent.com/arikallis-j/ns-atmospheres/refs/heads/main/spectra"
$configs = "https://raw.githubusercontent.com/arikallis-j/ns-atmospheres/refs/heads/main/configs"
wget --directory-prefix="spectra" -N $spectra/spec_ph.yml $spectra/fcol_He.dat $spectra/fcol_S1.dat $spectra/fcol_S001.dat
wget --directory-prefix="configs" -N $configs/__config__.yml  $configs/__grid__.yml
```

## Crash Tutorial: make own neutron star
1. Import this library
```python
from atmons import *
```

2. Setup your neutron star parameters 
```python
config, grid = loader()
config['spec_key'] = 'wfc'
config['chem'] = 'he' 

grid['n_nu'] = 500
grid['n_theta'] = 10
grid['n_phi'] = 10

dumper('my_neutron_star', config, grid)
config, grid = loader("my_neutron_star")
```

3. Create your neutron star and generator of burst
```python
ns = NeurtonStar(config,grid)
burst = ns.burst()
```

4. Model burst in neutron star's atmosphere!
```python
for shot in burst:
    print(shot.lum)
```

## Explore neutron star: touch its parameters
```python
ns = NeurtonStar(config,grid)
burst = ns.burst()
```

1. Neutron star parameters
```python
param = ns.param()
for key, val in param.items():
    print(f"{key} | {val.unit}")
```

2. Neutron star's surface parameters
```python
surf = ns.surface()
for key, val in surf.items():
    print(f"{key} | {val.unit}")
```

3. Neutron star's shot parameters
```python
shot = ns.shot()
for key, val in shot.items():
    print(f"{key} | {val.unit}")
```