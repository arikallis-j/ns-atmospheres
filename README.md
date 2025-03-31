# ns-atmospheres

Simple python library to work with neutron star's atmospheres

## Set-up library
```bash
> pip install git+https://github.com/arikallis-j/ns-atmospheres.git#egg=atmons
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
