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
wget --directory-prefix="spectra" -N $spectra/spec_he.json $spectra/spec_s1.json $spectra/spec_s001.json $spectra/fcol_he.json $spectra/fcol_s1.json $spectra/fcol_s001.json 
```

## Crash Tutorial: make own neutron star
1. Import this library
```python
from atmons import *
```

2. Setup your neutron star parameters 
```python
config = {
    'v_rot': 600,
    'i_ang': 45,
    'chem': 's1',
    'spec_key': 'wfc',
}
```

3. Create your neutron star
```python
ns = build_Neutron_Star(**config)
```

4. Model burst in neutron star's atmosphere!
```python
atm = ns.atmosphere
print(atm.B_real)
```

## Explore neutron star: touch its parameters
```python
ns = build_Neutron_Star()
```

1. Neutron star parameters
```python
body = ns.body
for key, val in body().items():
    print(f"{key} | {val.unit}")
```

2. Neutron star's surface parameters
```python
surf = ns.surface
for key, val in surf().items():
    print(f"{key} | {val.unit}")
```

3. Neutron star's spread layer parameters
```python
sp_layer = ns.sp_layer
for key, val in sp_layer().items():
    print(f"{key} | {val.unit}")
```

4. Neutron star's atmosphere parameters
```python
atm = ns.atmosphere
for key, val in atm().items():
    print(f"{key} | {val.unit}")
```

## Do you own experiments!
```python
class BaseModel(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const()

    def do_experiment(self, experiment='test'):
        ns = build_Neutron_Star()
        shot = ns.atmosphere
        sp_layer = ns.sp_layer

        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w")
        print(f"{'wfc':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f}")

        return shot.lum

model = BaseModel(name = 'discussion')

model()
```