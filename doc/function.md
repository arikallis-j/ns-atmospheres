# function

## basic
> Description of basic relationships
 
### Mass and Radius in CGS system
```python
def r_min(m_ns: Q[M_SUN]) -> Q[KM]:
    """Minimum neutron star radius (from ???)"""
```

```python
def R_NS(r_ns: Q[KM], m_ns: Q[M_SUN]) -> Q[CM]:
    """Neutron star radius (by definition)"""
```

```python
def M_NS(m_ns: Q[M_SUN]) -> Q[GRAM]:
    """Neutron star mass (by definition)"""
```

### Critical frequency of rotation
```python
def omega(nu: Q[HZ]) -> Q[SEC**(-1)]:
    """Cyclical frequency (by definition)"""
```

```python
def nu_crit(r_ns: Q[KM], m_ns: Q[M_SUN]) -> Q[HZ]:
    """Maximum possible rotation frequency (II.eq.1)"""
```

```python
def nu_relative(nu_rot: Q[HZ], nu_crit: Q[HZ]) -> Q[0]:
    """Relative rotation frequency (by definition)"""
```

### Velocity of rotation 
```python
def omega_rot(V: Q[CM / SEC], R: Q[CM]) -> Q[SEC**(-1)]:
    """Angular frequency of rotation (by definition)"""
```

```python
def V_rot(omega: Q[SEC**(-1)], R: Q[CM]) -> Q[CM / SEC]:
    """Linear velocity of rotation (by definition)"""
```

```python
def V_kep(g: Q[CM/SEC**2], R: Q[CM]) -> Q[CM / SEC]:
    """Linear velocity of keplerian movement"""
```

### Mass and Radius in relative units
```python
def r_eq(r_ns: Q[KM], m_ns: Q[M_SUN], v_rot: Q[HZ]) -> Q[KM]:
    """Equatorial radius of Neutron Star (II.eq.2)"""
```

```python
def m_cor(m_ns: Q[M_SUN], r_ns: Q[KM], v_rot: Q[HZ]) -> Q[M_SUN]:
    """Corrected mass of Neutron Star (II.eq.3)"""
```

### Inverse Mass and Radius in relative units
```python
def rel_r_eq(r_ns: Q[KM], m_ns: Q[M_SUN], v_rot: Q[HZ]) -> Q[KM]:   
    """Inverse of corrected mass of Neutron Star (II.eq.3)"""
```

```python
def rel_m_cor(m_ns: Q[M_SUN], r_ns: Q[KM], v_rot: Q[HZ]) -> Q[M_SUN]:
    """Inverse of equatorial radius of Neutron Star (II.eq.2)"""
```

## math


