# Pipeline

*Создаем нейтронную звезду с нуля*

## Нейтронная звезда как тело

Для начала, создадим саму нейтронную звезду.  
Нам необходимо задать следующие параметры:

`name` : название нейтронной звезды (по дефолту: 'J0000+0000')  
`chem` : химический состав  

    - s1: солнечная металличность  
    - s001: металичность 1/100 солнечной  
    - he: гелиевая атмосфера  

`rel` : задаем ли мы массу и радиус сразу с поправками ОТО, или же нет  
`r_ns` : радиус нейтронной звезды в километрах (или же экваториальный радиус)  
`m_ns` : масса нейтронной звезды в массах солнца (или же скорректированная масса)  
`v_rot` : частота вращения нейтронной звезды в герцах  
`i_ang` : угол наклона оси вращения к лучу зрения в градусах  

```python
config_body = {
    'name': 'J0000+0000',
    'chem': 's1',
    'rel': False,
    'r_ns': 12.0,
    'm_ns': 1.5,
    'v_rot': 600.0,
    'i_ang': 60.0,
}

bc = BodyConfig(**config_body)
body = Body(bc)
```

В этот момент создается следующий объект:

```python
class Body(Phenomenon):
    def __init__(self, cfg=None):
        if cfg is None:
            return None
        
        self.name = cfg.name

        # basic
        ...

        # physical
        ...

        # chemical
        ...

        # photometrical
        ...

        # rotational
        ...

        # relativical
        ...

        # metrical
        ...

        # keplerian
        ...
```

Разберем по порядку. 

### basic

```python
r_ns = cfg.r_ns * KM
m_ns = cfg.m_ns * M_SUN
v_rot =  cfg.v_rot * HZ
chem = cfg.chem

if cfg.rel:
    r = rel_r_eq(r_ns, m_ns, v_rot)
    m = rel_m_cor(m_ns, r_ns, v_rot)
    r_eq, m_cor = r_ns, m_ns
else:
    r, m = r_ns, m_ns
    r_eq = r_eq(r, m, v_rot)
    m_cor = m_cor(m, r, v_rot)

R = R_NS(r, m)
M = M_NS(m)
```
Помимо инциализации, здесь также рассчитываются экваториальный радиус и скорректированная масса нейтронной звезды по следующим формулам:
$$
R_e = R \left[ 0.9766 + \frac{0.025}{1.07 - \bar{\nu}}  + 0.07 M^{3/2}_{1.4} \bar{\nu}^2 \right] 
$$

$$
M' = M \left[ a_0 + \frac{a_1}{1.1 - \bar{\nu}}  + a_2 \bar{\nu}^2 \right]
$$

Подробнее смотрите в Suleimanov et al. (2020), выражения (2) и (3).


### physical

```python
R_sch = R_sch(M)
zsch = zsch(R, R_sch)
area_0 = Surf(R, zsch)
g = g(R, M, zsch)
log_g = log(g / g.unit)
```

$$
R_{s} = 2 \frac{G M}{c^2} 
$$

$$
z + 1 = \frac{1}{\sqrt{1 - R_s/R}}
$$

$$
S = (R \cdot (1 + z))^2
$$

$$
g = \frac{GM}{R^2} \cdot (1 + z)
$$


### chemical

```python
X_hyd = X_hyd(chem)
kappa_e = kappa_e(X_hyd)
```
$$
X = 
\begin{cases}
0.7374 & \text{если } \text{'s1'} \text{ или } \text{'s001'} \\
0 & \text{если } \text{'he'} 
\end{cases}
$$

$$
\kappa_{e} =  0.2 \cdot (1.0 + X) \cdot \text{cm}^2 \text{g}^{-1}
$$

### photometrical
```python
Flux_edd = Flux_edd(g, kappa_e)
T_edd = T_obs(Flux_edd, zsch)
Theta_edd = Theta(T_edd)
Epsilon_edd = Epsilon(T_edd)
Lum_edd = Lumen(Flux_edd, R)
Lum_obs = Lumen_obs(Flux_edd, R, zsch)
```
$$
F_{Edd} =  \frac{gc}{\kappa}
$$

$$
T_{Edd} =  \left(\frac{F_{Edd}}{\sigma_{SB}}\right)^{1/4} \cdot \frac{1}{1 + z}
$$

$$
\Theta_{Edd} =  k_B \cdot T_{Edd} 
$$

$$
L_{Edd} =   4 \pi R^2 \cdot F_{Edd} 
$$

$$
L_{obs} =   4 \pi R^2 \cdot \frac{F_{Edd}}{(1 + z)^2}
$$

### rotational

```python
nu_rot = v_rot
incl_ang = (cfg.i_ang * DEG).to(RAD)
sin_i = sin(incl_ang)
cos_i = cos(incl_ang)
omega_rot = omega(nu_rot)
```
$$
\omega_{rot} = 2 \pi \nu_{rot}
$$

### relativical

```python
v_cr = nu_crit(r, m)
v_rel = nu_relative(nu_rot, v_cr)
M_cor = M_NS(m_cor)
R_eq = R_NS(r_eq, m_cor)
R_sch_cor = R_sch(M_cor)
```

$$
\nu_{cr} = 1278 \cdot \left(\frac{10\text{ км}}{r_{ns}}\right)^{1.5} \cdot \sqrt{\frac{m_{ns}}{1.4M_{\odot}}}
$$

$$
\nu_{rel} = \frac{\nu_{rot}}{\nu_{crit}}
$$


### metrical

```python
chi = chi_metric(R_eq, M_cor)
Omega = Omega_metric(R_eq, M_cor, omega_rot)
q_c = q_c_metric(chi, Omega)
b_c = b_c_metric(chi, Omega)
i_bar = i_bar_metric(chi)
I = I_NS(i_bar, R_eq, M_cor)
J = J_NS(I, omega_rot)
g_0 = g_0_metric(R_eq, M_cor, chi)
```

$$
\chi = \frac{G M_{cor}}{R_{eq} c^2}
$$

$$
\Omega = \omega_{rot} \cdot \sqrt{\frac{R_{eq}^3}{G M_{cor}}}
$$

$$
q_c = -0.11 \cdot \left(\frac{\Omega}{\chi}\right)^2 
$$

$$
b_c = 0.4454 \cdot \Omega^2 \cdot \chi
$$


$$
i_m = \sqrt{\chi} \cdot (1.136 - 2.53 \chi + 5.6 \chi^2)
$$

$$
I = i_m M_{cor} R_{eq}^2 
$$

$$
J = I \omega_{rot} 
$$

$$
g_0 = \frac{G M_{cor}}{R_{eq}^2} \cdot \frac{1}{\sqrt{1 - 2 \chi}}
$$

### keplerian

```python
V_rot = V_rot(omega_rot, R_eq)
V_kep = V_kep(g_0, R_eq)
omega_kep = omega_rot(V_kep, R_eq)
omega_cr = omega(v_cr)
```
$$
V_{rot} = \omega_{rot} R_{eq}
$$

$$
V_{kep} = \sqrt{g_0 \cdot R_{eq}}
$$

$$
\omega_{kep} = \frac{V_{kep}}{R_{eq}}
$$

$$
\omega_{cr} = 2 \pi \nu_{cr}
$$

## Сетка координат

Далее, нам нужно создать сетку координат на нейтронной звезде, которая будет учитывать неоднородность радиуса на поверхности.  
Нам необходимо задать следующие параметры:

`n_phi` : количество точек по долготе   
`n_theta` : количество точек по широте  
`rng_phi` : диапазон значений для долготы  
`rng_theta` : диапазон значений для широты  
`unnull` : корректировка сетки для ненулевых значений на границе (по дефолту True)  

```python
config_grid = {
    'n_phi': 30,
    'n_theta': 30,
    'rng_phi': (0.0, 360.0),
    'rng_theta': (0.0, 180.0),
    'unnull': True,
}

gc = GridConfig(**config_grid)
grid = Grid(gc, body)
```

В этот момент создается следующий объект:

```python
class Grid(Phenomenon):
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
```

Из примечательного здесь: вычисление радиуса $R$ и производной радиуса $dR$ по соответсвующим формулам. 

Главное, что здесь задаются основные матрицы `phi`, `theta`,`R` с которыми мы в дальнейшем будем работать.

P.S. За счёт такого задания сетки мы можем реализовать вычисление спектров на ограниченной части поверхности - например на кольце по широте.

## Слой растекания

Далее, нам нужно создать параметры слоя растекания, который мы будем расчитывать на поверхности нейтронной звезды.   
Нам необходимо задать следующие параметры:

`w_func` : единичный профиль скорости (название функции)  
`th_star` : ширина слоя растекания по широте  
`w_par` : параметр слоя растекания (максимальная скорость на экваторе)  


```python
config_splayer = {
    'w_func': 'base',
    'th_star': 45.0,
    'w_par': 1,
}

slc = SpreadLayerConfig(**config_splayer)
sp_layer = SpreadLayer(slc, body)
```

В этот момент создается следующий объект:

```python
class SpreadLayer(Phenomenon):
    """Description of SpreadLayer parameters"""
    def __init__(self, cfg=None,  body=None):
        if cfg is None or body is None:
            return None
        
        self.w_func = cfg.w_func
        self.w_par = cfg.w_par
        self.th_star = cfg.th_star * DEG << RAD

        # different factors
        kep_part = self.w_par
        chi_omega = body.omega_rot / body.omega_kep
        kep_part = max(kep_part, chi_omega)

        self.kep_part = kep_part
        self.omega_kep_local = body.omega_kep * kep_part
        self.rel_omega = self.omega_kep_local/body.omega_rot
```
Здесь помимо инициализации проверяется, что максимальная угловая скорость на экваторе не меньше чем скорость угловая скорость вращения нейтронной звезды.

## Поверхность нейтронной звезды

Далее, нам нужно расчитать параметры поверхности нейтронной звезды, на созданной нами сетке координат с учетом слоя растекания.  
Нам необходимо задать следующие параметры:

`kep_part_fn` : ограничение на максимальную скорость (название функции)  
`log_g_cr` : минимально возможный логарифм гравитации 

```python
config_surface = {
    'kep_part_fn': 'none',
    'log_g_cr': 13.7,
}

sf = SurfaceConfig(**config_surface)
surface = Surface(sf, grid, body, sp_layer)
```


В этот момент создается следующий объект:

```python
class Surface(Phenomenon):
    def __init__(self, cfg=None, grid = None, body = None, sp_layer = None):
        if cfg is None or grid is None or body is None or sp_layer is None:
            return None
        
        phi, theta, R, dR = grid.phi, grid.theta, grid.R, grid.dR
        sin_ph, cos_ph = grid.sin_ph, grid.cos_ph
        sin_th, cos_th = grid.sin_th, grid.cos_th
        th_star = sp_layer.th_star

        # keplerian part
        ...

        # spread layer
        ...
            
        # metrical
        ...

        # gravity
        ...

        # radiational
        ...

        # rotational
        ...

        # integration
        ...
```

Разберем по порядку.

### keplerian part

```python
Omega_bs =  Omega_metric(body.R_eq, body.M_cor, body.omega_rot) 
g_th_0 = g_metric(1, 0, body.chi, Omega_bs, Omega_bs, body.i_bar)
kappa_teory = np.sqrt(G_GRAV*body.M_cor/(body.R_eq * body.V_kep**2))
chi_i = (1 + body.chi * (-1 + 2*body.i_bar) + body.chi**2 * (-2 + 4*body.i_bar - 8*body.i_bar**2))
kappa_teory = kappa_teory * np.sqrt(g_th_0/chi_i + Omega_bs**2)
self.kappa_max = kappa_teory

Omega_bs =  Omega_metric(body.R_eq, body.M_cor, body.omega_rot) 
g_th_null = 10.0**13.7 / body.g_0.value
g_th_0 = g_metric(1, 0, body.chi, Omega_bs, Omega_bs, body.i_bar)
kappa_teory = np.sqrt(G_GRAV*body.M_cor/(body.R_eq * body.V_kep**2))
chi_i = (1 + body.chi * (-1 + 2*body.i_bar) + body.chi**2 * (-2 + 4*body.i_bar - 8*body.i_bar**2))
kappa_teory = kappa_teory * np.sqrt((g_th_0 - g_th_null)/chi_i + Omega_bs**2)
self.kappa_int = kappa_teory

self.kappa_cr = body.omega_cr / body.omega_kep

kep_part_local = sp_layer.kep_part
if cfg.kep_part_fn=='max':
    kep_part_local = min(kep_part_local, self.kappa_max)
if cfg.kep_part_fn=='int':
    kep_part_local = min(kep_part_local, self.kappa_int)
if cfg.kep_part_fn=='cr':
    kep_part_local = min(kep_part_local, self.kappa_cr)
if cfg.kep_part_fn=='none':
    kep_part_local = sp_layer.kep_part

omega_kep_local = body.omega_kep * kep_part_local
self.kep_part = kep_part_local
```

### spread layer
```python
self.W_model = W_model(R, body.R_eq, theta, th_star, sp_layer.w_func, sp_layer.w_par, omega_kep_local, body.omega_rot)
self.W_base = np.ones(self.W_model.shape) * self.W_model.unit
self.omega_model = body.omega_rot * self.W_model
self.omega_base = body.omega_rot * self.W_base

self.psi = abs(90 * DEG - theta) << RAD
self.spread_layer = self.psi <= th_star
self.spread_layer_base = self.psi < 0.0
self.spread_layer_true = self.spread_layer
self.Omega_model = Omega_metric(body.R_eq, body.M_cor, self.omega_model)
self.Omega_base = Omega_metric(body.R_eq, body.M_cor, self.omega_base)

if sp_layer.w_func=='base':
    self.W_model = self.W_base
    self.omega_model = self.omega_base
    self.Omega_model = self.Omega_base
    self.spread_layer = self.spread_layer_base
```

### metrical
```python
self.u = u_metric(R, body.R_sch_cor)
self.r_bar, self.u_bar = r_u_metric(R, cos_th, body.q_c, body.b_c, body.R_sch_cor)
self.nu, self.B, self.zeta = nu_B_dzeta_metric(cos_th, self.u_bar, body.q_c, body.b_c)
self.omega_bar = omega_bar_metric(self.r_bar, self.u_bar, body.J)
self.beta_ph = beta_ph_metric(R, sin_th, self.omega_bar, self.nu)
self.g_th = g_metric(sin_th, cos_th, body.chi, self.Omega_base, self.Omega_model, body.i_bar)
```

### gravity

```python
g = self.g_th * body.g_0.value
g = g * hs(g - 1) + 1 * hs(1 - g)
log_g = log(g)
log_g_crit = cfg.log_g_cr
g_th_null = 10.0**log_g_crit / body.g_0.value
if (log_g <= log_g_crit).any():
    print("Incorrect gravity")
self.g_th = np.where(log_g > log_g_crit, self.g_th, self.g_th * 0.0 + g_th_null)
self.g_th_base = g_metric(sin_th, cos_th, body.chi, self.Omega_base, self.Omega_base, body.i_bar)

self.f_th = f_theta(R, dR, self.nu, self.B, self.zeta)
self.sin_eta, self.cos_eta = eta_metric(self.f_th)
self.beta = beta_metric(R, sin_th, self.nu, self.omega_bar, self.omega_model)
self.gamma = gamma_metric(self.beta)

self.grv = grv_metric(theta, self.g_th, body.g_0)
self.grv_base = grv_metric(theta, self.g_th_base, body.g_0)
self.log_g = log(self.grv / self.grv.unit)
self.log_g_base = log(self.grv_base / self.grv_base.unit)
```

### radiational
```python
self.Flux_edd = Flux_edd(self.grv, body.kappa_e)
self.Flux_edd_base = Flux_edd(self.grv_base, body.kappa_e)
```

### rotational
```python
self.sin_psi, self.cos_psi = psi_rot(sin_th, cos_th, cos_ph, body.sin_i, body.cos_i)
self.G_yu = G_yu_rot(self.cos_psi, self.u)
self.D = D_rot(self.cos_psi, self.u)
self.sin_a, self.cos_a = alpha_rot(self.cos_psi, self.u, self.G_yu)
self.cos_chi = chi_rot(sin_th, cos_th, self.sin_psi, self.cos_psi, body.cos_i)
self.cos_sig = sigma_rot(self.sin_eta, self.cos_eta, self.sin_a, self.cos_a, self.cos_chi, cos_th)
self.cos_xi = xi_rot(self.sin_a, self.sin_psi, sin_ph, body.sin_i)
self.delta = delta_rot(self.beta, self.gamma, self.cos_xi)
self.cos_sig_1 = sigma_1_rot(self.cos_sig, self.delta)
```

### integration
```python
self.dS = dS_metric_1(theta, self.cos_eta, R, grid.n_phi, grid.n_theta, grid.ph_range, grid.th_range)

self.dOmega = self.dS
self.dOmega_obs = dOmega_rot(self.dS, self.cos_sig, self.D)
self.dOmega_obs_real = np.where(np.logical_not(self.cos_sig < 0.0), self.dOmega_obs, np.zeros(self.dOmega_obs.shape))

self.area = np.sum(self.dOmega)
self.area_real = np.sum(self.dOmega_obs_real)

ph_min, ph_max = grid.ph_range
self.l_phi = (ph_max - ph_min)

self.R_pr = sqrt(self.area / (self.l_phi / RAD))
```