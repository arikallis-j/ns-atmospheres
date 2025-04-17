from atmons import *
## Два стандартных графика

# model
w_func = 'line'
th_star = 89
w_par = 2

# grids
grids = (35, 35, 500)

# neutron star
v_rot = 700
i_ang = 45 # 0 # 90
lum = 0.1
chem = 's1'
flux_key = 1 # 2


# plt.style.use('seaborn-whitegrid')
_, ax = plt.subplots(figsize=(7,7))

# ax.set_title(f"$L_E(\\varepsilon)$, $i = {i_ang}^o$, $l={lum}$, $\\nu = {v_rot}$, $chem = {chem}$", loc='center', fontsize=20)
ax.set_title("$L_E(\\varepsilon)$, $w = line$ | $\\theta_{\\star} = 89^{\\circ}$ ", loc='center', fontsize=20)
ax.set_xlabel("$\\varepsilon, keV$")
ax.set_ylabel("$L_E, 10^{36} erg s^{-1} keV^{-1} sr^{-1} $")
ax.grid(True, which='minor')

if lum==0.1:
    ax.set_ylim(0.003, 0.5)
    ax.set_xlim(1.0, 20)

if lum==0.9:
    ax.set_ylim(0.09, 2)
    ax.set_xlim(1.0, 20)

config, grid = loader()

config['spec_key'] = 'be'
config['chem'] = chem
config['rel'] = True

config['m_ns'] = 1.519
config['r_ns'] = 15.48
config['v_rot'] = v_rot

n_phi, n_theta, n_nu = grids
grid['n_nu'] = n_nu
grid['n_theta'] = n_theta
grid['n_phi'] = n_phi
grid['rng_erg'] = [1.0, 50.0]
config['i_ang'] = i_ang
config['flux_key'] = flux_key
print("N  |  L/L_Edd | T_eff (keV)  | f_c      | w")

dumper('paper', config, grid)
config, grid = loader("paper")
ns = NeutronStar(config,grid)
burst = ns.burst()
counter = 1
Lum = []

for shot in burst:
    if FLUX_REL[counter-1] == lum:
        print(f"{counter:<3}| {shot.lum:.6f} | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue')
        counter += 1
        break
    else:
        counter += 1
        continue


config, grid = loader("paper")
config['spec_key'] = 'wfc'
dumper('paper', config, grid)
config, grid = loader("paper")
ns = NeutronStar(config,grid)
burst = ns.burst()
counter = 1
Lum = []
for shot in burst:
    if FLUX_REL[counter-1] == lum:
        print(f"{counter:<3}| {shot.lum:.6f} | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, linestyle='dashed', color='blue')
        counter += 1
        break
    else:
        counter += 1
        continue


config, grid = loader("paper")
config['spec_key'] = 'be'
config['w_func'] = w_func
config['th_star'] = th_star
config['w_par'] = w_par

dumper('paper', config, grid)
config, grid = loader("paper")
ns = NeutronStar(config,grid)
burst = ns.burst()
W_model = ns.surface.W_model
grv_real = ns.surface.grv_real
omega_kep = ns.param.omega_kep
counter = 1
Lum = []

for shot in burst:
    if FLUX_REL[counter-1] == lum:
        print(f"{counter:<3}| {shot.lum:.6f} | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red')
        counter += 1
        break
    else:
        counter += 1
        continue


config, grid = loader("paper")
config['spec_key'] = 'wfc'

dumper('paper', config, grid)
config, grid = loader("paper")
ns = NeutronStar(config,grid)
burst = ns.burst()
counter = 1
Lum = []
for shot in burst:
    if FLUX_REL[counter-1] == lum:
        print(f"{counter:<3}| {shot.lum:.6f} | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36,linestyle='dashed', color='red')
        counter += 1
        break
    else:
        counter += 1
        continue

print(f"v_kep = {omega_kep / 2 / PI << u.Hz}")
print(f"grv_sign = {np.round(grv_real.value[0,:]/np.abs(grv_real.value[0,:]))}")
print(f"model = {np.round(W_model.value[0,n_theta//2::], decimals=2)}")
plt.savefig('graph/REPORT/line/' + "fig_7b" + f"_line_th{th_star}" +  f"_i{i_ang}" + f"_l{lum}" + f"_{chem}" +  f"_{v_rot}hz" + '.pdf')
plt.show()


## Example
# config, grid = loader()
# config['spec_key'] = 'be'
# config['chem'] = 'he' 

# grid['n_nu'] = 500
# grid['n_theta'] = 10
# grid['n_phi'] = 10

# dumper('paper', config, grid)
# config, grid = loader("paper")

# ns = NeurtonStar(config,grid)
# burst = ns.burst()
# for shot in burst:
#     print(shot.lum)

## Поток по новой интерполяции
# config, grid = loader()
# config['spec_key'] = 'be'
# config['chem'] = 's1'
# config['flux_key'] = 1

# grid['n_nu'] = 50
# grid['n_theta'] = 5
# grid['n_phi'] = 5

# # config['w_func'] = 'const'
# # config['w_par'] = 0.5

# dumper('paper', config, grid)
# config, grid = loader("paper")

# ns = NeutronStar(config,grid)
# burst = ns.burst()
# counter = 1
# Lum = []

# # print(ns.surface.E_dop[0,4,:])
# for shot in burst:
#     print(f"{counter:<3}| {shot.lum:.6f} | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
#     if counter==1:
#         print(f"{np.sum(shot.flux_real , axis=(0,1))}")
#         print(f"{shot.lum}")
#         break
#     counter += 1

