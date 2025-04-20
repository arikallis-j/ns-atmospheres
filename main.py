from atmons import *

## TODO:
# - Как именно мы подгоняем (какой принцип минимизации)?
# - Отрицательная гравитация на экваторе (как отонситься)?
# - Как соотносятся f_c видимые и теоретические (как выправление графика решает проблему низких f_c)?
# - Реальные данные (не пора ли)? 

class BaseLuminosity(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.base = OneStage(mode='base', inter=['wfc', 'be'])
        self.model = OneStage(mode='model', inter=['wfc', 'be'])

    def do_experiment(self, 
            w_func='line', th_star=45, 
            grids = (20, 20, 500),
            i_ang = 45, v_rot=700,
            chem='s1', flux_key = 1,
            mode='base', inter='wfc',
            show=True, save=False, experiment='test',
        ):

        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))
        print(w_func)
        ax.set_title("$L_E(\\varepsilon)$, $w = " + w_func + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        ax.set_xlabel("$\\varepsilon, keV$")
        ax.set_ylabel("$L_E, 10^{36} erg s^{-1} keV^{-1} sr^{-1} $")
        ax.grid(True, which='minor')

        ax.set_ylim(0.003, 2)
        ax.set_xlim(1.0, 20)

        config, grid = loader()
        config['spec_key'] = 'be'
        config['chem'] = chem
        config['rel'] = True

        config['m_ns'] = 1.519
        config['r_ns'] = 15.48
        config['v_rot'] = v_rot

        config['w_func'] = w_func
        config['th_star'] = th_star

        n_phi, n_theta, n_nu = grids
        grid['n_nu'] = n_nu
        grid['n_theta'] = n_theta
        grid['n_phi'] = n_phi
        grid['rng_erg'] = [1.0, 50.0]
        config['i_ang'] = i_ang
        config['flux_key'] = flux_key
        dumper(experiment, config, grid)
        config, grid = loader(experiment)
        ns = NeutronStar(config,grid)
        burst = ns.burst(mode=mode, inter=inter)

        print("N  | L/L_Edd  | L_sl/L_Edd | T_eff (keV)  | f_c      | w")
        counter = 1
        for shot in burst:
            print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
            ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36)
            counter += 1

            if save:
                if not os.path.isdir(f'graph/{experiment}'):
                        os.mkdir(f'graph/{experiment}')
                name = f'graph/{experiment}/'
                name += f"{w_func}"
                name += f"_th{th_star}"
                name += f"_{mode}"
                name += f"_{inter}"

                if i_ang != 45:
                    name += f"_i{i_ang}" 
                if chem != 's1':
                    name += f"_{chem}"
                if v_rot != 700:
                    name += f"_{v_rot}hz"

                name += '.pdf'

                plt.savefig(name)

            if show:
                plt.show()

class BaseSpectra(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const(th_star=45, i_ang=[0,45,89.99])
        self.sqrt = Sqrt(th_star=[45, 89])
        self.line = Line(th_star=[45, 89])
        self.pow_2 = Power(n=2, th_star=[45, 60])
        self.pow_3 = Power(n=3, th_star=[45, 55])
        self.pow_4 = Power(n=4, th_star=[45, 50])

    def do_experiment(self,
            w_func='const', th_star=45, w_par=1.0, 
            grids=(20, 20, 500), 
            v_rot=700, i_ang=45, lum=0.1, 
            chem='s1', flux_key=1, 
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))
        print(w_func)
        ax.set_title("$L_E(\\varepsilon)$, $w = " + w_func + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
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

        config['w_func'] = w_func
        config['th_star'] = th_star
        config['w_par'] = w_par

        n_phi, n_theta, n_nu = grids
        grid['n_nu'] = n_nu
        grid['n_theta'] = n_theta
        grid['n_phi'] = n_phi
        grid['rng_erg'] = [1.0, 50.0]
        config['i_ang'] = i_ang
        config['flux_key'] = flux_key
        print("N  | L/L_Edd  | L_sl/L_Edd | T_eff (keV)  | f_c      | w")

        dumper(experiment, config, grid)
        config, grid = loader(experiment)
        ns = NeutronStar(config,grid)
        burst = ns.burst(mode='base', inter='be')
        counter = 1

        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue')
                counter += 1
                break
            else:
                counter += 1
                continue

        burst = ns.burst(mode='base', inter='wfc')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, linestyle='dashed', color='blue')
                counter += 1
                break
            else:
                counter += 1
                continue

        burst = ns.burst(mode='model', inter='be')
        W_model = ns.surface_model.W_model
        grv_real = ns.surface_model.grv_real
        omega_kep = ns.param.omega_kep
        counter = 1

        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                lambda_lum_be = shot.lambda_lum
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red')
                counter += 1
                break
            else:
                counter += 1
                continue

        burst = ns.burst(mode='model', inter='wfc')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                lambda_lum_wfc = shot.lambda_lum
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36,linestyle='dashed', color='red')
                counter += 1
                break
            else:
                counter += 1
                continue

        print(f"v_kep = {omega_kep / 2 / PI << u.Hz}")
        print(f"grv_sign = {np.round(grv_real.value[0,:]/np.abs(grv_real.value[0,:]))}")
        print(f"model = {np.round(W_model.value[0,n_theta//2::], decimals=2)}")
        print(f"lambda_be = {np.round(lambda_lum_be, decimals=2)}")
        print(f"lambda_wfc = {np.round(lambda_lum_wfc, decimals=2)}")

        if save:
            if not os.path.isdir(f'graph/{experiment}'):
                    os.mkdir(f'graph/{experiment}')
            name = f'graph/{experiment}/'

            name += f"{w_func}"

            name += f"_th{th_star}"

            if i_ang != 45:
                name += f"_i{i_ang}" 
            if lum != 0.1:
                name += f"_l{lum}"
            if chem != 's1':
                name += f"_{chem}"
            if v_rot != 700:
                name += f"_{v_rot}hz"

            name += '.pdf'

            plt.savefig(name)

        if show:
            plt.show()


base_luminosity = BaseLuminosity(
    name = 'test', 
    grids = (2, 2, 5), 
    save = True, 
    show = False,
    w_func = 'line', 
    th_star = 45,
)

base_spectra = BaseSpectra(
    name = 'test', 
    grids = (2, 2, 5), 
    save = True, 
    show = False,
)

# print(base_luminosity)
# base_luminosity()

# print(base_spectra)
# base_spectra()