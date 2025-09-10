from atmons import *

<<<<<<< HEAD
class BaseModel(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const()

    def do_experiment(self, experiment='test'):
        config = {
            'v_rot': 600,
            'i_ang': 45,
            'chem': 's1',
            'spec_key': 'wfc',
        }
        ns = build_Neutron_Star(**config)
        shot = ns.atmosphere
        sp_layer = ns.sp_layer

        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        print(f"{'wfc':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")

        return 0


model = BaseModel(name = 'discussion')

model()
=======
## TODO:
# - Как именно мы подгоняем (какой принцип минимизации)?
# - Отрицательная гравитация на экваторе (как отонситься)?
# - Как соотносятся f_c видимые и теоретические (как выправление графика решает проблему низких f_c)?
# - Реальные данные (не пора ли)? 



class TestModel(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.theta_1 = OneStage(th_star=1)
        self.theta_15 = OneStage(th_star=15)
        self.theta_30 = OneStage(th_star=30)
        self.theta_45 = OneStage(th_star=45)
        self.theta_60 = OneStage(th_star=60)
        self.theta_75 = OneStage(th_star=75)
        self.theta_89 = OneStage(th_star=89)

    def do_experiment(self,
            th_star=45, 
            grids=(20, 20, 500), 
            v_rot=700, i_ang=45, lum=0.1, 
            chem='s1', flux_key=1, 
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))

        ax.set_title("$L_E(\\varepsilon)$, $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
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

        print("N  | L/L_Edd  | L_sl/L_Edd | T_eff (keV)  | f_c      | w")

        config['w_func'] = 'const'
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='purple', label=config['w_func'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['w_func'] = 'power-4'
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue', label=config['w_func'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['w_func'] = 'power-3'
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='cyan', label=config['w_func'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['w_func'] = 'power-2'
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='green', label=config['w_func'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['w_func'] = 'line'
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='yellow', label=config['w_func'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['w_func'] = 'sqrt'
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='orange', label=config['w_func'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['w_func'] = 'none'
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red', label='be')
                counter += 1
                break
            else:
                counter += 1
                continue

        burst = ns.burst(inter='wfc')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='black', linestyle='dashed', label='wfc')
                counter += 1
                break
            else:
                counter += 1
                continue

        ax.legend()

        if save:
            if not os.path.isdir(f'graph/{experiment}'):
                    os.mkdir(f'graph/{experiment}')
            name = f'graph/{experiment}/'

            name += f"th{th_star}"

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

test_model = TestModel(
    name = 'test-model', 
    grids = (20, 20, 500), 
    save = True, 
    show = False,
)

# print(test_model)
# test_model()

class TestTheta(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        # self.none = OneStage(w_func='none')
        self.sqrt = OneStage(w_func='sqrt')
        self.line = OneStage(w_func='line')
        self.pow_2 = OneStage(w_func='power-2')
        self.pow_3 = OneStage(w_func='power-3')
        self.pow_4 = OneStage(w_func='power-4')
        self.pow_9 = OneStage(w_func='power-9')
        self.exp = OneStage(w_func='exp')
        self.const = OneStage(w_func='const')

    def do_experiment(self,
            w_func='none', 
            grids=(20, 20, 500), 
            v_rot=700, i_ang=45, lum=0.1, 
            chem='s1', flux_key=1, 
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))

        ax.set_title("$L_E(\\varepsilon)$, $w = " + w_func + "$", loc='center', fontsize=20)
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

        n_phi, n_theta, n_nu = grids
        grid['n_nu'] = n_nu
        grid['n_theta'] = n_theta
        grid['n_phi'] = n_phi
        grid['rng_erg'] = [1.0, 50.0]
        config['i_ang'] = i_ang
        config['flux_key'] = flux_key
        dumper(experiment, config, grid)
        config, grid = loader(experiment)

        print("N  | L/L_Edd  | L_sl/L_Edd | T_eff (keV)  | f_c      | w")
        config['th_star'] = 89
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{config['th_star']:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='purple', label=config['th_star'])
                counter += 1
                break
            else:
                counter += 1
                continue
        
        config['th_star'] = 75
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{config['th_star']:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue', label=config['th_star'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['th_star'] = 60
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{config['th_star']:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='cyan', label=config['th_star'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['th_star'] = 45
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{config['th_star']:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='green', label=config['th_star'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['th_star'] = 30
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{config['th_star']:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='yellow', label=config['th_star'])
                counter += 1
                break
            else:
                counter += 1
                continue

        config['th_star'] = 15
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{config['th_star']:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='orange', label=config['th_star'])
                counter += 1
                break
            else:
                counter += 1
                continue
        

        config['th_star'] = 1
        ns = NeutronStar(config, grid)
        burst = ns.burst(inter='be', mode='base')
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{"be":<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red', label='be')
                counter += 1
                break
            else:
                counter += 1
                continue

        burst = ns.burst(inter='wfc', mode='base')
        num = 0
        counter = 1
        for shot in burst:
            if FLUX_REL[counter-1] == lum:
                print(f"{"wf":<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
                ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='black', linestyle='dashed', label='wfc')
                counter += 1
                break
            else:
                counter += 1
                continue

        ax.legend()


        if save:
            if not os.path.isdir(f'graph/{experiment}'):
                    os.mkdir(f'graph/{experiment}')
            name = f'graph/{experiment}/'

            name += f"{w_func}"

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

test_theta = TestTheta(
    name = 'test-theta', 
    grids = (32, 32, 500), 
    save = True, 
    show = False,
)

# print(test_theta)
# test_theta()

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
        self.const = Const(w_func="const", th_star = [45,50,55])
        # self.sqrt = Sqrt(th_star=[45, 89])
        # self.line = Line(th_star=[45, 89])
        # self.pow_2 = Power(n=2, th_star=[45, 60])
        # self.pow_3 = Power(n=3, th_star=[45, 55])
        # self.pow_4 = Power(n=4, th_star=[45, 50])

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

        # else:
        #     ax.set_ylim(0.003, 2)
        #     ax.set_xlim(1.0, 20)

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
        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w")

        dumper(experiment, config, grid)
        config, grid = loader(experiment)
        ns = NeutronStar(config,grid)
        Flux_full_base = None
        Flux_full_sl_base = None
        burst = ns.burst(mode='base', inter='wfc', lum=lum)
        counter = 1
        for shot in burst:
            Flux_full_base = shot.Flux_full
            Flux_full_sl_base = shot.Flux_full_sl
            print(f"{"bwf":<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} ")
            ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, linestyle='dashed', color='blue')


        burst = ns.burst(mode='base', inter='be', lum=lum)
        counter = 1
        for shot in burst:
            print(f"{"bbe":<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} ")
            ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue')

        burst = ns.burst(mode='model', inter='be', lum=lum)
        W_model = ns.surface_model.W_model
        grv_real = ns.surface_model.grv_real
        omega_kep = ns.param.omega_kep
        B_real = None
        counter = 1

        for shot in burst:
            lambda_lum_be = shot.lambda_lum
            print(f"{"mbe":<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} ")
            ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red')
            B_real = shot.B_real


        # Flux_full = None
        # Flux_full_sl = None
        # burst = ns.burst(mode='model', inter='wfc', lum=lum)
        # counter = 1
        # for shot in burst:
        #     Flux_full = shot.Flux_full
        #     Flux_full_sl = shot.Flux_full_sl
        #     print(f"{"mwf":<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
        #     ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, linestyle='dashed', color='red')

            # print(ns.surface.spread_layer)
            # print(shot.flux[0,:])
            # print(shot.B_real)
            # print(ns.surface.log_g)
            # print(ns.surface.E_real)
            # print(shot.rho)
        # print(f"v_kep = {omega_kep / 2 / PI << u.Hz}")
        # print(f"grv_sign = {np.round(grv_real.value[0,:]/np.abs(grv_real.value[0,:]))}")
        # print(f"model = {np.round(W_model.value[0,n_theta//2::], decimals=2)}")
        # print(f"lambda_be = {np.round(lambda_lum_be, decimals=2)}")
        # print(f"lambda_wfc = {np.round(lambda_lum_wfc, decimals=2)}")

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
    name = 'base-spectra', 
    grids = (32, 32, 500), 
    save = True, 
    show = True,
    lum = 0.1,
    w_func='const',
)

# print(base_luminosity)
# base_luminosity()

print(base_spectra)
base_spectra()


class BaseEnergy(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const(th_star=45, flux_key=1, i_ang=[45])
        # self.sqrt = Sqrt(th_star=[45, 89])
        # self.line = Line(th_star=[45, 89])
        # self.pow_2 = Power(n=2, th_star=[45, 60])
        # self.pow_3 = Power(n=3, th_star=[45, 55])
        # self.pow_4 = Power(n=4, th_star=[45, 50])

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

        # if lum==0.1:
        #     # ax.set_ylim(0.003, 0.5)
        #     # ax.set_xlim(1.0, 20)

        # if lum==0.9:
        #     # ax.set_ylim(0.09, 2)
        #     # ax.set_xlim(1.0, 20)

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
        print("N  | L/L_Edd  | F_sl/F_Edd | T_eff (keV)  | f_c      | w")

        dumper(experiment, config, grid)
        config, grid = loader(experiment)
        ns = NeutronStar(config,grid)
        burst = ns.burst(mode='base', inter='be', lum=lum)
        counter = 1

        for shot in burst:
            print(f"{counter:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
            ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real * shot.dE_null/shot.E_null.unit/shot.B_real.unit/10**36, color='blue')


        # burst = ns.burst(mode='base', inter='wfc')
        # counter = 1
        # for shot in burst:
        #     if FLUX_REL[counter-1] == lum:
        #         print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
        #         ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real * shot.E_null/shot.E_null.unit/shot.B_real.unit/10**36, linestyle='dashed', color='blue')
        #         counter += 1
        #         break
        #     else:
        #         counter += 1
        #         continue

        burst = ns.burst(mode='model', inter='be', lum=lum)
        W_model = ns.surface_model.W_model
        grv_real = ns.surface_model.grv_real
        omega_kep = ns.param.omega_kep
        counter = 1

        for shot in burst:
            lambda_lum_be = shot.lambda_lum
            Flux_edd_real = ns.shot.Flux_edd_real
            flux = ns.shot.flux
            Flux_full = ns.shot.Flux_full
            Flux_full_sl = ns.shot.Flux_full_sl
            xi_lum = ns.shot.xi_lum
            Flux = ns.shot.Flux
            Flux_sl = ns.shot.Flux_sl
            print(f"{counter:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
            ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real * shot.dE_null/shot.E_null.unit/shot.B_real.unit/10**36, color='red')
            # print(shot.B_real)

        # burst = ns.burst(mode='model', inter='wfc')
        # counter = 1
        # for shot in burst:
        #     if FLUX_REL[counter-1] == lum:
        #         lambda_lum_wfc = shot.lambda_lum
        #         print(f"{counter:<3}| {shot.lum:.6f} | {shot.lum_sl:.6f}   | {shot.Epsilon_eff:.6f} | {shot.fc:.6f} | {shot.w:.6f} ")
        #         ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real * shot.E_null/shot.E_null.unit/shot.B_real.unit/10**36,linestyle='dashed', color='red')
        #         counter += 1
        #         break
        #     else:
        #         counter += 1
        #         continue
        # print(f"xi_lum = \n{xi_lum}")
        # print(f"flux = \n{flux}")
        # print(f"Flux_full_my = \n{Flux_edd_real * flux}")
        # print(f"Flux_full = \n{Flux_full}")
        # print(f"Flux_full_sl = \n{Flux_full_sl}")
        # print(f"Flux = {Flux}")
        # print(f"Flux_sl = {Flux_sl}")
        # print(f"flux_full = \n{flux_full}")
        # print(f"flux_full_int = {np.sum(flux_full)}")
        print(f"v_kep = {omega_kep / 2 / PI << u.Hz}")
        print(f"grv_sign = {np.round(grv_real.value[0,:]/np.abs(grv_real.value[0,:]))}")
        print(f"model = {np.round(W_model.value[0,n_theta//2::], decimals=2)}")
        print(f"lambda_be = {np.round(lambda_lum_be, decimals=2)}")
        print(f"xi_be = {np.round(Flux_sl/Flux, decimals=2)}")
        # print(f"lambda_wfc = {np.round(lambda_lum_wfc, decimals=2)}")

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


base_energy = BaseEnergy(
    name = 'test-energy', 
    grids = (21, 21, 500), 
    save = False, 
    show = True,
    lum = 0.05,
    w_func = 'const',
)

# print(base_energy)
# base_energy()
>>>>>>> main
