from atmons import *

import matplotlib.pyplot as plt

class BaseKappa(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const(w_func="vconst",  w_par = 1.0, th_star=[90])

    def do_experiment(self,
            w_func='const', th_star=45, w_par=1.0, 
            grids=(30, 30), 
            v_rot=600, i_ang=90,
            chem='s1',
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))
        print(w_func)
        # ax.set_ylim(0.3, 1.2)
        # ax.set_xlim(1., 20)
        ax.set_ylim(0.0, 1.0)
        # ax.set_xlim(1., 20)
        # ax.set_title("$g_{eff}(\\theta)$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        # ax.set_xlabel("$\\cos \\theta$")
        # ax.set_ylabel("$g_{eff}/g_0$")
        ax.set_title("$\\kappa(\\nu_{\\star})$ | $M = 1.5 M_{\\odot}$ | $R = 12km$", loc='center', fontsize=20)
        ax.set_xlabel("$\\nu_{\\star}, Hz$")
        ax.set_ylabel("$\\kappa$")
        ax.grid(True, which='minor')
        kappa, kappa_int, kappa_cr, nu_rel, nu = [], [], [], [], []
        # for k in range(0, 1100, 1):
        ns = build_Surface(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1],
        th_star=th_star, w_par = w_par,
        w_func='base',
        )
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # print(f"{v_rot}", surf.kappa_teory) 

        
        # l_cos_th = len(grid.cos_th[0, ::]) // 2

        kappa.append(surf.kappa_teory)
        kappa_int.append(surf.kappa_int)
        kappa_cr.append(surf.kappa_cr)
        nu.append(v_rot)
        nu_rel.append(body.v_cr.value)

        print(kappa, kappa_int, kappa_cr)
        ax.plot(nu, kappa, color='blue', label="$\\kappa_{ng}$")
        ax.plot(nu, kappa_int, color='red', label="$\\kappa_{int}$")
        ax.plot(nu, kappa_cr, color='green', label="$\\kappa_{cr}$")
        ax.plot(nu_rel, np.linspace(0,1,len(nu_rel)), color='black', linestyle='dashed', label="$\\nu_{crit}$")

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func=w_func, th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer

        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # # print(f"theta_max = {surf.th_star}")


        # l_cos_th = len(grid.cos_th[0, ::]) // 2
        # # print(body.Omega)
        # # print(sp_layer.omega_kep_local / body.omega_rot * v_rot)
        # ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='blue', label=f'spread layer: const')
        



        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func=w_func, th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # # print(sp_layer.kep_part)   
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)

        # l_cos_th = len(grid.cos_th[0, ::]) // 2
        # print(body.Omega)
        # print(np.max(surf.Omega_model))
        # ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='green', label=f'{w_func}')
        
        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='vpower-2', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # print(sp_layer.kep_part)   
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)

        # l_cos_th = len(grid.cos_th[0, ::]) // 2
        # # print(body.Omega)
        # # print(np.max(surf.omega_model/(2*PI)))
        # ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='green', label=f'spread layer: power-2')

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='vline', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # # print(sp_layer.kep_part)   
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)

        # l_cos_th = len(grid.cos_th[0, ::]) // 2
        # # print(body.Omega)
        # # print(np.max(surf.omega_model/(2*PI)))
        # ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='red', label=f'spread layer: line')

        # v_rot = 954.8
        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='base', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # # print(sp_layer.kep_part)   
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)

        # l_cos_th = len(grid.cos_th[0, ::]) // 2
        # print(np.max(surf.Omega_model))
        # ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='cyan', label=f'v_rot = {v_rot} Hz')


        # v_rot = 600
        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='base', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], label=f'v_rot = {v_rot} Hz')

        # v_rot = 800
        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='base', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], label=f'v_rot = {v_rot} Hz')

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='const', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], color='red', label='const')

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='vkconst', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], color='yellow', label='vkconst')

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='vconst', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], color='green', label='vconst')

        plt.legend()
        if save:
            if not os.path.isdir(f'graph/{experiment}'):
                    os.mkdir(f'graph/{experiment}')
            name = f'graph/{experiment}/'

            name += f"{w_func}"

            name += f"_th{th_star}"

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

        return 0

class BaseGeff(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const(w_func="vconst",  w_par = 1.0, th_star=[90])

    def do_experiment(self,
            w_func='const', th_star=45, w_par=1.0, 
            grids=(30, 30), 
            v_rot=600, i_ang=90,
            chem='s1',
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))
        print(w_func)
        ax.set_ylim(0.3, 1.2)
        # ax.set_xlim(1., 20)
        ax.set_title("$g_{eff}(\\tau)$ | $\\tau_{\\star} = " + str(th_star) + "^{\\circ}$ | $ \\kappa = \\kappa_{int}$", loc='center', fontsize=20)
        ax.set_xlabel("$\\sin \\tau$")
        ax.set_ylabel("$g_{eff}/g_0$")
        ax.grid(True, which='minor')


        ns = build_Surface(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1],
        th_star=th_star, w_par = w_par,
        w_func='base',
        )
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body
        psi = np.abs(90*DEG - grid.theta_init)
        
        l_cos_th = len(grid.cos_th[0, ::]) // 2
        ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='black', label=f'original')


        ns = build_Surface(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1],
        w_func=w_func, th_star=th_star, w_par = w_par,
        )
        sp_layer = ns.sp_layer

        grid = ns.grid
        surf = ns.surface
        body = ns.body
        psi = np.abs(90*DEG - grid.theta_init)
        print("surf", surf.kappa_max)  
        print("surf", surf.kappa_int)   
        print("surf", surf.kappa_cr)  

        l_cos_th = len(grid.cos_th[0, ::]) // 2
        # print(body.Omega)
        # print(sp_layer.omega_kep_local / body.omega_rot * v_rot)
        ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='blue', label='spread layer: $W_{\\infty}$')
        



        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func=w_func, th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # # print(sp_layer.kep_part)   
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)

        # l_cos_th = len(grid.cos_th[0, ::]) // 2
        # print(body.Omega)
        # print(np.max(surf.Omega_model))
        # ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='green', label=f'{w_func}')
        
        ns = build_Surface(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1],
        w_func='vpower-2', th_star=th_star, w_par = w_par,
        )
        sp_layer = ns.sp_layer
        # print(sp_layer.kep_part)   
        grid = ns.grid
        surf = ns.surface
        body = ns.body
        psi = np.abs(90*DEG - grid.theta_init)

        l_cos_th = len(grid.cos_th[0, ::]) // 2
        # print(body.Omega)
        # print(np.max(surf.omega_model/(2*PI)))
        ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='green', label=f'spread layer: $W_2$')

        ns = build_Surface(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1],
        w_func='vline', th_star=th_star, w_par = w_par,
        )
        sp_layer = ns.sp_layer
        # print(sp_layer.kep_part)   
        grid = ns.grid
        surf = ns.surface
        body = ns.body
        psi = np.abs(90*DEG - grid.theta_init)

        l_cos_th = len(grid.cos_th[0, ::]) // 2
        # print(body.Omega)
        # print(np.max(surf.omega_model/(2*PI)))
        ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='red', label=f'spread layer: $W_1$')

        # v_rot = 954.8
        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='base', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # # print(sp_layer.kep_part)   
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)

        # l_cos_th = len(grid.cos_th[0, ::]) // 2
        # print(np.max(surf.Omega_model))
        # ax.plot(grid.cos_th[0,:l_cos_th-1:], surf.g_th[0,:l_cos_th-1:], color='cyan', label=f'v_rot = {v_rot} Hz')


        # v_rot = 600
        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='base', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], label=f'v_rot = {v_rot} Hz')

        # v_rot = 800
        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='base', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], label=f'v_rot = {v_rot} Hz')

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='const', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], color='red', label='const')

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='vkconst', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], color='yellow', label='vkconst')

        # ns = build_Surface(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1],
        # w_func='vconst', th_star=th_star, w_par = w_par,
        # )
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body
        # psi = np.abs(90*DEG - grid.theta_init)
        # ax.plot(psi[0], surf.g_th[0], color='green', label='vconst')

        plt.legend()
        if save:
            if not os.path.isdir(f'graph/{experiment}'):
                    os.mkdir(f'graph/{experiment}')
            name = f'graph/{experiment}/'

            name += f"{w_func}"

            name += f"_th{th_star}"

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

        return 0

class BaseSpectra(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        #self.const = Const(w_func="vkpower-2.0", flux_key = 'rel',  w_par = 0.8, th_star=[25])
        self.vkline = Line(w_func="base", flux_key = 'rel')
        # self.vkpow2 = Power(w_func="vkpower-2", flux_key = 'rel', i_ang=[45,75,90])
        # self.vkconst = Const(w_func="vkconst", flux_key = 'rel', i_ang=[45,75,90])
        # self.sqrt = Sqrt(w_func="sqrt", flux_key=1, w_par = 0.9, th_star=[45])
        # self.line = Line(w_func="line", flux_key=1, w_par = 1.0, th_star=[45])
        # self.pow_2 = Power(n=2, w_func="power-2", flux_key=1, w_par = 0.9, th_star=[45])
        # self.pow_3 = Power(n=3, th_star=[45, 55])
        # self.pow_4 = Power(n=4, th_star=[45, 50])
        # self.pow_1t3 = Power(n=0.01, w_func="power-0.01", w_par=0.9, th_star=[45])

    def do_experiment(self,
            w_func='const', th_star=45, w_par=1.0, 
            grids=(30, 30, 500), 
            v_rot=600, i_ang=90, lum=0.1, 
            chem='s1', flux_key='rel',
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))
        print(w_func)
        if w_func == 'vkpower-2':
            ax.set_title("$ v(\\theta) \\sim " + "\\frac{1 - (\\theta/\\theta_{\\star})^2}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        elif w_func == 'vkline':
            ax.set_title("$ v(\\theta) \\sim " + "\\frac{1 - (\\theta/\\theta_{\\star})}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        elif w_func == 'vkconst':
            ax.set_title("$ v(\\theta) \\sim " + "\\frac{1}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        else: 
            ax.set_title("$B(\\varepsilon), f = 0.1, \\nu = 600$, chem=s1", loc='center', fontsize=20)

        ax.set_xlabel("$\\varepsilon, keV$")
        ax.set_ylabel("$B(\\varepsilon), 10^{36} erg s^{-1} keV^{-1} sr^{-1} $")
        ax.grid(True, which='minor')

        # if lum==0.1:
        ax.set_ylim(0.003, 0.5)
        ax.set_xlim(1.0, 20)

        # if lum==0.9:
        #     ax.set_ylim(0.09, 2)
        #     ax.set_xlim(1.0, 20)

        # else:
        #     ax.set_ylim(0.003, 2)
        #     ax.set_xlim(1.0, 20)


        # config['rel'] = False#True
        # config['m_ns'] = 1.5#1.519
        # config['r_ns'] = 12#15.48
        # config['v_rot'] = 600#700
        # config['rel'] = True

        # config['m_ns'] = 1.519
        # config['r_ns'] = 15.48
        # config['v_rot'] = v_rot

        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")

        ns = build_Neutron_Star(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        th_star=th_star, w_par = w_par,
        flux_key='rel', lum=lum, 
        spec_key = 'wfc', w_func='base',
        )

        shot = ns.atmosphere
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body
        print(f"{'wfc':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red', label='diluted BB') #linestyle='dashed'
        B_dbb = rho_rad(shot.E_null, shot.Epsilon_eff* shot.fc , shot.w) * surf.area_real 
        ax.loglog(shot.E_null/shot.E_null.unit, B_dbb/B_dbb.unit/10**36, color='red', linestyle='dashed', label='diluted BB') #


        ns = build_Neutron_Star(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        th_star=th_star, w_par = w_par,
        lum=lum, flux_key=flux_key,
        spec_key = 'be', w_func='base',
        )
        shot = ns.atmosphere
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body

        counter = 1
        th = grid.theta.value[0][len(surf.g_th[0])//2::]
        g_th = surf.g_th.value[0][len(surf.g_th[0])//2::]
        flux_base = shot.flux
        print(f"{'s1 ':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue', label='model spectra')
        
        B_dbb = rho_rad(shot.E_null, shot.Epsilon_eff* shot.fc , shot.w) * surf.area_real 
        ax.loglog(shot.E_null/shot.E_null.unit, B_dbb/B_dbb.unit/10**36, color='blue', linestyle='dashed', label='model spectra') #

        # ns = build_Neutron_Star(
        # v_rot=v_rot, i_ang=i_ang, chem='s001',
        # N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        # th_star=th_star, w_par = w_par,
        # lum=lum, flux_key=flux_key,
        # spec_key = 'be', w_func='base',
        # )
        # shot = ns.atmosphere
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body


        # th = grid.theta.value[0][len(surf.g_th[0])//2::]
        # g_th = surf.g_th.value[0][len(surf.g_th[0])//2::]
        # flux_base = shot.flux
        # print(f"{"s0 ":<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        # ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='green', label='without Fe')


        # ns = build_Neutron_Star(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        # w_func=w_func, th_star=th_star, w_par = w_par,
        # lum=lum, spec_key = 'be', flux_key='sl',
        # )
        # shot = ns.atmosphere
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body

        # W_model = surf.W_model
        # grv_real = surf.grv
        # omega_kep = body.omega_kep
        # th_m = grid.theta.value[0][len(surf.g_th[0])//2::]
        # g_th_m = surf.g_th.value[0][len(surf.g_th[0])//2::]
        # B_real = None
        # counter = 1

        # flux = shot.flux
        # W_model = surf.W_model
        # print(f"{"sl ":<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        # ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red', label='spread-layer')
        # B_real = shot.B_real
        
        # l_flux = len(flux[0, ::]) // 2
        # w_mod = len(W_model[0, ::]) // 2
        # print(flux_base[0, l_flux::])
        # print(flux[0, l_flux::])
        # # print(W_model[0, w_mod::])


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
        plt.legend()

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

        # th_zero = 0
        # f_th = g_th_m - g_th
        # for k in range(len(f_th)):
        #     if f_th[k]>0:
        #         y1, y2 = f_th[k-1], f_th[k]
        #         x1, x2 = th[k-1], th[k]
        #         a = (x2 - x1) / (y2 - y1)
        #         b = y1 - a*x1
        #         th_zero = - b / a
        #         break
    

        # return th_zero, th, g_th, g_th_m

class BaseModel(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        #self.const = Const(w_func="vkpower-2.0", flux_key = 'rel',  w_par = 0.8, th_star=[25])
        self.const = Const()
        # self.vkpow2 = Power(w_func="vkpower-2", flux_key = 'rel', i_ang=[45,75,90])
        # self.vkconst = Const(w_func="vkconst", flux_key = 'rel', i_ang=[45,75,90])

    def do_experiment(self,
            w_func='const', th_star=45, w_par=1.0, 
            grids=(30, 30, 500), 
            v_rot=600, i_ang=90, lum=0.1, 
            chem='s1', flux_key='rel',
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))

        if w_func == 'vkpower-2':
            ax.set_title("$ v(\\theta) \\sim " + "\\frac{1 - (\\theta/\\theta_{\\star})^2}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        elif w_func == 'vkline':
            ax.set_title("$ v(\\theta) \\sim " + "\\frac{1 - (\\theta/\\theta_{\\star})}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        elif w_func == 'vkconst':
            ax.set_title("$ v(\\theta) \\sim " + "\\frac{1}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        else:
            ax.set_title("$\\kappa = \\kappa_{int}$" + " | $\\tau_{\\star} = " + str(th_star) + "^{\\circ} " + "$ | $i = " + str(i_ang) + "^{\\circ}$", loc='center', fontsize=20)

            #ax.set_title("$B(\\varepsilon), f = 0.1, \\nu = 600$, chem=s1", loc='center', fontsize=20)

        ax.set_xlabel("$\\varepsilon, keV$")
        ax.set_ylabel("$B(\\varepsilon), 10^{36} erg s^{-1} keV^{-1} sr^{-1} $")
        ax.grid(True, which='minor')
        ax.grid(True, which='minor')

        ax.set_xlim(1.0, 20)
        if lum==0.1:
            ax.set_ylim(0.003, 0.5)

        if lum==0.9:
            ax.set_ylim(0.09, 2)


        # else:
        #     ax.set_ylim(0.003, 2)
        #     ax.set_xlim(1.0, 20)


        # config['rel'] = False#True
        # config['m_ns'] = 1.5#1.519
        # config['r_ns'] = 12#15.48
        # config['v_rot'] = 600#700
        # config['rel'] = True

        # config['m_ns'] = 1.519
        # config['r_ns'] = 15.48
        # config['v_rot'] = v_rot

        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")

        ns = build_Neutron_Star(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        th_star=th_star, w_par = w_par,
        flux_key=flux_key, lum=lum, 
        spec_key = 'wfc', w_func='base',
        )

        shot = ns.atmosphere
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body
        print(f"{'wfc':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, linestyle='dashed', color='black', label='diluted BB')


        ns = build_Neutron_Star(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        th_star=th_star, w_par = w_par,
        lum=lum, flux_key=flux_key,
        spec_key = 'be', w_func='base',
        )
        shot = ns.atmosphere
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body

        counter = 1
        th = grid.theta.value[0][len(surf.g_th[0])//2::]
        g_th = surf.g_th.value[0][len(surf.g_th[0])//2::]
        flux_base = shot.flux
        print(f"{'s1 ':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='black', label='original')

        ns = build_Neutron_Star(
        v_rot=v_rot, i_ang=i_ang, chem='s001',
        N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        th_star=th_star, w_par = w_par,
        lum=lum, flux_key=flux_key,
        spec_key = 'be', w_func='base',
        )
        shot = ns.atmosphere
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body


        th = grid.theta.value[0][len(surf.g_th[0])//2::]
        g_th = surf.g_th.value[0][len(surf.g_th[0])//2::]
        flux_base = shot.flux
        print(f"{'s0 ':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='green', label='original: $0.01 z_{\\odot}$')

        w_func = 'vline'
        ns = build_Neutron_Star(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        w_func=w_func, th_star=th_star, w_par = w_par,
        lum=lum, spec_key = 'be', flux_key='sl',
        )
        shot = ns.atmosphere
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body

        W_model = surf.W_model
        grv_real = surf.grv
        omega_kep = body.omega_kep
        th_m = grid.theta.value[0][len(surf.g_th[0])//2::]
        g_th_m = surf.g_th.value[0][len(surf.g_th[0])//2::]
        B_real = None
        counter = 1

        flux = shot.flux
        W_model = surf.W_model
        print(f"{'sll':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='red', label='spread-layer: $W_1$')
        B_real = shot.B_real
        
        l_flux = len(flux[0, ::]) // 2
        w_mod = len(W_model[0, ::]) // 2
        # print(flux_base[0, l_flux::])
        print(flux[0, l_flux::])

        # w_func = 'vkpower-2'
        # ns = build_Neutron_Star(
        # v_rot=v_rot, i_ang=i_ang, chem=chem,
        # N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        # w_func=w_func, th_star=th_star, w_par = w_par,
        # lum=lum, spec_key = 'be', flux_key='sl',
        # )
        # shot = ns.atmosphere
        # sp_layer = ns.sp_layer
        # grid = ns.grid
        # surf = ns.surface
        # body = ns.body

        # W_model = surf.W_model
        # grv_real = surf.grv
        # omega_kep = body.omega_kep
        # th_m = grid.theta.value[0][len(surf.g_th[0])//2::]
        # g_th_m = surf.g_th.value[0][len(surf.g_th[0])//2::]
        # B_real = None
        # counter = 1

        # flux = shot.flux
        # W_model = surf.W_model
        # print(f"{"slp":<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        # ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue', label='spread-layer: pow-2')
        # B_real = shot.B_real
        
        # l_flux = len(flux[0, ::]) // 2
        # w_mod = len(W_model[0, ::]) // 2
        # print(flux_base[0, l_flux::])
        # print(flux[0, l_flux::])

        w_func = 'vconst'
        ns = build_Neutron_Star(
        v_rot=v_rot, i_ang=i_ang, chem=chem,
        N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        w_func=w_func, th_star=th_star, w_par = w_par,
        lum=lum, spec_key = 'be', flux_key='sl',
        )
        shot = ns.atmosphere
        sp_layer = ns.sp_layer
        grid = ns.grid
        surf = ns.surface
        body = ns.body

        W_model = surf.W_model
        grv_real = surf.grv
        omega_kep = body.omega_kep
        th_m = grid.theta.value[0][len(surf.g_th[0])//2::]
        g_th_m = surf.g_th.value[0][len(surf.g_th[0])//2::]
        B_real = None
        counter = 1

        flux = shot.flux
        W_model = surf.W_model
        print(f"{'slc':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue', label='spread-layer: $W_{\\infty}$')
        B_real = shot.B_real
        
        l_flux = len(flux[0, ::]) // 2
        w_mod = len(W_model[0, ::]) // 2
        print(flux_base[0, l_flux::])
        print(flux[0, l_flux::])

        # print(W_model[0, w_mod::])


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
        plt.legend()

        if save:
            if not os.path.isdir(f'graph/{experiment}'):
                    os.mkdir(f'graph/{experiment}')
            name = f'graph/{experiment}/'

            name += f"model"

            name += f"_{w_par}vkep"

            name += f"_th{th_star}"

            name += f"_i{i_ang}"

            # if i_ang != 45:
            #     name += f"_i{i_ang}" 
            # if lum != 0.1:
            #     name += f"_l{lum}"
            # if chem != 's1':
            #     name += f"_{chem}"
            # if v_rot != 700:
            #     name += f"_{v_rot}hz"

            name += '.pdf'

            plt.savefig(name)

        if show:
            plt.show()

        th_zero = 0
        f_th = g_th_m - g_th
        for k in range(len(f_th)):
            if f_th[k]>0:
                y1, y2 = f_th[k-1], f_th[k]
                x1, x2 = th[k-1], th[k]
                a = (x2 - x1) / (y2 - y1)
                b = y1 - a*x1
                th_zero = - b / a
                break
    

        return th_zero, th, g_th, g_th_m 

class BaseFc(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        #self.const = Const(w_func="vkpower-2.0", flux_key = 'rel',  w_par = 0.8, th_star=[25])
        self.const = Const()
        # self.vkpow2 = Power(w_func="vkpower-2", flux_key = 'rel', i_ang=[45,75,90])
        # self.vkconst = Const(w_func="vkconst", flux_key = 'rel', i_ang=[45,75,90])

    def do_experiment(self,
            w_func='const', th_star=45, w_par=1.0, 
            grids=(30, 30, 500), 
            v_rot=600, i_ang=90, lum=0.1, 
            chem='s1', flux_key='rel',
            show=True, save=False, experiment='test'
        ):
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))
        print(w_func)
        # if w_func == 'vkpower-2':
        #     ax.set_title("$ v(\\theta) \\sim " + "\\frac{1 - (\\theta/\\theta_{\\star})^2}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        # elif w_func == 'vkline':
        #     ax.set_title("$ v(\\theta) \\sim " + "\\frac{1 - (\\theta/\\theta_{\\star})}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        # elif w_func == 'vkconst':
        #     ax.set_title("$ v(\\theta) \\sim " + "\\frac{1}{(r/r_{eq})}" + "$ | $ v(0) = " + str(w_par) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(th_star) + "^{\\circ}$ ", loc='center', fontsize=20)
        # else:

        w_par_name = w_par
        
        ax.set_title("$\\kappa = \\kappa_{int}$" + " | $\\tau_{\\star} = " + str(th_star) + "^{\\circ} " + "$ | $i = " + str(i_ang) + "^{\\circ}$", loc='center', fontsize=20)

        ax.set_xlabel("$L/L_{Edd}$")
        ax.set_ylabel("$f_c$")
        ax.grid(True, which='minor')

        # # if lum==0.1:
        # ax.set_ylim(0.003, 0.5)
        # ax.set_xlim(1.0, 20)

        # if lum==0.9:
        #     ax.set_ylim(0.09, 2)
        #     ax.set_xlim(1.0, 20)

        # Fc, Lum = [], []
        # print("\nwfc\n")
        # print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        # for k in range(len(FLUX_REL)):
        #     lum = FLUX_REL[k]
        #     ns = build_Neutron_Star(
        #     v_rot=v_rot, i_ang=i_ang, chem=chem,
        #     N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        #     th_star=th_star, w_par = w_par,
        #     flux_key=flux_key, lum=lum, 
        #     spec_key = 'wfc', w_func='base',
        #     )
        #     shot = ns.atmosphere
        #     sp_layer = ns.sp_layer
        #     print(f"{k+1:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        #     Lum.append(shot.lum)
        #     Fc.append(shot.fc)
        
        # ax.plot(Lum, Fc, linestyle='-', color='blue', label='diluted BB')

        Fc, Lum = [], []
        print("\ns1\n")
        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        for k in range(len(FLUX_REL)):
            lum = FLUX_REL[k]
            ns = build_Neutron_Star(
            v_rot=v_rot, i_ang=i_ang, chem=chem,
            N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
            th_star=th_star, w_par = w_par,
            flux_key=flux_key, lum=lum, 
            spec_key = 'be', w_func='base',
            )
            shot = ns.atmosphere
            sp_layer = ns.sp_layer
            print(f"{k+1:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
            Lum.append(shot.lum)
            Fc.append(shot.fc)
        l_flux = len(shot.flux[0, ::]) // 2
        print(shot.flux[0, l_flux::])  
        ax.plot(Lum, Fc, linestyle='-', color='black', label='original')
        
        Fc, Lum = [], []
        print("\ns0\n")
        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        for k in range(len(FLUX_REL)):
            lum = FLUX_REL[k]
            ns = build_Neutron_Star(
            v_rot=v_rot, i_ang=i_ang, chem='s001',
            N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
            th_star=th_star, w_par = w_par,
            flux_key=flux_key, lum=lum, 
            spec_key = 'be', w_func='base',
            )
            shot = ns.atmosphere
            sp_layer = ns.sp_layer
            print(f"{k+1:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
            Lum.append(shot.lum)
            Fc.append(shot.fc)
        l_flux = len(shot.flux[0, ::]) // 2
        print(shot.flux[0, l_flux::])   
        
        ax.plot(Lum, Fc, linestyle='-', color='green', label='original: $0.01z_{\\odot}$')
        
        Fc, Lum = [], []
        print("\nsl - line\n")
        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        for k in range(len(FLUX_REL)):
            lum = FLUX_REL[k]
            ns = build_Neutron_Star(
            v_rot=v_rot, i_ang=i_ang, chem=chem,
            N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
            th_star=th_star, w_par = w_par,
            flux_key='sl', lum=lum, 
            spec_key = 'be', w_func='vline',
            )
            shot = ns.atmosphere
            sp_layer = ns.sp_layer
            print(f"{k+1:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
            Lum.append(shot.lum)
            Fc.append(shot.fc)
        l_flux = len(shot.flux[0, ::]) // 2
        print(shot.flux[0, l_flux::])
        ax.plot(Lum, Fc, linestyle='-', color='red', label='spread-layer: $W_1$')

        # Fc, Lum = [], []
        # print("\nsl - pow-2\n")
        # print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        # for k in range(len(FLUX_REL)):
        #     lum = FLUX_REL[k]
        #     ns = build_Neutron_Star(
        #     v_rot=v_rot, i_ang=i_ang, chem=chem,
        #     N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
        #     th_star=th_star, w_par = w_par,
        #     flux_key='sl', lum=lum, 
        #     spec_key = 'be', w_func='vpower-2',
        #     )
        #     shot = ns.atmosphere
        #     sp_layer = ns.sp_layer
        #     print(f"{k+1:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        #     Lum.append(shot.lum)
        #     Fc.append(shot.fc)
        # l_flux = len(shot.flux[0, ::]) // 2
        # print(shot.flux[0, l_flux::])
        # ax.plot(Lum, Fc, linestyle='-', color='yellow', label='spread-layer: $W_2$')

        Fc, Lum = [], []
        print("\nsl - const\n")
        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        for k in range(len(FLUX_REL)):
            lum = FLUX_REL[k]
            ns = build_Neutron_Star(
            v_rot=v_rot, i_ang=i_ang, chem=chem,
            N_ph=grids[0], N_th=grids[1], N_nu=grids[2],
            th_star=th_star, w_par = w_par,
            flux_key='sl', lum=lum, 
            spec_key = 'be', w_func='vconst',
            )
            shot = ns.atmosphere
            sp_layer = ns.sp_layer
            print(f"{k+1:<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
            Lum.append(shot.lum)
            Fc.append(shot.fc)
        l_flux = len(shot.flux[0, ::]) // 2
        print(shot.flux[0, l_flux::])
        ax.plot(Lum, Fc, linestyle='-', color='blue', label='spread-layer: $W_{\\infty}$')
        plt.legend()

        if save:
            if not os.path.isdir(f'graph/{experiment}'):
                    os.mkdir(f'graph/{experiment}')
            name = f'graph/{experiment}/'

            name += f"fc"

            name += f"_i{i_ang}"

            name += f"_{w_par_name}k"

            name += f"_th{th_star}"


            name += '.pdf'

            plt.savefig(name)

        if show:
            plt.show()

        return 0
    
class DiscussModel(Experiment):
    def __init__(self, name, **param):
        super().__init__(name, **param)
        self.const = Const()

    def do_experiment(self, experiment='test'):
        config = {
            'v_rot': 600,
            'm_ns': 1.4,
            'r_ns': 12,
            'i_ang': 90,
            'chem': 's1',
            'spec_key': 'be',
            'th_star': 60,
            'w_func': 'const',
            'w_par': 0.932,
        }
        ns = build_Neutron_Star(**config)
        config['w_func'] = 'base'
        ns_base = build_Neutron_Star(**config)
    
        shot = ns.atmosphere
        surf = ns.surface
        sp_layer = ns.sp_layer

        shot_base = ns_base.atmosphere
        
        g_th = surf.log_g.value[0][len(surf.g_th[0])//2::]
        flux = shot.flux.value[0][len(surf.g_th[0])//2::]
        plt.style.use('seaborn-v0_8-whitegrid')
        _, ax = plt.subplots(figsize=(7,7))
        # line
        #ax.set_title("$ v(\\theta) \\sim " + "1 - \\theta/\\theta_{\\star}" + "$ | $ v(0) = " + str(config['w_par']) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(config['th_star']) + "^{\\circ}$ ", loc='center', fontsize=20)
        # const
        ax.set_title("$ v(\\theta) \\sim " + "const" + "$ | $ v(0) = " + str(config['w_par']) + "\\cdot v_{kep}" + "$ | $\\theta_{\\star} = " + str(config['th_star']) + "^{\\circ}$ ", loc='center', fontsize=20)
        ax.set_xlabel("$\\varepsilon, keV$")
        ax.set_ylabel("$B(\\varepsilon), 10^{36} erg s^{-1} keV^{-1} sr^{-1} $")
        ax.grid(True, which='minor')
        ax.set_ylim(0.003, 0.5)
        ax.set_xlim(1.0, 20)

        print("N  | L/L_Edd  | L_sl/L_Edd | f_c      | w        | kep/rot")
        print(f"{'be':<3}| {shot.lum:.6f} | {shot.xi_sl:.6f}   | {shot.fc:.6f} | {shot.w:.6f} | {sp_layer.rel_omega:.6f} ")
        print(g_th)
        print(flux)
        
        ax.loglog(shot.E_null/shot.E_null.unit, shot.B_real/shot.B_real.unit/10**36, color='blue', label='spread layer')
        ax.loglog(shot_base.E_null/shot_base.E_null.unit, shot_base.B_real/shot_base.B_real.unit/10**36, color='red', label='base spectrum')
        ax.legend()
        plt.show()
        return 0


model = DiscussModel(name = 'discussion')

model()

spectra_rel = BaseSpectra(
    name = 'ended-spectra-new', 
    grids = (30, 30, 500), 
    w_par = 0.5,
    th_star = 50,
    save = True, 
    show = True,
    flux_key = 'rel',
    lum = 0.1,
    w_func='const',
)

g_eff_rel = BaseGeff(
    name = 'rep-geff-new', 
    grids = (10, 10000), 
    save = True, 
    show = True,
    w_func='const',
)

model_rel = BaseModel(
    name = 'discussion', 
    grids = (20, 20, 500), 
    i_ang = 90,
    w_par = 1.0,
    th_star = 45,
    save = True, 
    show = True,
    flux_key = 'rel',
    lum = 0.9,
    w_func='const',
)

# model_rel()


# model_fc = BaseFc(
#     name = 'rep-model-fc', 
#     grids = (10, 10, 50), 
#     i_ang = 60,
#     w_par = 1.0,
#     th_star = 30,
#     save = True, 
#     show = True,
# )
# model_fc()

# g_eff_rel()

# for lum in [0.9]:
#     for i in [90]:
#         for th in [90]:
#             print(f"\n Check ({i}i, {th}th)\n")
#             model_fc = BaseModel(
#                 name = 'rep-model-spec', 
#                 grids = (30, 30, 500), 
#                 i_ang = i,
#                 w_par = 1.0,
#                 th_star = th,
#                 lum=lum,
#                 save = True, 
#                 show = False,
#             )
#             model_fc()

# for i in [45]:
#     for th in [60]:
#         print(f"\n Check ({i}i, {th}th)\n")
#         model_fc = BaseModel(
#             name = 'last', 
#             grids = (20, 20, 500), 
#             i_ang = i,
#             w_par = 1.0,
#             th_star = th,
#             save = True, 
#             show = False,
#         )
#         model_fc()
#         model_fc = BaseFc(
#             name = 'last', 
#             grids = (20, 20, 500), 
#             i_ang = i,
#             w_par = 1.0,
#             th_star = th,
#             save = True, 
#             show = False,
#         )
#         model_fc()
