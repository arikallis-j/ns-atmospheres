from time import time
from .classes import *
from .nstar import *
from .graph import *

name = "J0000+0000"
system = "base"
chem = "s1"
rel = False
fc_key = 1 
flux_key = 1
w_key = "null"

m_ns = 1.5
r_ns = 12.0
v_rot = 700.0
i_ang = 60.0
lum_0 = 0.1

n_phi = 120
n_theta = 60
n_nu = 500


class Calculator:
    def test(self):
        lum_k = FLUX_REL[0]
        NS = NeurtonStarNew(name=name, sys=system, chem=chem, rel=rel,
                        fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                        m_ns=m_ns, r_ns=r_ns, v_rot=v_rot, i_ang=i_ang, lum=lum_k)
        
        # NS.calc(n_phi=130, n_theta=70, n_nu=n_nu)
        # print(f"lum:{NS.lum} energy: {NS.Epsilon_eff}")
        # print(f"w:{NS.w} fc: {NS.fc}")
        NS.init_paramters()
        # print(NS.NS_Keys)
        # print(NS.NS_Init)
        # print(NS.NS_Phis)
        # print(NS.NS_Chem)
        # print(NS.NS_wbTc)
        # print(NS.NS_Phot)
        # print(NS.NS_Rot)
        # print(NS.NS_Rel)
        # print(NS.NS_Met)
        # print(NS.NS_Met)
        # print(NS.NS_Rad)

        
        start = time()
        NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
        mid = time()
        end = time()
        print(f'{mid-start:.4f}, {end-mid:.4f}, {end-start:.4f}')
        return "Tests ended"
    


    def print(self, name, kind="base"):
        return f"Printing of NeurtonStar: {name} with kind={kind}"
    
    def build(self, name):
        NS = NeurtonStar(name)
        print(f"Initicialisation of NeurtonStar: {NS.name} with n_phi={NS.n_phi}, n_theta={NS.n_theta}")
        print(f"Calculation of NeurtonStar: {name}")
        return f"Building done. "
    
    def init(self, name):
        NS = NeurtonStar(name)
        return f"Initicialisation of NeurtonStar: {NS.name} with n_phi={NS.n_phi}, n_theta={NS.n_theta}"
    
    def calc(self, name):
        return f"Calculation of NeurtonStar: {name}"
    
    def lines(self):
        I_arr = [90, 60, 40, 20, 10, 5]
        Colors = ['black', 'orange', 'violet', 'green', 'blue', 'red']
        xlabel = "Relative energy $E/(1-u)^(1/2)$"
        ylabel = "Normalized flux"
        N, M, k = 5, 20, 0#60, 1000, 0
        rel = True
        m_cor = 1.5
        r_eq = 14
        v_rot = 600
        n_phi = 2*N
        n_theta = N
        n_nu = M
        lum_k = FLUX_REL[k]
        g_database = GraphDataBase()
        for k in range(6):
            print(f"The {k} iteration...")
            NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
                            fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                            m_ns=m_cor, r_ns=r_eq, v_rot=v_rot, i_ang=I_arr[k], lum=lum_k)
            NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
            factor = sqrt(1 - NS.R_sch/NS.R_eq)
            Norm_energy = np.array(NS.E) / factor
            Norm_flux = np.array(NS.B_real) / max(NS.B_real)
            #draw_graph(Norm_energy, Norm_flux, color=Colors[k], xrange=(0.7, 1.3), title="$L/L_{edd} = 0.1$", xlabel=xlabel, ylabel=ylabel)
            data = GraphData(Norm_energy, Norm_flux)
            g_database.DataBase.append(data)

        BasePalette = GraphPalette(['red', 'blue', 'green'])
        LinePalette = GraphPalette(['black', 'orange', 'violet', 'green', 'blue', 'red'])

        g_size = GraphSize()
        g_style = GraphStyle(LinePalette)
        graph = Graph(g_database, g_size, g_style)
        graph.draw_picture(save=True)

    def spectra(self):
        Colors = ['red', 'blue', 'green', 'cyan']
        BasePalette = GraphPalette(['red', 'blue', 'green'])
        I_arr = [0, 45, 90]
        N, M = 60, 1000
        rel = True
        m_cor = 1.5
        r_eq = 15.48
        v_rot = 700
        n_phi = 2*N
        n_theta = N
        n_nu = M
        lum_k = 0.1
        xlabel = "Photon energy (keV)"
        ylabel = "$L_E (10^{36} erg^{-1} keV^{-1} sr^{-1})$"
        title="$L/L_{edd} = " + f"{lum_k}$"
        title_ascii = "L/L_{edd} " + f"{lum_k}"
        ylabel_ascii = "L_E, 10^36"
        file = "sp_01" #"sp90"
        DataBaseOld = [[],[],[],[]]
        with open(f'ns_atm/spectra/{file}.dat', 'r') as f:
            for line in f:
                arr_line = line.split(" ")[1::]
                arr_line[-1] = arr_line[-1].split("\n")[0]
                for k in range(len(arr_line)):
                    list_line = list(arr_line[k])
                    list_norm = ""
                    for j in range(len(list_line)):
                        if list_line[j] == 'D':
                            list_norm += 'e'
                        else:
                            list_norm += list_line[j]

                    arr_line[k] = list_norm
                    arr_line[k] = float(arr_line[k])
                    if k!=0:
                        arr_line[k] = arr_line[k]/10**36
                    DataBaseOld[k].append(arr_line[k])

        g_database_old = GraphDataBase(tex_mode=True,filename="old", latex_labels=(title, xlabel, ylabel), ascii_labels=(title_ascii, xlabel, ylabel_ascii))
        
        for n in range(3):
            data = GraphData(DataBaseOld[0], DataBaseOld[n+1], ascii_label=f"i = {I_arr[n]} deg")
            g_database_old.DataBase.append(data)
        g_size = GraphSize(x_range=(1,20), y_range=(0.003,0.4))
        g_style = GraphStyle(BasePalette, loglog=True)
        graph = Graph(g_database_old, g_size, g_style)
        graph.draw_picture(save=True, draw=False)   
        graph.draw_table() 

        g_database = GraphDataBase(tex_mode=True, latex_labels=(title, xlabel, ylabel), ascii_labels=(title_ascii, xlabel, ylabel_ascii))
        for k in range(3):
            NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
                            fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                            m_ns=m_cor, r_ns=r_eq, v_rot=v_rot, i_ang=I_arr[k], lum=lum_k)
            NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
            Energy = NS.E
            Flux = NS.B_real / 10**36
            data = GraphData(Energy, Flux, ascii_label=f"i = {I_arr[k]} deg")
            g_database.DataBase.append(data)

        g_size = GraphSize(x_range=(1,20), y_range=(0.003,0.4))
        g_style = GraphStyle(BasePalette, loglog=True)
        graph = Graph(g_database, g_size, g_style)
        graph.draw_picture(save=True, draw=False)   
        graph.draw_table() 

        g_database = GraphDataBase(tex_mode=True, filename="diff")
        for k in range(3):
            NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
                            fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                            m_ns=m_cor, r_ns=r_eq, v_rot=v_rot, i_ang=I_arr[k], lum=lum_k)
            NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
            Energy = NS.E
            # LogNew = log(NS.B_real / 10**36 + 0.000000001)
            # LogOld = log(np.array(DataBaseOld[k+1]) + 0.000000001)
            # Flux = 100 * (LogOld - LogNew) / (LogOld + 0.000000001)
            Old = np.array(DataBaseOld[k+1])
            New = np.array(NS.B_real) / 10**36
            Flux = []
            for l in range(len(New)):
                #koeff = -abs(log(New[l]) - log(Old[l])) / log(Old[l]) 
                koeff = abs(New[l] - Old[l]) / (Old[l]) 
                Flux.append(100 * koeff)
            data = GraphData(Energy, Flux, ascii_label=f"i = {I_arr[k]} deg")
            g_database.DataBase.append(data)

        g_size = GraphSize()
        g_style = GraphStyle(BasePalette)
        graph = Graph(g_database, g_size, g_style)
        graph.draw_picture(save=True, draw=False)   
        graph.draw_table() 


        #     draw_graph(Energy, Flux, color=Colors[k], loglog=True, xrange=(1,20), yrange=(0.003,0.4), title="$L/L_{edd} = 0.1$", xlabel=xlabel, ylabel=ylabel)
        # show_graph()

        # lum_k = FLUX_REL[13]
        # print(lum_k)
        # for k in range(2):
        #     NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
        #                     fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
        #                     m_ns=m_cor, r_ns=r_eq, v_rot=v_rot, i_ang=I_arr[k], lum=lum_k)
        #     NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
        #     Energy = np.array(NS.E)
        #     Flux = np.array(NS.B_real) / 10**36
        #     draw_graph(Energy, Flux, color=Colors[k], loglog=True, xrange=(1,20), yrange=(0.1, 2), title="$L/L_{edd} = 0.9$", xlabel=xlabel, ylabel=ylabel)
        # show_graph()

    def lfc(self, name=name, sys=system, chem=chem, rel=rel,
                  fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                  m_ns=m_ns, r_ns=r_ns, v_rot=v_rot, i_ang=i_ang, lum=lum_0,
                  n_phi=n_phi, n_theta=n_theta, n_nu=n_nu):
        start = time()
        for k in range(0, N_MODEL-2):
            lum_k = FLUX_REL[k]
            NS = NeurtonStarNew(name=name, sys=system, chem=chem, rel=rel,
                             fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                             m_ns=m_ns, r_ns=r_ns, v_rot=v_rot, i_ang=i_ang, lum=lum_k)
            NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
            if k==0:
                print("\n  surf/R_eq^2     rpr             R_eq:")
                print(f"  {sci(NS.surf / NS.R_eq**2, ndig=6)}     {sci(NS.R_pr, ndig=6)}     {sci(NS.R_eq, ndig=6)}")
                print("\n         L/L_Edd         T_eff (keV)     f_c             w")

            print(f"  {' ' if k<10 else ''}{k+1}     {sci(NS.lum)}      {sci(NS.Epsilon_eff)}      {sci(NS.fc)}      {sci(NS.w)}")
        mid = time()
        for k in range(0, N_MODEL-2):
            lum_k = FLUX_REL[k]
            NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
                             fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                             m_ns=m_ns, r_ns=r_ns, v_rot=v_rot, i_ang=i_ang, lum=lum_k)
            NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
            if k==0:
                print("\n  surf/R_eq^2     rpr             R_eq:")
                print(f"  {sci(NS.surf / NS.R_eq**2, ndig=6)}     {sci(NS.R_pr, ndig=6)}     {sci(NS.R_eq, ndig=6)}")
                print("\n         L/L_Edd         T_eff (keV)     f_c             w")

            print(f"  {' ' if k<10 else ''}{k+1}     {sci(NS.lum)}      {sci(NS.Epsilon_eff)}      {sci(NS.fc)}      {sci(NS.w)}")
        end = time()
        print(f'{mid-start:.4f}, {end-mid:.4f}, {end-start:.4f}')