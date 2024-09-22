from time import time

from .classes import *
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
        Colors = ['black', 'orange', 'violet', 'magenta', 'blue', 'red']
        N, M, k = 5, 200, 0
        rel = True
        m_cor = 1.5
        r_eq = 14
        v_rot = 600
        n_phi = 2*N
        n_theta = N
        n_nu = M
        lum_k = FLUX_REL[k]
        FWHM = 1
        for k in range(1):
            NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
                            fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                            m_ns=m_cor, r_ns=r_eq, v_rot=v_rot, i_ang=I_arr[k], lum=lum_k)
            NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
            Energy = np.array(NS.E)
            Flux = np.array(NS.B_real) / max(NS.B_real)
            for i in range(len(Energy)):
                if Energy[i]>1.0:
                    norm = (Flux[i-1] + Flux[i])/2
                    break
            sig = 1/(FWHM * 2 * sqrt(2 * ln(2)))
            #sig = 1 / (norm * sqrt(2*PI))
            Gauss = gaussian(Energy, mu = 1.0, sig=sig)
            Cut_spectrum = cutting(Flux, Gauss)
            #Norm_flux = np.array(Cut_spectrum) / max(Cut_spectrum)
            Norm_flux = np.array(Gauss) / max(Gauss)
            factor = sqrt(1 - NS.R_sch/NS.R_eq)
            print(k, factor)
            #raw_graph(Energy / factor, Norm_flux, color=Colors[k])
            draw_graph(Energy, Norm_flux, color=Colors[k])
        
        show_graph()

    def spectra(self):
        I_arr = [0, 90]
        Colors = ['blue', 'red']
        xlabel = "Photon energy (keV)"
        ylabel = "$L_E (10^{36} erg^{-1} keV^{-1} sr^{-1})$"
        N, M = 10, 500
        rel = True
        m_cor = 1.519
        r_eq = 15.48
        v_rot = 700
        n_phi = 2*N
        n_theta = N
        n_nu = M
        lum_k = FLUX_REL[0]
        print(lum_k)
        for k in range(2):
            NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
                            fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                            m_ns=m_cor, r_ns=r_eq, v_rot=v_rot, i_ang=I_arr[k], lum=lum_k)
            NS.calc(n_phi=n_phi, n_theta=n_theta, n_nu=n_nu)
            Energy = np.array(NS.E)
            Flux = np.array(NS.B_real) / 10**36
            draw_graph(Energy, Flux, color=Colors[k], loglog=True, xrange=(1,20), yrange=(0.003,0.4), title="$L/L_{edd} = 0.1$", xlabel=xlabel, ylabel=ylabel)
        show_graph()

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
        
        for k in range(0, N_MODEL):
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
    
    def test(self):
        lum_k = FLUX_REL[0]
        NS = NeurtonStar(name=name, sys=system, chem=chem, rel=rel,
                        fc_key=fc_key, flux_key=flux_key, w_key=w_key, 
                        m_ns=m_ns, r_ns=r_ns, v_rot=v_rot, i_ang=i_ang, lum=lum_k)
        
        NS.calc(n_phi=130, n_theta=70, n_nu=n_nu)
        print(f"lum:{NS.lum} energy: {NS.Epsilon_eff}")
        print(f"w:{NS.w} fc: {NS.fc}")

        t1 = time()
        print("Testing")
        t2 = time()
        print(t2 - t1)
        return "Tests ended"
    