from time import time

from .classes import *

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
        # println("\n  surf/R_eq^2     rpr             R_eq:")
        # println(@sprintf(" %13.6e   %13.6e   %13.6e", surf / R_eq^2, rpr, R_eq))
        # println("\n         L/L_Edd         T_eff (keV)     f_c             w")

        # t1 = time()
        # print("Testing")
        # t2 = time()
        # print(t2 - t1)
        return "Tests ended"
    