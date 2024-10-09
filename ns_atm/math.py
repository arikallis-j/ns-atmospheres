import numpy as np

from .const import *

def R_const(r_0, phi, theta, star):
    return r_0

class Point:
    """
    Class of point
    par: phi, theta
    opt: r, r_pt[r_func, args]
    """
    def __init__(self, sphere, coord, phi, theta, r_0=1.0, r_func=R_const):
        self.th, self.ph = coord
        self.sphere = sphere
        self.phi = phi * RAD
        self.theta = theta * RAD
        self.r_func = r_func
        self.r_0, self.R_0 = r_0, r_0 * KM 

        self.r = self.r_func(r_0, self.phi, self.theta, self.sphere)
        self.R = self.r * KM

    def __str__(self):
        phi_0 = round(self.phi * DEG,1)
        theta_0 = round(self.theta * DEG,1)
        r_0 = round(self.r, 1)
        return str((phi_0, theta_0, r_0))

    def get_par(self, par):
        dict_coord = {
            "r": self.R,
            "phi": self.phi,
            "theta": self.theta,
        }
        return dict_coord[par]
 
class Surface:
    """
    Class of point
    par: PointClass
    opt: r, r_pt[r_func, args]
    """
    def __init__(self, sphere, PointClass, r=1.0, r_func=R_const):
        self.sphere = sphere
        self.Point = PointClass
        self.size = (0, 0)
        self.ranges = ((0.0, 0.0), (0.0, 0.0))
        self.Phi, self.Theta = [], []

        self.Grid = []
        self.r, self.r_func = r, r_func
        self.R = r * KM
        self.is_filled = False

    def __str__(self):
        grid_str = ""
        n_phi, n_theta = self.size
        for th in range(n_theta):
            for ph in range(n_phi):
                pnt = self.get_point(th, ph)
                grid_str += f"{pnt} "
            grid_str += "\n"
        return grid_str

    def add_point(self, phi, theta, coord = (0, 0), is_filling = False):
        if self.is_filled:
            self.clear()
        if is_filling:
            th, ph = coord
        else:
            th, ph = len(self.Grid), 0
        r, r_func = self.r, self.r_func
        self.Grid.append(self.Point(self.sphere, coord, phi, theta, r, r_func))

    def get_point(self, th, ph):
        if self.is_filled:
            n_phi, n_theta = self.size
            k = th*n_phi + ph
            pnt = self.Grid[k]
            return pnt
        else:
            raise IndexError("The Grid is not filled")

    def clear(self):
        self.Grid = []
        self.is_filled = False

    def fill(self, size=(10, 10), ranges=((0, 360), (0, 180))):
        self.clear()
        self.size = size
        
        n_phi, n_theta = size
        ph_min, ph_max = ranges[0]
        th_min, th_max = ranges[1]

        for th in range(n_theta):
            for ph in range(n_phi):
                phi = ph_min + ph/(n_phi-1) * (ph_max - ph_min)
                theta = th_min + th/(n_theta-1) * (th_max - th_min)
                self.add_point(phi, theta, coord=(th, ph), is_filling=True)
                if ph==0:
                    self.Theta.append(theta * RAD)
                if th==0:
                    self.Phi.append(phi * RAD)
        self.is_filled = True
        self.ranges = ((ph_min * RAD, ph_max * RAD), (th_min * RAD, th_max * RAD))
        
class Star:
    """
    Class of star
    par: name
    opt: r
    """
    def __init__(self, name, r=12):
        self.name = name
        self.r, self.R = r, r * KM
        self.R_func = R_const
        self.args_pnt = ()
        self.Surface = Surface(self, Point, self.r, r_func=self.R_func)

        self.par = {
            "name": self.name,
            "r": self.r,
        }
 
    def init_grid(self, n_phi=10, n_theta=10, rng_phi=(0.0, 360.0), rng_theta=(0.0, 180.0), unnull=True):
        if unnull:
            dphi =  (rng_phi[1] - rng_phi[0]) / (n_phi)
            dtheta = (rng_theta[1] - rng_theta[0]) / (n_theta)
            rng_phi = (rng_phi[0] + dphi/2, rng_phi[1] - dphi/2)
            rng_theta = (rng_theta[0] + dtheta/2, rng_theta[1] - dtheta/2)

        self.Surface.fill(size=(n_phi, n_theta), ranges=(rng_phi, rng_theta))
    
    def init_point(self, phi=0.0, theta=0.0):
        self.Surface.add_point(phi, theta)

    def __str__(self):
        sphere_str = "Star with:\n"
        for key in self.par:
            value = self.par[key]
            sphere_str += f"{key} = {value} \n"
        return sphere_str

def sqrt(x):
    return np.sqrt(x)

def exp(x):
    return np.exp(x)

def exp10(x):
    return np.exp(x * np.log(10))

def log(x):
    return np.log10(x)

def ln(x):
    return np.log(x)

def base(x):
    return 10**(log(abs(x))//1)

def rnd(x, ndig=N_DIG_ERR):
    return round(x/base(x),ndig)

def sci(x, ndig=5):
    x_base = int(log(x)//1)
    x_norm = round(x/base(x), ndigits=ndig)
    x_str = f"{x_norm}e{'+' + str(x_base) if x_base>=0 else str(x_base)}"
    return x_str

def sin(x):
    return np.sin(x)

def cos(x):
    return np.cos(x)

def gaussian(x, mu=1.0, sig=1.0):
    return 1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)

def norm_gauss(x, mu=1.0, sig=1.0):
    return gaussian(x, mu=mu, sig=sig) / gaussian(mu, mu=mu, sig=sig)

def cutting(Spectrum, Gauss):
    CutSpectrum = []
    N = len(Spectrum)
    for k in range(N):
        if Spectrum[k]>Gauss[k]:
            CutSpectrum.append(Gauss[k])
        else:
            CutSpectrum.append(Spectrum[k])
    return CutSpectrum

def P_0(cos_th):
    return 1.0

def P_2(cos_th):
    return (3.0 * cos_th**2 - 1.0) / 2.0

def P_4(cos_th):
    return ((35.0 * cos_th**2 - 30.0) * cos_th**2 + 3.0) / 8.0

def dP_2(cos_th, sin_th):
    return -3.0 * cos_th * sin_th

def dP_4(cos_th, sin_th):
    return 2.5 * (3.0 - 7.0 * cos_th**2) * cos_th * sin_th

def hs(x):
    return np.heaviside(x)

def inverse(f, y0, args=(), pw=5, base=0.0):
    eps = 10.0**(-pw)
    h = 1.0
    x0, x1, x2 = base, base, base
    while abs(f(x0,*args) - y0) > (eps):
        x2 = x1 + h
        if f(x2,*args)>y0:
            h /= 2.0
        elif f(x2,*args)<y0:
            x1 = x2
        x0 = (x1 + x2)/2
    return round(x0, pw-1)

def map1(X, F, M, Y):  
    X = np.array(X)
    Y = np.array(Y)
    Y = Y.reshape(np.size(Y))

    F = np.array(F)
    G = np.zeros(Y.shape)
    G_full = np.full(Y.shape, False)


    N = np.full(Y.shape, M) - 1
    l, ll, l1, l2 = np.full(Y.shape, 1), np.full(Y.shape, -1), np.full(Y.shape, 0), np.full(Y.shape, 0)
    
    a, b, c, d = np.full(Y.shape, 0.0), np.full(Y.shape, 0.0), np.full(Y.shape, 0.0), np.full(Y.shape, 0.0)
    a1, b1, c1 = np.full(Y.shape, 0.0), np.full(Y.shape, 0.0), np.full(Y.shape, 0.0)
    af, bf, cf = np.full(Y.shape, 0.0), np.full(Y.shape, 0.0), np.full(Y.shape, 0.0)

    flag_30 = np.full(Y.shape, False)

    cycle_flag = np.logical_not(Y < X[l])
    while cycle_flag.any():
        l = np.where(cycle_flag, l+1, l)
        flag_30 =  np.where(l > N, True, flag_30)
        cycle_flag = np.logical_not(Y < X[l])
        cycle_flag = np.where(l >= N, False, cycle_flag)

    flag_31 = np.logical_not(l==ll)
    flag_3031 = np.logical_and(flag_30, flag_31)
    # flag 30
    l_3031 = np.where(l < N, l, N) 
    l = np.where(flag_3031, l_3031, l)
    
    flag = flag_3031
    
    c_3031 = np.zeros(c.shape) 
    c = np.where(flag, c_3031, c)

    b_3031 = (F[l]-F[l-1])/(X[l]-X[l-1])

    b = np.where(flag, b_3031, b)

    a_3031 = F[l]-X[l]*b 
    a = np.where(flag, a_3031, a)

    ll_3031 = l
    ll = np.where(flag_3031, ll_3031, ll)

    # flag 50
    flag = flag_30

    G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
    G_full = np.where(flag, True, G_full)
    if G_full.all():
        return G

    flag_20 = np.logical_not(l==ll)
    flag_21 = l==1
    flag_2021 = np.logical_and(flag_20, flag_21)
    flag_31 = np.logical_not(l==ll)
    flag_202131 = np.logical_and(flag_2021, flag_31)

    # flag 30
    l_202131 = np.where(l < N, l, N) 
    l = np.where(flag_202131, l_202131, l) 

    flag = flag_202131

    c_202131 = np.zeros(c.shape)  
    c = np.where(flag, c_202131, c)

    b_202131 = (F[l]-F[l-1])/(X[l]-X[l-1])
    b = np.where(flag, b_202131, b)

    a_202131 = F[l]-X[l]*b 
    a = np.where(flag, a_202131, a)

    ll_202131 = l
    ll = np.where(flag_202131, ll_202131, ll)

    # flag 50l
    flag = flag_2021
    G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
    G_full = np.where(flag, True, G_full)
    if G_full.all():
        return G

    # # #  
    flag_m21 = np.logical_not(flag_21)
    flag_20m21 = np.logical_and(flag_20, flag_m21)
    flag_311 = np.logical_or((l > (ll + 1)), l==2)
    flag_20m2131 = np.logical_and(flag_20m21, flag_311)

    # flag 31
    l1_20m2131 = l - 1
    l1 = np.where(flag_20m2131, l1_20m2131, l1)

    l2_20m2131 = l - 2
    l2 = np.where(flag_20m2131, l2_20m2131, l2)

    flag = flag_20m2131
    d_20m2131 = (F[l1] - F[l2]) / (X[l1] - X[l2])
    d = np.where(flag, d_20m2131, d)

    c1_20m2131 = F[l]/((X[l]-X[l1])*(X[l]-X[l2])) + (F[l2]/(X[l]-X[l2]) - F[l1]/(X[l]-X[l1]))/(X[l1]-X[l2]) 
    c1 = np.where(flag, c1_20m2131, c1)

    b1_20m2131 = d - (X[l1] + X[l2]) * c1                                 
    b1 = np.where(flag, b1_20m2131, b1)

    a1_20m2131 = F[l2] - X[l2] * d + X[l1] * X[l2] * c1 
    a1 = np.where(flag, a1_20m2131, a1)

    # flag_42
    flag_42 = np.logical_not(l < N)
    flag_20m213142 = np.logical_and(flag_20m2131, flag_42)
    flag = flag_20m213142

    c = np.where(flag, c1, c)                                                      
    b = np.where(flag, b1, b)                                                         
    a = np.where(flag, a1, a)                                                          
    ll= np.where(flag_20m213142, l, ll)      

    # flag 50
    G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
    G_full = np.where(flag, True, G_full)
    if G_full.all():
        return G

    # flag 25
    flag = flag_20m2131

    d_20m2131 = (F[l]-F[l1])/(X[l]-X[l1])  
    d = np.where(flag, d_20m2131, d)
    cf_20m2131 = F[l+1]/((X[l+1]-X[l])*(X[l+1]-X[l1])) + (F[l1]/(X[l+1]-X[l1]) - F[l]/(X[l+1]-X[l]))/(X[l]-X[l1])      
    cf = np.where(flag, cf_20m2131, cf)
    
    bf_20m2131 = d - (X[l] + X[l1]) * cf                                   
    bf = np.where(flag, bf_20m2131, bf)
    
    af_20m2131 = F[l1]-X[l1]*d + X[l]*X[l1]*cf 
    af = np.where(flag, af_20m2131, af)

    # flag 26
    eps = 1e-15
    WT = np.zeros(cf.shape)
    WT = np.where((abs(cf)) == WT, WT, abs(cf) / (abs(cf) + abs(c1) + eps))

    a_20m2131 = af + WT*(a1-af)  
    a = np.where(flag, a_20m2131, a)

    b_20m2131 = bf + WT*(b1-bf) 
    b = np.where(flag, b_20m2131, b) 

    c_20m2131 = cf + WT*(c1-cf) 
    c = np.where(flag, c_20m2131, c) 

    ll = np.where(flag_20m2131, l, ll) 

    # flag 50
    G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
    G_full = np.where(flag, True, G_full)
    if G_full.all():
        return G

    flag_m31 = np.logical_not(flag_31)
    flag_20m21m31 = np.logical_and(flag_20m21, flag_m31)
    flag_42 = (l==N)
    flag_20m21m3142 = np.logical_and(flag_20m21m31, flag_42)
    flag = flag_20m21m3142
    
    c = np.where(flag, c1, c)                                                      
    b = np.where(flag, b1, b)                                                         
    a = np.where(flag, a1, a)                                                          
    ll = np.where(flag, l, ll)  

    # flag 50
    flag = flag_20m213142
    G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
    G_full = np.where(flag, True, G_full)
    if G_full.all():
        return G

    flag_m42 = np.logical_not(flag_42)
    flag_20m21m31m42 = np.logical_and(flag_20m21m31, flag_m42)

    l1_20m21m31m42 = l - 1
    l1 = np.where(flag_20m21m31m42, l1_20m21m31m42, l1)

    flag = flag_20m21m31m42

    c1 = np.where(flag, cf, c1)                                                      
    b1 = np.where(flag, bf, b1)                                                         
    a1 = np.where(flag, af, a1)   

    # flag 25
    L = np.max(l)
    d_20m21m31m42 = (F[l]-F[l1])/(X[l]-X[l1])  
    d = np.where(flag, d_20m21m31m42, d)

    cf_20m21m31m42 = F[l+1]/((X[l+1]-X[l])*(X[l+1]-X[l1])) + (F[l1]/(X[l+1]-X[l1]) - F[l]/(X[l+1]-X[l]))/(X[l]-X[l1])      
    cf = np.where(flag, cf_20m21m31m42, cf)
    
    bf_20m21m31m42 = d - (X[l] + X[l1]) * cf                                   
    bf = np.where(flag, bf_20m21m31m42, bf)
    
    af_20m21m31m42 = F[l1]-X[l1]*d + X[l]*X[l1]*cf 
    af = np.where(flag, af_20m21m31m42, af)

    # flag 26
    eps = 1e-15
    WT = np.zeros(cf.shape)
    WT = np.where((abs(cf)) == WT, WT, abs(cf) / (abs(cf) + abs(c1) + eps))

    a_20m21m31m42 = af + WT*(a1-af)  
    a = np.where(flag, a_20m21m31m42, a)

    b_20m21m31m42 = bf + WT*(b1-bf) 
    b = np.where(flag, b_20m21m31m42, b) 

    c_20m21m31m42 = cf + WT*(c1-cf) 
    c = np.where(flag, c_20m21m31m42, c) 

    ll = np.where(flag, l, ll) 

    # flag 50
    G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
    G_full = np.where(flag, True, G_full)
    if G_full.all():
        return G

    flag_m20 = np.logical_not(flag_20) 
    flag = flag_m20

    G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
    G_full = np.where(flag, True, G_full)
    if G_full.all():
        return G

    return G


# def map1(X, F, M, Y):  
#     F = np.array(F)
#     X0 = np.array(X)
#     X = np.full(F.T.shape, X0).T
#     U_shape = Y.shape
#     Y = Y.reshape(Y.shape[0]*Y.shape[1],)
#     Y1 = np.full((F.T.shape[0], *Y.shape), Y).T

#     print(X.shape)
#     print(F.shape)
#     print(Y.shape)

#     G = np.zeros(Y1.shape)
#     G_full = np.full(Y1.shape, False)
#     l, ll, l1, l2 = np.full(Y.shape, 1), np.full(Y.shape, -1), np.full(Y.shape, 0), np.full(Y.shape,0)
#     print(l.shape)
#     print(f"F[l] {F[l].shape}")
#     print(f"X[l] {X[l].shape}")
#     print(f"F[l] - X[l] {(F[l] - X[l]).shape}")
       
#     N = np.full(Y.shape, M) - 1

#     a, b, c, d = np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0)
#     a1, b1, c1 = np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0)
#     af, bf, cf = np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0)

#     flag_30 = np.full(Y.shape, False)

#     cycle_flag = np.logical_not(Y < X[l][:,0])
#     while cycle_flag.any():
#         l = np.where(cycle_flag, l+1, l)
#         flag_30 =  np.where(l > N, True, flag_30)
#         cycle_flag = np.logical_not(Y < X[l][:,0])
#         cycle_flag = np.where(l >= N, False, cycle_flag)

#     flag_31 = np.logical_not(l==ll)
#     flag_3031 = np.logical_and(flag_30, flag_31)
#     # flag 30
#     l_3031 = np.where(l < N, l, N) 
#     l = np.where(flag_3031, l_3031, l)
    
#     flag = np.full((F.T.shape[0], *flag_3031.shape), flag_3031).T
    
#     c_3031 = np.zeros(c.shape) 
#     c = np.where(flag, c_3031, c)

#     b_3031 = (F[l]-F[l-1])/(X[l]-X[l-1])

#     b = np.where(flag, b_3031, b)

#     a_3031 = F[l]-X[l]*b 
#     a = np.where(flag, a_3031, a)

#     ll_3031 = l
#     ll = np.where(flag_3031, ll_3031, ll)

#     # flag 50
#     flag = np.full((F.T.shape[0], *flag_30.shape), flag_30).T 

#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_20 = np.logical_not(l==ll)
#     flag_21 = l==1
#     flag_2021 = np.logical_and(flag_20, flag_21)
#     flag_31 = np.logical_not(l==ll)
#     flag_202131 = np.logical_and(flag_2021, flag_31)

#     # flag 30
#     l_202131 = np.where(l < N, l, N) 
#     l = np.where(flag_202131, l_202131, l) 

#     flag = np.full((F.T.shape[0], *flag_202131.shape), flag_202131).T 

#     c_202131 = np.zeros(c.shape)  
#     c = np.where(flag, c_202131, c)

#     b_202131 = (F[l]-F[l-1])/(X[l]-X[l-1])
#     b = np.where(flag, b_202131, b)

#     a_202131 = F[l]-X[l]*b 
#     a = np.where(flag, a_202131, a)

#     ll_202131 = l
#     ll = np.where(flag_202131, ll_202131, ll)

#     # flag 50
#     flag = np.full((F.T.shape[0], *flag_2021.shape), flag_2021).T 
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     # # #  
#     flag_m21 = np.logical_not(flag_21)
#     flag_20m21 = np.logical_and(flag_20, flag_m21)
#     flag_311 = np.logical_or((l > (ll + 1)), l==2)
#     flag_20m2131 = np.logical_and(flag_20m21, flag_311)

#     # flag 31
#     l1_20m2131 = l - 1
#     l1 = np.where(flag_20m2131, l1_20m2131, l1)

#     l2_20m2131 = l - 2
#     l2 = np.where(flag_20m2131, l2_20m2131, l2)

#     flag = np.full((F.T.shape[0], *flag_20m2131.shape), flag_20m2131).T 
#     d_20m2131 = (F[l1] - F[l2]) / (X[l1] - X[l2])
#     d = np.where(flag, d_20m2131, d)

#     c1_20m2131 = F[l]/((X[l]-X[l1])*(X[l]-X[l2])) + (F[l2]/(X[l]-X[l2]) - F[l1]/(X[l]-X[l1]))/(X[l1]-X[l2]) 
#     c1 = np.where(flag, c1_20m2131, c1)

#     b1_20m2131 = d - (X[l1] + X[l2]) * c1                                 
#     b1 = np.where(flag, b1_20m2131, b1)

#     a1_20m2131 = F[l2] - X[l2] * d + X[l1] * X[l2] * c1 
#     a1 = np.where(flag, a1_20m2131, a1)

#     # flag_42
#     flag_42 = np.logical_not(l < N)
#     flag_20m213142 = np.logical_and(flag_20m2131, flag_42)
#     flag = np.full((F.T.shape[0], *flag_20m213142.shape), flag_20m213142).T 

#     c = np.where(flag, c1, c)                                                      
#     b = np.where(flag, b1, b)                                                         
#     a = np.where(flag, a1, a)                                                          
#     ll= np.where(flag_20m213142, l, ll)      

#     # flag 50
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     # flag 25
#     flag = np.full((F.T.shape[0], *flag_20m2131.shape), flag_20m2131).T 

#     d_20m2131 = (F[l]-F[l1])/(X[l]-X[l1])  
#     d = np.where(flag, d_20m2131, d)
#     cf_20m2131 = F[l+1]/((X[l+1]-X[l])*(X[l+1]-X[l1])) + (F[l1]/(X[l+1]-X[l1]) - F[l]/(X[l+1]-X[l]))/(X[l]-X[l1])      
#     cf = np.where(flag, cf_20m2131, cf)
    
#     bf_20m2131 = d - (X[l] + X[l1]) * cf                                   
#     bf = np.where(flag, bf_20m2131, bf)
    
#     af_20m2131 = F[l1]-X[l1]*d + X[l]*X[l1]*cf 
#     af = np.where(flag, af_20m2131, af)

#     # flag 26
#     eps = 1e-15
#     WT = np.zeros(cf.shape)
#     WT = np.where((abs(cf)) == WT, WT, abs(cf) / (abs(cf) + abs(c1) + eps))

#     a_20m2131 = af + WT*(a1-af)  
#     a = np.where(flag, a_20m2131, a)

#     b_20m2131 = bf + WT*(b1-bf) 
#     b = np.where(flag, b_20m2131, b) 

#     c_20m2131 = cf + WT*(c1-cf) 
#     c = np.where(flag, c_20m2131, c) 

#     ll = np.where(flag_20m2131, l, ll) 

#     # flag 50
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_m31 = np.logical_not(flag_31)
#     flag_20m21m31 = np.logical_and(flag_20m21, flag_m31)
#     flag_42 = (l==N)
#     flag_20m21m3142 = np.logical_and(flag_20m21m31, flag_42)
#     flag = np.full((F.T.shape[0], *flag_20m21m3142.shape), flag_20m21m3142).T 
    
#     c = np.where(flag, c1, c)                                                      
#     b = np.where(flag, b1, b)                                                         
#     a = np.where(flag, a1, a)                                                          
#     ll = np.where(flag, l, ll)  

#     # flag 50
#     flag = np.full((F.T.shape[0], *flag_20m213142.shape), flag_20m213142).T 
#     G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_m42 = np.logical_not(flag_42)
#     flag_20m21m31m42 = np.logical_and(flag_20m21m31, flag_m42)

#     l1_20m21m31m42 = l - 1
#     l1 = np.where(flag_20m21m31m42, l1_20m21m31m42, l1)

#     flag = np.full((F.T.shape[0], *flag_20m21m31m42.shape), flag_20m21m31m42).T

#     c1 = np.where(flag, cf, c1)                                                      
#     b1 = np.where(flag, bf, b1)                                                         
#     a1 = np.where(flag, af, a1)   

#     # flag 25
#     L = np.max(l)
#     d_20m21m31m42 = (F[l]-F[l1])/(X[l]-X[l1])  
#     d = np.where(flag, d_20m21m31m42, d)

#     cf_20m21m31m42 = F[l+1]/((X[l+1]-X[l])*(X[l+1]-X[l1])) + (F[l1]/(X[l+1]-X[l1]) - F[l]/(X[l+1]-X[l]))/(X[l]-X[l1])      
#     cf = np.where(flag, cf_20m21m31m42, cf)
    
#     bf_20m21m31m42 = d - (X[l] + X[l1]) * cf                                   
#     bf = np.where(flag, bf_20m21m31m42, bf)
    
#     af_20m21m31m42 = F[l1]-X[l1]*d + X[l]*X[l1]*cf 
#     af = np.where(flag, af_20m21m31m42, af)

#     # flag 26
#     eps = 1e-15
#     WT = np.zeros(cf.shape)
#     WT = np.where((abs(cf)) == WT, WT, abs(cf) / (abs(cf) + abs(c1) + eps))

#     a_20m21m31m42 = af + WT*(a1-af)  
#     a = np.where(flag, a_20m21m31m42, a)

#     b_20m21m31m42 = bf + WT*(b1-bf) 
#     b = np.where(flag, b_20m21m31m42, b) 

#     c_20m21m31m42 = cf + WT*(c1-cf) 
#     c = np.where(flag, c_20m21m31m42, c) 

#     ll = np.where(flag, l, ll) 

#     # flag 50
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)  
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_m20 = np.logical_not(flag_20) 
#     flag = np.full((F.T.shape[0], *flag_m20.shape), flag_m20).T

#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     return G.reshape(*U_shape, G.shape[1]) 


# def map2(X, F, M, Y):  
#     F = np.array(F)
#     X0 = np.array(X)
#     X = np.full(F.T.shape, X0).T
#     U_shape = Y.shape
#     Y = Y.reshape(Y.shape[0]*Y.shape[1],)
#     Y1 = np.full((F.T.shape[0], *Y.shape), Y).T

#     print(X.shape)
#     print(F.shape)
#     print(Y.shape)
#     print(Y1.shape)

#     G = np.zeros(Y1.shape)
#     G_full = np.full(Y1.shape, False)
#     l, ll, l1, l2 = np.full(Y.shape, 1), np.full(Y.shape, -1), np.full(Y.shape, 0), np.full(Y.shape,0)
#     print(l.shape)
#     print(f"F[l] {F[l].shape}")
#     print(f"X[l] {X[l].shape}")
#     print(f"F[l] - X[l] {(F[l] - X[l]).shape}")
       
#     N = np.full(Y.shape, M) - 1

#     a, b, c, d = np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0)
#     a1, b1, c1 = np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0)
#     af, bf, cf = np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0), np.full(Y1.shape, 0.0)

#     flag_30 = np.full(Y.shape, False)

#     cycle_flag = np.logical_not(Y < X[l][:,0])
#     while cycle_flag.any():
#         l = np.where(cycle_flag, l+1, l)
#         flag_30 =  np.where(l > N, True, flag_30)
#         cycle_flag = np.logical_not(Y < X[l][:,0])
#         cycle_flag = np.where(l >= N, False, cycle_flag)

#     flag_31 = np.logical_not(l==ll)
#     flag_3031 = np.logical_and(flag_30, flag_31)
#     # flag 30
#     l_3031 = np.where(l < N, l, N) 
#     l = np.where(flag_3031, l_3031, l)
    
#     flag = np.full((F.T.shape[0], *flag_3031.shape), flag_3031).T
    
#     c_3031 = np.zeros(c.shape) 
#     c = np.where(flag, c_3031, c)

#     b_3031 = (F[l]-F[l-1])/(X[l]-X[l-1])

#     b = np.where(flag, b_3031, b)

#     a_3031 = F[l]-X[l]*b 
#     a = np.where(flag, a_3031, a)

#     ll_3031 = l
#     ll = np.where(flag_3031, ll_3031, ll)

#     # flag 50
#     flag = np.full((F.T.shape[0], *flag_30.shape), flag_30).T 

#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_20 = np.logical_not(l==ll)
#     flag_21 = l==1
#     flag_2021 = np.logical_and(flag_20, flag_21)
#     flag_31 = np.logical_not(l==ll)
#     flag_202131 = np.logical_and(flag_2021, flag_31)

#     # flag 30
#     l_202131 = np.where(l < N, l, N) 
#     l = np.where(flag_202131, l_202131, l) 

#     flag = np.full((F.T.shape[0], *flag_202131.shape), flag_202131).T 

#     c_202131 = np.zeros(c.shape)  
#     c = np.where(flag, c_202131, c)

#     b_202131 = (F[l]-F[l-1])/(X[l]-X[l-1])
#     b = np.where(flag, b_202131, b)

#     a_202131 = F[l]-X[l]*b 
#     a = np.where(flag, a_202131, a)

#     ll_202131 = l
#     ll = np.where(flag_202131, ll_202131, ll)

#     # flag 50
#     flag = np.full((F.T.shape[0], *flag_2021.shape), flag_2021).T 
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     # # #  
#     flag_m21 = np.logical_not(flag_21)
#     flag_20m21 = np.logical_and(flag_20, flag_m21)
#     flag_311 = np.logical_or((l > (ll + 1)), l==2)
#     flag_20m2131 = np.logical_and(flag_20m21, flag_311)

#     # flag 31
#     l1_20m2131 = l - 1
#     l1 = np.where(flag_20m2131, l1_20m2131, l1)

#     l2_20m2131 = l - 2
#     l2 = np.where(flag_20m2131, l2_20m2131, l2)

#     flag = np.full((F.T.shape[0], *flag_20m2131.shape), flag_20m2131).T 
#     d_20m2131 = (F[l1] - F[l2]) / (X[l1] - X[l2])
#     d = np.where(flag, d_20m2131, d)

#     c1_20m2131 = F[l]/((X[l]-X[l1])*(X[l]-X[l2])) + (F[l2]/(X[l]-X[l2]) - F[l1]/(X[l]-X[l1]))/(X[l1]-X[l2]) 
#     c1 = np.where(flag, c1_20m2131, c1)

#     b1_20m2131 = d - (X[l1] + X[l2]) * c1                                 
#     b1 = np.where(flag, b1_20m2131, b1)

#     a1_20m2131 = F[l2] - X[l2] * d + X[l1] * X[l2] * c1 
#     a1 = np.where(flag, a1_20m2131, a1)

#     # flag_42
#     flag_42 = np.logical_not(l < N)
#     flag_20m213142 = np.logical_and(flag_20m2131, flag_42)
#     flag = np.full((F.T.shape[0], *flag_20m213142.shape), flag_20m213142).T 

#     c = np.where(flag, c1, c)                                                      
#     b = np.where(flag, b1, b)                                                         
#     a = np.where(flag, a1, a)                                                          
#     ll= np.where(flag_20m213142, l, ll)      

#     # flag 50
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     # flag 25
#     flag = np.full((F.T.shape[0], *flag_20m2131.shape), flag_20m2131).T 

#     d_20m2131 = (F[l]-F[l1])/(X[l]-X[l1])  
#     d = np.where(flag, d_20m2131, d)
#     cf_20m2131 = F[l+1]/((X[l+1]-X[l])*(X[l+1]-X[l1])) + (F[l1]/(X[l+1]-X[l1]) - F[l]/(X[l+1]-X[l]))/(X[l]-X[l1])      
#     cf = np.where(flag, cf_20m2131, cf)
    
#     bf_20m2131 = d - (X[l] + X[l1]) * cf                                   
#     bf = np.where(flag, bf_20m2131, bf)
    
#     af_20m2131 = F[l1]-X[l1]*d + X[l]*X[l1]*cf 
#     af = np.where(flag, af_20m2131, af)

#     # flag 26
#     eps = 1e-15
#     WT = np.zeros(cf.shape)
#     WT = np.where((abs(cf)) == WT, WT, abs(cf) / (abs(cf) + abs(c1) + eps))

#     a_20m2131 = af + WT*(a1-af)  
#     a = np.where(flag, a_20m2131, a)

#     b_20m2131 = bf + WT*(b1-bf) 
#     b = np.where(flag, b_20m2131, b) 

#     c_20m2131 = cf + WT*(c1-cf) 
#     c = np.where(flag, c_20m2131, c) 

#     ll = np.where(flag_20m2131, l, ll) 

#     # flag 50
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_m31 = np.logical_not(flag_31)
#     flag_20m21m31 = np.logical_and(flag_20m21, flag_m31)
#     flag_42 = (l==N)
#     flag_20m21m3142 = np.logical_and(flag_20m21m31, flag_42)
#     flag = np.full((F.T.shape[0], *flag_20m21m3142.shape), flag_20m21m3142).T 
    
#     c = np.where(flag, c1, c)                                                      
#     b = np.where(flag, b1, b)                                                         
#     a = np.where(flag, a1, a)                                                          
#     ll = np.where(flag, l, ll)  

#     # flag 50
#     flag = np.full((F.T.shape[0], *flag_20m213142.shape), flag_20m213142).T 
#     G = np.where(np.logical_not(G_full),  a + b*Y + c * Y**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_m42 = np.logical_not(flag_42)
#     flag_20m21m31m42 = np.logical_and(flag_20m21m31, flag_m42)

#     l1_20m21m31m42 = l - 1
#     l1 = np.where(flag_20m21m31m42, l1_20m21m31m42, l1)

#     flag = np.full((F.T.shape[0], *flag_20m21m31m42.shape), flag_20m21m31m42).T

#     c1 = np.where(flag, cf, c1)                                                      
#     b1 = np.where(flag, bf, b1)                                                         
#     a1 = np.where(flag, af, a1)   

#     # flag 25
#     L = np.max(l)
#     d_20m21m31m42 = (F[l]-F[l1])/(X[l]-X[l1])  
#     d = np.where(flag, d_20m21m31m42, d)

#     cf_20m21m31m42 = F[l+1]/((X[l+1]-X[l])*(X[l+1]-X[l1])) + (F[l1]/(X[l+1]-X[l1]) - F[l]/(X[l+1]-X[l]))/(X[l]-X[l1])      
#     cf = np.where(flag, cf_20m21m31m42, cf)
    
#     bf_20m21m31m42 = d - (X[l] + X[l1]) * cf                                   
#     bf = np.where(flag, bf_20m21m31m42, bf)
    
#     af_20m21m31m42 = F[l1]-X[l1]*d + X[l]*X[l1]*cf 
#     af = np.where(flag, af_20m21m31m42, af)

#     # flag 26
#     eps = 1e-15
#     WT = np.zeros(cf.shape)
#     WT = np.where((abs(cf)) == WT, WT, abs(cf) / (abs(cf) + abs(c1) + eps))

#     a_20m21m31m42 = af + WT*(a1-af)  
#     a = np.where(flag, a_20m21m31m42, a)

#     b_20m21m31m42 = bf + WT*(b1-bf) 
#     b = np.where(flag, b_20m21m31m42, b) 

#     c_20m21m31m42 = cf + WT*(c1-cf) 
#     c = np.where(flag, c_20m21m31m42, c) 

#     ll = np.where(flag, l, ll) 

#     # flag 50
#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)  
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     flag_m20 = np.logical_not(flag_20) 
#     flag = np.full((F.T.shape[0], *flag_m20.shape), flag_m20).T

#     G = np.where(np.logical_not(G_full),  a + b*Y1 + c * Y1**2, G)
#     G_full = np.where(flag, True, G_full)
#     if G_full.all():
#         return G.reshape(*U_shape, G.shape[1])

#     return G.reshape(*U_shape, G.shape[1]) 


# def f_30(X, F, N, l):
#     l = min(N,l) - 1                                                  
#     c = 0.0                                                             
#     b = (F[l]-F[l-1])/(X[l]-X[l-1])                        
#     a = F[l]-X[l]*b                                            
#     ll = l
#     return a, b, c, (l+1), (ll+1)

# def f_21(X, F, N, l, l1):
#     l, l1 = l - 1, l1 - 1
#     l2 = l - 2                                                          
#     d = (F[l1] - F[l2]) / (X[l1] - X[l2])
#     c1 = F[l]/((X[l]-X[l1])*(X[l]-X[l2])) + (F[l2]/(X[l]-X[l2]) - F[l1]/(X[l]-X[l1]))/(X[l1]-X[l2]) 
#     b1 = d - (X[l1] + X[l2]) * c1                                 
#     a1 = F[l2] - X[l2] * d + X[l1] * X[l2] * c1  
#     return a1, b1, c1, l2+1

# def f_25(X, F, N, l, l1):
#     l, l1 = l - 1, l1 - 1
#     d = (F[l]-F[l1])/(X[l]-X[l1])                          
#     cf = F[l+1]/((X[l+1]-X[l])*(X[l+1]-X[l1])) + (F[l1]/(X[l+1]-X[l1]) - F[l]/(X[l+1]-X[l]))/(X[l]-X[l1])   
#     bf = d - (X[l] + X[l1]) * cf                                   
#     af = F[l1]-X[l1]*d + X[l]*X[l1]*cf 
#     return af, bf, cf

# def f_26(a1, b1, c1, af, bf, cf, l):
#     WT = 0.0
#     if (abs(cf)!=0.0):
#         WT = abs(cf) / (abs(cf)+abs(c1))
#     a = af + WT*(a1-af)                                            
#     b = bf + WT*(b1-bf)                                            
#     c = cf + WT*(c1-cf)                                            
#     ll = l  
#     return a, b, c, ll

# def f_22(a, b, c, l, ll ):
#     return a, b, c, l, ll 

# def map1(X, F, N, Y):
#     G = 0.0
#     l, ll, l1 = 2, 0, 0

#     a, b, c = 1.0, 1.0, 1.0
#     a1, b1, c1 = 1.0, 1.0, 1.0
#     af, bf, cf = 1.0, 1.0, 1.0

#     y = Y #Y[1]

#     while not(y < X[l-1]):
#         l += 1
#         if (l > N): 
#             break
#     if (l > N) or (l == 2):
#         a, b, c, l, ll = f_30(X, F, N, l)
#     else:
#         l1 = l - 1
#         a1, b1, c1, l2 = f_21(X, F, N, l, l1)
#         if (l != N): 
#             af, bf, cf = f_25(X, F, N, l, l1)    
#             a, b, c, ll = f_26(a1, b1, c1, af, bf, cf, l)       
#         else:
#             a, b, c, ll = a1, b1, c1, l
#     G = a + b*y + c*y**2
#     return G 

# def map1(X, F, N, Y):
#     l, ll, l1, l2 = 2, 0, 0, 0
#     a, b, c = 1.0, 1.0, 1.0
#     a1, b1, c1 = 1.0, 1.0, 1.0
#     af, bf, cf = 1.0, 1.0, 1.0
#     for k in range(1):
#         y = Y
#         while not(y < X[l]):
#             l += 1
#             if (l > N):
#                 if (l != ll):
#                     a, b, c, l, ll = f_30(X, F, N, l)
#                 break

#         if (l != ll) and (l != 2):
#             l1 = l - 1 
#             if (l > (ll + 1) or l == 3):
#                 a1, b1, c1, l2 = f_21(X, F, N, l, l1)

#                 if (l >= N):    
#                     a, b, c, l, ll = a1, b1, c1, l, l    
#                 else: 
#                     af, bf, cf = f_25(X, F, N, l, l1) 
#                     a, b, c, ll = f_26(a1, b1, c1, af, bf, cf, l)                                               
#             else: 
#                 c1 = cf                                                       
#                 b1 = bf                                                       
#                 a1 = af

#                 if (l==N):
#                     a, b, c, ll = f_22(a1, b1, c1, l)
#                 else:
#                     af, bf, cf = f_25(X, F, N, l, l1)    
#                     a, b, c, ll = f_26(a1, b1, c1, af, bf, cf, l)   


#         elif (l != ll) and (l == 2):
#             a, b, c, l, ll = f_30(X, F, N, l)

#         G = a + b*y + c*y**2
#     return G 
