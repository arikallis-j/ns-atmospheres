import numpy as np
from icecream import ic
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

#оригинальная функция map1 с переводом на Python
def map0(xold, fold, xnew):
    # Определение переменных
    l = 1
    N = len(xold)

    # Определяем точку относительно массива xold
    while l < N and xnew >= xold[l]:
        l += 1
        
    # Краевой случай для l == 1 or l == N
    if l == 1 or l == N:
        ic("Border case: line interpolation")
        l = min(N - 1, l)
        x1, x0 = xold[l], xold[l-1]
        f1, f0 = fold[l], fold[l-1]
        c_l = 0.0
        b_l = (f1 - f0) / (x1 - x0)
        a_l = f1 - x1 * b_l
        return a_l + (b_l + c_l * xnew) * xnew

    # Основная квадратичная интерполяция

    x2, x1, x0 = xold[l], xold[l-1], xold[l-2]
    f2, f1, f0 = fold[l], fold[l-1], fold[l-2]

    d  = (f1 - f0) / (x1 - x0)
    c = (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0)
    b = d - (x1 + x0) * c
    a = f0 - x0 * d + x1 * x0 * c

    cm , bm, am = c, b, a

    # Дополнительная квадратичная интерполяция, если не (почти) краевой случай
    if l != N - 1:
        ic("General case: quadric interpolation with correction")
        x2, x1, x0 = xold[l+1], xold[l], xold[l-1]
        f2, f1, f0 = fold[l+1], fold[l], fold[l-1]

        d  = (f1 - f0) / (x1 - x0)
        c = (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0)
        b = d - (x1 + x0) * c
        a = f0 - x0 * d + x1 * x0 * c

        cp , bp, ap = c, b, a

        wt = abs(cp) / (abs(cp) + abs(cm)) if abs(cp) != 0 else 0.0

        a = ap + wt * (am - ap)
        b = bp + wt * (bm - bp)
        c = cp + wt * (cm - cp)

    else:
        ic("Specific case: quadric interpolation without correction")

    return a + (b + c * xnew) * xnew

def map1(X, F, Y):     
    X = np.array(X)
    F = np.array(F)
    Y = np.array(Y)
    #print(X.shape, F.shape, Y.shape)

    l_shape = X.shape[0]
    p_shape = F.shape[1]
    g_shape_a = Y.shape[0]
    g_shape_b = Y.shape[1]

    X = np.full((g_shape_b, g_shape_a, p_shape, l_shape), X.T).T
    Y = np.full((p_shape, g_shape_a, g_shape_b), Y)

    X = X.reshape((l_shape, p_shape*g_shape_a*g_shape_b))
    F = F.reshape((l_shape, p_shape*g_shape_a*g_shape_b))
    
    Y = Y.reshape((p_shape*g_shape_a*g_shape_b))

    d_shape = F[0,:].shape

    R = np.full(d_shape, 0.0)

    N = l_shape - 1

    for j in range(X.shape[1]):
        #ic(100*j/X.shape[1], "%")
        a, b, c = 0.0, 0.0, 0.0
        l = 1
        #ic(Y[j], X[l,j], Y[j] >= X[l,j])
        while l < N and Y[j] >= X[l,j]:
            l += 1
            
        # Краевой случай для l == 1 or l == N
        if l == 1 or l == N:
            #ic("Border case: line interpolation")
            l = min(N - 1, l)
            x1, x0 = X[l,j], X[l-1,j]
            f1, f0 = F[l,j], F[l-1,j]
            c = 0.0
            b = (f1 - f0) / (x1 - x0)
            a = f1 - x1 * b
            #ic(1,l,a,b,c)
        
        # Основная квадратичная интерполяция
        else:
            x2, x1, x0 = X[l,j], X[l-1,j], X[l-2,j]
            f2, f1, f0 = F[l,j], F[l-1,j], F[l-2,j]

            d  = (f1 - f0) / (x1 - x0)
            c = (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0)
            b = d - (x1 + x0) * c
            a = f0 - x0 * d + x1 * x0 * c

            cm , bm, am = c, b, a
            
            # Дополнительная квадратичная интерполяция, если не (почти) краевой случай
            if l != N - 1:
                #ic("General case: quadric interpolation with correction")
                x2, x1, x0 = X[l+1,j], X[l,j], X[l-1,j]
                f2, f1, f0 = F[l+1,j], F[l,j], F[l-1,j]

                d  = (f1 - f0) / (x1 - x0)
                c = (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0)
                b = d - (x1 + x0) * c
                a = f0 - x0 * d + x1 * x0 * c

                cp , bp, ap = c, b, a

                wt = abs(cp) / (abs(cp) + abs(cm)) if abs(cp) != 0 else 0.0

                a = ap + wt * (am - ap)
                b = bp + wt * (bm - bp)
                c = cp + wt * (cm - cp)
                #ic(3,l,a,b,c)
            else:
                pass
                #ic(2,l,a,b,c)
                #ic("Specific case: quadric interpolation without correction")


        R[j] = a + (b + c * Y[j]) * Y[j]
        #ic(R[j])
    return R.reshape(p_shape, g_shape_a, g_shape_b)


def map2(X, F, Y):
    X = np.array(X)
    F = np.array(F)
    Y = np.array(Y)
    #print(X.shape, F.shape, Y.shape)

    l_shape = X.shape[0]
    p_shape = F.shape[1]
    g_shape_a = Y.shape[0]
    g_shape_b = Y.shape[1]

    X = np.full((g_shape_b, g_shape_a, p_shape, l_shape), X.T).T.reshape((l_shape, -1))
    Y = np.full((p_shape, g_shape_a, g_shape_b), Y)

    F = F.reshape((l_shape, -1))
    
    Y = Y.ravel()

    d_shape = F[0,:].shape

    R = np.zeros(d_shape)

    N = l_shape - 1

    a, b, c, d = np.zeros(d_shape), np.zeros(d_shape), np.zeros(d_shape), np.zeros(d_shape)

    cp, bp, ap = np.zeros(d_shape), np.zeros(d_shape), np.zeros(d_shape)
    cm, bm, am = np.zeros(d_shape), np.zeros(d_shape), np.zeros(d_shape)

    x2, x1, x0 = np.full(d_shape, 2.0), np.full(d_shape, 1.0), np.full(d_shape, 0.0)
    f2, f1, f0 = np.full(d_shape, 2.0), np.full(d_shape, 1.0), np.full(d_shape, 0.0)

    j = np.arange(d_shape[0])
    l = np.full(d_shape, 1)

    mask_cycle = (l < N) & (Y >= X[l, j])

    while mask_cycle.any():
        l = np.where(mask_cycle, l + 1, l)
        mask_cycle = (l < N) & (Y >= X[l, j])


    # Краевой случай для l == 1 or l == N
    mask_board = (l == 1) | (l == N)
    l = np.where(mask_board, np.minimum(N - 1, l), l)

    x1 = np.where(mask_board, X[l, j], x1)
    x0 = np.where(mask_board, X[l-1, j], x0)
    f1 = np.where(mask_board, F[l, j], f1)
    f0 = np.where(mask_board, F[l-1, j], f0)

    c = np.where(mask_board, 0.0, c)
    b = np.where(mask_board, (f1 - f0) / (x1 - x0), b)
    a = np.where(mask_board, f1 - x1 * b  , a)

    # Основная квадратичная интерполяция
    mask_main = ~mask_board
    
    x2 = np.where(mask_main, X[l,j], x2)
    x1 = np.where(mask_main, X[l-1,j], x1)
    x0 = np.where(mask_main, X[l-2,j], x0)

    f2 = np.where(mask_main, F[l,j], f2)
    f1 = np.where(mask_main, F[l-1,j], f1)
    f0 = np.where(mask_main, F[l-2,j], f0)

    d = np.where(mask_main, (f1 - f0) / (x1 - x0), d)
    c = np.where(mask_main, (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0), c)
    b = np.where(mask_main,  d - (x1 + x0) * c, b) 
    a = np.where(mask_main, f0 - x0 * d + x1 * x0 * c, a) 


    cm, bm, am = c, b, a

    # Дополнительная квадратичная интерполяция, если не (почти) краевой случай
    mask_center = mask_main & (l != N - 1)

    x2 = np.where(mask_center, X[l+1,j], x2)
    x1 = np.where(mask_center, X[l,j], x1)
    x0 = np.where(mask_center, X[l-1,j], x0)

    f2 = np.where(mask_center, F[l+1,j], f2)
    f1 = np.where(mask_center, F[l,j], f1)
    f0 = np.where(mask_center, F[l-1,j], f0)

    d = np.where(mask_center, (f1 - f0) / (x1 - x0), d)
    c = np.where(mask_center, (f2 / ((x2 - x1) * (x2 - x0))) + ((f0 / (x2 - x0)) - (f1 / (x2 - x1))) / (x1 - x0), c)
    b = np.where(mask_center,  d - (x1 + x0) * c, b) 
    a = np.where(mask_center, f0 - x0 * d + x1 * x0 * c, a) 

    cp, bp, ap = c, b, a
    determinator = np.where(np.abs(cp) != 0, np.abs(cp) + np.abs(cm), 1.0)
    wt = np.where(np.abs(cp) != 0, np.abs(cp) / determinator, 0.0)
    wt = np.where(mask_center, wt, 0.0)


    a = np.where(mask_center, ap + wt * (am - ap), a)
    b = np.where(mask_center, bp + wt * (bm - bp), b)
    c = np.where(mask_center, cp + wt * (cm - cp), c)

    R = a + (b + c * Y) * Y

    return R.reshape(p_shape, g_shape_a, g_shape_b)