import unittest
import numpy as np

import astropy.units as u
import astropy.constants as astro_const

from atmons.ns.const import *
from atmons.ns.phisical import *
from atmons.ns.metric import *

class PaperTest(unittest.TestCase):
    def test_one(self):
        self.assertEqual(1, 1)

    def test_line(self):
        for k in range(0, 10):
            with self.subTest(i=k):
                self.assertEqual(1, 1)

class ConstTest(unittest.TestCase):
    def test_math_const(self):
        self.assertEqual(np.pi, PI)
        self.assertEqual(np.e, EULER)

    def test_phisical_const(self):
        self.assertEqual(astro_const.c.cgs, C_LIGHT)
        self.assertEqual(astro_const.G.cgs, G_GRAV)
        self.assertEqual(astro_const.h.cgs, H_PL)
        self.assertEqual(astro_const.k_B.cgs, K_B)
        self.assertEqual(astro_const.sigma_sb.cgs, SIGMA_SB)

    def test_phisical_units(self):
        self.assertEqual(u.deg, DEG)
        self.assertEqual(u.rad, RAD)
        self.assertEqual(u.km, KM)
        self.assertEqual(u.cm, CM)
        self.assertEqual(u.M_sun, M_SUN)
        self.assertEqual(u.erg, ERG)
        self.assertEqual(u.keV, KEV)

class NS_Parameters_Test(unittest.TestCase):
    def test_phisical_param(self):
        self.assertEqual(r_min((1*u.Msun).cgs).unit, u.km)
        self.assertEqual(R_NS(10*u.km, 1*u.Msun).unit, u.cm)
        self.assertEqual(M_NS(1*u.Msun).unit, u.g)
        self.assertEqual(R_sch(1*u.Msun).unit, u.cm)
        self.assertEqual(zsch(1*u.km, 1*u.cm).unit, u.Unit())
        self.assertEqual(g(1*u.km, 1*u.Msun, 1*u.Unit()).unit, u.cm / u.s**2)
        self.assertEqual(log((1*u.cm/u.s**2) / (1*u.cm/u.s**2).unit).unit, u.Unit)
        self.assertEqual(Surf(1*u.km, 1*u.Unit()).unit, u.cm**2)

    def test_chemical_param(self):
        self.assertEqual(X_hyd("s1").unit, u.Unit())
        self.assertEqual(kappa_e(0*u.Unit()).unit, u.cm**2 / u.g)
        self.assertEqual(w_b("s1", 1).unit, u.Unit())
        self.assertEqual(T_c("s1", 1).unit, u.keV)
        
    def test_photometrical_param(self):
        self.assertEqual(Flux_edd(1*u.m/u.s**2, 1*u.cm**2/u.g).unit, u.erg/(u.s*u.cm**2))
        self.assertEqual(Flux_SB(1*u.K).unit, u.erg/(u.s*u.cm**2))
        self.assertEqual(T_SB(1*u.erg/(u.s*u.cm**2)).unit, u.K)
        self.assertEqual(T_obs(1*u.erg/(u.s*u.cm**2), 1*u.Unit()).unit, u.K)
        self.assertEqual(Theta(1*u.K).unit, u.erg)
        self.assertEqual(Epsilon(1*u.K).unit, u.keV)
        self.assertEqual(Lumen(1*u.erg/(u.s*u.cm**2), 1*u.km).unit, u.erg/u.s)
        self.assertEqual(Lumen_obs(1*u.erg/(u.s*u.cm**2), 1*u.km, 1*u.Unit()).unit, u.erg/u.s)
    
    def test_rotatinal_param(self):
        self.assertEqual(((30*u.deg).to(u.rad)).unit, u.rad)
        self.assertEqual((sin((30*u.deg).to(u.rad))).unit, u.Unit())
        self.assertEqual(omega(1*u.Hz).unit, u.s**(-1))

    def test_relativical_param(self):
        self.assertEqual(nu_crit(1*u.km, 1*u.Msun).unit, u.Hz)
        self.assertEqual(nu_relative(1*u.Hz, 1*u.Hz).unit, u.Unit())
        self.assertEqual(r_eq(1*u.km, 1*u.Msun, 1*u.Hz).unit, u.km)
        self.assertEqual(m_cor(1*u.Msun, 1*u.km, 1*u.Hz).unit, u.Msun)
    
    def test_metrical_param(self):
        self.assertEqual(chi_metric(1*u.km, 1*u.Msun).unit, u.Unit())
        self.assertEqual(Omega_metric(1*u.km, 1*u.Msun, 1*u.Hz).unit, u.Unit())
        self.assertEqual(q_c_metric(1*u.Unit(), 1*u.Unit()).unit, u.Unit())
        self.assertEqual(b_c_metric(1*u.Unit(), 1*u.Unit()).unit, u.Unit())
        self.assertEqual(i_bar_metric(1*u.Unit()).unit, u.Unit())
        self.assertEqual(I_NS(1*u.Unit(), 1*u.cm, 1*u.g).unit, u.cm**2 * u.g)
        self.assertEqual(J_NS(1*u.cm**2*u.g, 1/u.s).unit, u.cm**2 * u.g / u.s)
        self.assertEqual(g_0_metric(1*u.km, 1*u.Msun, 1*u.Unit()).unit, u.cm / u.s**2)
        
class NS_Surface_Test(unittest.TestCase):
    def test_one(self):
        self.assertEqual(1, 1)


if __name__ == "__main__":
    unittest.main()