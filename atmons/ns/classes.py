"""
Description of model classes
"""

import numpy as np
import warnings
import toml
import os

from .src.metric import * 
from .src.radiacional import *     
from .src.rotational import *       
from .src.model import * 


class ConfigNeurtonStar:
    """
    Main class of the Neurton Star's Config
    """
    def __init__(self, name='base', config=None, **kwargs):
        # descript parameters
        self.name = name
        self.description = "The base config for model test."

        # base parameters
        self.r_ns = 12.0
        self.m_ns = 1.5
        self.v_rot = 700.0
        self.i_ang = 60.0
        self.th_star = 0.0
        self.w_par = [0.0]

        # key parameters
        self.spec_key = 'wfc'
        self.fc_key = 1
        self.flux_key = 1
        self.chem = 's1'
        self.rel = False
        self.w_func = 'const'

        # grid parameters
        self.n_phi = 120
        self.n_theta = 60
        self.n_nu = 500
        self.rng_phi = [0, 360]
        self.rng_theta = [0, 180]
        self.rng_erg = [0.1, 50.0]
        self.unnull = True
        self.zsch_key = False

        if config is not None:
            self.__set_config(config)

        if kwargs is not None:
            self.__set_config(kwargs)   

    def __set_config(self, config):
        if config is not None:
            for key, value in config.items():
                if key in self.__dict__:
                    self.__dict__[key] = value
                else:
                    warnings.warn(f"key '{key}' is not in config", UserWarning)
            self.__update_config()

    def __get_config(self, name):
        with open(f'configs/{name}.toml', 'r') as f:
            config = toml.load(f)
    
        for part, content in config.items():
            if part in self.config:
                for key, value in config[part].items():
                    if key in self.__dict__:
                        self.__dict__[key] = value
                    else:
                        warnings.warn(f"key '{key}' is not in config", UserWarning)
            else:
                warnings.warn(f"part '{part}' is not in config", UserWarning)
        self.__update_config()
        
    def __update_config(self):
        self.descript = {
            'name': self.name,
            'description': self.description,
        }

        self.base = {
            'r_ns': self.r_ns,
            'm_ns': self.m_ns,
            'v_rot': self.v_rot,
            'i_ang': self.i_ang,
            'th_star': self.th_star,
            'w_par': self.w_par,
        }

        self.keys = {
            'spec_key': self.spec_key,
            'fc_key': self.fc_key,
            'flux_key': self.flux_key,
            'chem': self.chem,
            'rel': self.rel,
            'w_func': self.w_func,
        }

        self.grid = {
            'n_phi': self.n_phi,
            'n_theta': self.n_theta,
            'n_nu': self.n_nu,
            'rng_phi': self.rng_phi,
            'rng_theta': self.rng_theta,
            'rng_erg': self.rng_erg,
            'unnull': self.unnull,
            'zsch_key': self.zsch_key,
        }

        self.config = {
            'descript': self.descript,
            'base': self.base,
            'keys': self.keys,
            'grid': self.grid,
        }
        if not os.path.isdir('configs'):
            os.mkdir('configs')
        with open(f'configs/{self.name}.toml', 'w') as f:
            toml.dump(self.config, f)

    def __call__(self, name=None, config=None, **kwargs):
        if name is not None:
            self.__get_config(name)
        if config is not None:
            self.__set_config(config)
        if kwargs is not None:
            self.__set_config(kwargs)

        return self

    def __str__(self):
        config = ""
        for part in self.config.keys():
            config += f"{part}:\n"
            for key, value in self.config[part].items():
                config += f"\t{key}: {value}\n"
        return config

    def __repr__(self):
        return self.__str__()
    
    def output(self):
        out = {}
        for part in self.config.keys():
            if part in self.config:
                for key, value in self.config[part].items():
                    if key in self.__dict__:
                        out[key] = value
        return out
    
class ParametersNeurtonStar:
    """
    Main class of the Neurton Star's Parameters
    """
    def __init__(self, name=None, config=None, save=True):
        # its config
        if config is None:
            config = ConfigNeurtonStar()
            if name is not None:
                config = config(name=name)

        self.config = config
        self.name = self.config.name,
        self.description = self.config.description,

        # model parameters
        self.th_star = config.th_star
        self.th_star <<= DEG
        # self.w_func, self.w_par = config.w_func, config.w_par

        # phisical parameters
        if config.rel:
            self.r = inverse(r_eq, (config.r_ns * KM), ((config.m_ns * M_SUN), (config.v_rot * u.Hz)), base=1.0, pw=5)
            self.m = inverse(m_cor, (config.m_ns * M_SUN), (config.r_ns * KM, config.v_rot * u.Hz), base=1.0, pw=5)
            self.r_eq, self.m_cor = config.r_ns, config.m_ns
        else:
            self.r, self.m = config.r_ns, config.m_ns
            self.r_eq, self.m_cor = 0.0, 0.0

        self.r <<= KM
        self.m <<= M_SUN 
        self.r_eq <<= KM
        self.m_cor <<= M_SUN
        self.R = R_NS(self.r, self.m)
        self.M = M_NS(self.m)
        self.R_sch = R_sch(self.M) # TODO: ВОПРОС ПРО НЕСКОРРЕКТИРОВАННЫЙ РАДИУС ШВРАЦШИЛЬДА
        self.zsch = zsch(self.R, self.R_sch)
        self.g = g(self.R, self.M, self.zsch)
        self.log_g = log(self.g / self.g.unit)
        self.area_0 = Surf(self.R, self.zsch)

        # chemical parameters
        self.X_hyd = X_hyd(config.chem)
        self.kappa_e = kappa_e(self.X_hyd)
        # self.w_b = w_b(config.chem, config.fc_key)
        # self.T_c = T_c(config.chem, config.fc_key)

        # photometrical parameters
        self.Flux_edd = Flux_edd(self.g, self.kappa_e)
        self.T_edd = T_obs(self.Flux_edd, self.zsch)
        self.Theta_edd = Theta(self.T_edd)
        self.Epsilon_edd = Epsilon(self.T_edd)
        self.Lum_edd = Lumen(self.Flux_edd, self.R)
        self.Lum_obs = Lumen_obs(self.Flux_edd, self.R, self.zsch)

        # rotatinal parameters
        self.nu_rot = config.v_rot * u.Hz
        self.incl_ang = (config.i_ang * u.deg).to(u.rad)
        self.sin_i = sin(self.incl_ang)
        self.cos_i = cos(self.incl_ang)
        self.omega_rot = omega(self.nu_rot)   

        # relativical parameters
        self.v_cr = nu_crit(self.r, self.m)
        self.v_rel = nu_relative(self.nu_rot, self.v_cr)
        self.m_cor = m_cor(self.m, self.r, self.nu_rot, define=config.rel, m=self.m_cor)
        self.M_cor = M_NS(self.m_cor)
        self.r_eq = r_eq(self.r, self.m, self.nu_rot, define=config.rel, r=self.r_eq)
        self.R_eq = R_NS(self.r_eq, self.m_cor)
        self.R_sch_cor = R_sch(self.M_cor)

        # metrical parameters
        self.chi = chi_metric(self.R_eq, self.M_cor)
        self.Omega = Omega_metric(self.R_eq, self.M_cor, self.omega_rot)
        self.q_c = q_c_metric(self.chi, self.Omega)
        self.b_c = b_c_metric(self.chi, self.Omega)
        self.i_bar = i_bar_metric(self.chi)
        self.I = I_NS(self.i_bar, self.R_eq, self.M_cor)
        self.J = J_NS(self.I, self.omega_rot)
        self.g_0 = g_0_metric(self.R_eq, self.M_cor, self.chi)
        self.V_kep = np.sqrt(self.g_0 * self.R_eq)
        self.V_rot = self.omega_rot * self.R_eq
    
        if save:
            self.__update_param()

    def __set_param(self, param):
        if param is not None:
            for key, value in param.items():
                if key in self.__dict__:
                    self.__dict__[key] = value
                else:
                    warnings.warn(f"key '{key}' is not in param", UserWarning)
            self.__update_param()

    def __get_param(self, name):
        self.config = self.config(name=name)
        
        with open(f'calc/{self.config.name}_param.toml', 'r') as f:
            param = toml.load(f)

        for part, content in param.items():
            if part in self.param:
                for key, value in param[part].items():
                    if key in self.__dict__:
                        if isinstance(value, tuple):
                            self.__dict__[key] = float(value[0]) * u.Unit(value[1])
                        else:
                            self.__dict__[key] = value
                    else:
                        warnings.warn(f"key '{key}' is not in param", UserWarning)
            else:
                warnings.warn(f"part '{part}' is not in param", UserWarning)
        
        self.__update_param()

    def __update_param(self):
        self.descript = {
            'name': self.config.name,
            'description': self.config.description,
        }

        self.model = {
            'th_star': (self.th_star.value, str(self.th_star.unit)),
        }

        self.phisical = {
            'r': (self.r.value, str(self.r.unit)),
            'm': (self.m.value, str(self.m.unit)),
            'r_eq': (self.r_eq.value, str(self.r_eq.unit)),
            'm_cor': (self.m_cor.value, str(self.m_cor.unit)),
            'R': (self.R.value, str(self.R.unit)),
            'M': (self.M.value, str(self.M.unit)),
            'R_sch': (self.R_sch.value, str(self.R_sch.unit)),
            'zsch': (self.zsch.value, str(self.zsch.unit)),
            'g': (self.g.value, str(self.g.unit)),
            'log_g': (self.log_g.value, str(self.log_g.unit)),
            'area_0': (self.area_0.value, str(self.area_0.unit)),
        }

        self.chemical = {
            'X_hyd': (self.X_hyd.value, str(self.X_hyd.unit)),
            'kappa_e': (self.kappa_e.value, str(self.kappa_e.unit)),
        }


        self.photometrical = {
            'Flux_edd': (self.Flux_edd.value, str(self.Flux_edd.unit)),
            'T_edd': (self.T_edd.value, str(self.T_edd.unit)),
            'Theta_edd': (self.Theta_edd.value, str(self.Theta_edd.unit)),
            'Epsilon_edd': (self.Epsilon_edd.value, str(self.Epsilon_edd.unit)),
            'Lum_edd': (self.Lum_edd.value, str(self.Lum_edd.unit)),
            'Lum_obs': (self.Lum_obs.value, str(self.Lum_obs.unit)),
        }

        self.rotatinal = {
            'nu_rot': (self.nu_rot.value, str(self.nu_rot.unit)),
            'incl_ang': (self.incl_ang.value, str(self.incl_ang.unit)),
            'sin_i': (self.sin_i.value, str(self.sin_i.unit)),
            'cos_i': (self.cos_i.value, str(self.cos_i.unit)),
            'omega_rot': (self.omega_rot.value, str(self.omega_rot.unit)),
        }

        self.relativical = {
            'v_cr': (self.v_cr.value, str(self.v_cr.unit)),
            'v_rel': (self.v_rel.value, str(self.v_rel.unit)),
            'm_cor': (self.m_cor.value, str(self.m_cor.unit)),
            'M_cor': (self.M_cor.value, str(self.M_cor.unit)),
            'r_eq': (self.r_eq.value, str(self.r_eq.unit)),
            'R_eq': (self.R_eq.value, str(self.R_eq.unit)),
            'R_sch_cor': (self.R_sch_cor.value, str(self.R_sch_cor.unit)),
        }

        self.metrical = {
            'chi': (self.chi.value, str(self.chi.unit)),
            'Omega': (self.Omega.value, str(self.Omega.unit)),
            'q_c': (self.q_c.value, str(self.q_c.unit)),
            'b_c': (self.b_c.value, str(self.b_c.unit)),
            'i_bar': (self.i_bar.value, str(self.i_bar.unit)),
            'I': (self.I.value, str(self.I.unit)),
            'J': (self.J.value, str(self.J.unit)),
            'g_0': (self.g_0.value, str(self.g_0.unit)),
            'V_kep': (self.V_kep.value, str(self.V_kep.unit)),
            'V_rot': (self.V_rot.value, str(self.V_rot.unit)),
        }

        self.param = {
            'descript': self.descript,
            'model': self.model,
            'phisical': self.phisical,
            'chemical': self.chemical,
            'photometrical': self.photometrical,
            'rotatinal': self.rotatinal,
            'relativical': self.relativical,
            'metrical': self.metrical,
        }   

        if not os.path.isdir('calc'):
            os.mkdir('calc')

        with open(f'calc/{self.config.name}_param.toml', 'w') as f:
            toml.dump(self.param, f)

    def __call__(self, name=None, config=None, **kwargs):
        if name is not None:
            self.__get_param(name)

        if config is not None:
            self.__init__(config=config)

        if kwargs is not None:
            self.__set_param(kwargs)

        return self
    
    def __str__(self):
        param = ""
        for part in self.param.keys():
            param += f"{part}:\n"
            for key, value in self.param[part].items():
                if isinstance(value, tuple):
                    val = float(value[0]) * u.Unit(value[1])
                    param += f"\t{key}: {val:.3e}\n"
                else:
                    val = value
                    param += f"\t{key}: {val}\n"
        return param

    def __repr__(self):
        return self.__str__()
    
    def output(self):
        out = {}
        for part in self.param.keys():
            if part in self.param:
                for key, value in self.param[part].items():
                    if key in self.__dict__:
                        if isinstance(value, tuple):
                            out[key] = float(value[0]) * u.Unit(value[1])
                        else:
                            out[key] = value
        return out


class NeurtonStar:
    """
    Main class of the Neurton Star
    """ 
    def __init__(self, config):
        self.config = config
        self.param = ParametersNeurtonStar(config=config)

