from .math import *

def W_null(theta, args):
    return 0.0

def W_const(theta, args):
    return 0.75

def W_line(theta, args):
    theta_max = args
    th_max = theta_max * RAD
    W = (1 - theta/th_max) * hs(1 - theta/th_max)
    return W
