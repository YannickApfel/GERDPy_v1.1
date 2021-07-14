# -*- coding: utf-8 -*-
""" Thermal load generator

    Author: Yannick Apfel
"""

import numpy as np
from scipy.constants import pi


def load(x):
    """
    Synthetic load profile of Bernier et al. (2004).

    Returns load y (in watts) at time x (in hours).
    """
    A = 2000.0
    B = 2190.0
    C = 80.0
    D = 2.0
    E = 0.01
    F = 0.0
    G = 0.95

    func = (168.0-C)/168.0
    for i in [1, 2, 3]:
        func += 1.0/(i*pi) * (np.cos(C*pi*i/84.0)-1.0) \
                          * (np.sin(pi*i/84.0*(x-B)))
    func = func*A*np.sin(pi/12.0*(x-B))  \
           *np.sin(pi/4380.0*(x-B))

    y = func + (-1.0)**np.floor(D/8760.0*(x-B))*abs(func) \
      + E*(-1.0)**np.floor(D/8760.0*(x-B))/np.sign(np.cos(D*pi/4380.0*(x-F))+G)
    return abs(-np.array([y])) * 3