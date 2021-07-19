# -*- coding: utf-8 -*-
""" Import und Aufbereitung der Wetterdaten

    Autor: Yannick Apfel
"""
import numpy as np


def get_u_inf(Nt):

    # wind speed vector of Nt random floats between 0 and 5 ([m/s])
    u_inf = np.random.uniform(0, 5, Nt)

    return u_inf


def get_T_inf(Nt):

    # temperature vector of Nt random floats between -2 and 5 (°C)
    T_inf = np.random.uniform(-2, 3, Nt)

    return T_inf
