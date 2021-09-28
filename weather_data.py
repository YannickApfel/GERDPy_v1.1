# -*- coding: utf-8 -*-
""" Autor: Yannick Apfel
"""
import numpy as np


def load_u_inf(Nt):

    # Windgeschwindigkeit (Vektor mit Nt random floats zwischen 0 und 5) [m/s]
    u_inf = np.random.uniform(0, 5, Nt)
    # u_inf = np.ones(Nt)

    return u_inf


def load_T_inf(Nt):

    # Umgebungstemperatur (Vektor mit Nt random floats zwischen -2 und 3) [°C]
    T_inf = np.random.uniform(-2, 3, Nt)
    # T_inf = np.ones(Nt) * -1
    T_inf[50:60] = 12

    return T_inf


def load_S_w(Nt):

    # Schneefallrate (Vektor mit Nt random floats zwischen 0 und 1.5) [mm/h]
    S_w = np.random.uniform(0, 1.5, Nt)
    # S_w = np.ones(Nt)

    return S_w


def load_B(Nt):
    
    # Bewölkungsgrad (Vektor mit Nt ints zwischen 0 und 8) [-]
    B = np.random.randint(9, size=Nt)
    # B = np.ones(Nt) * 3
    
    return B


def load_Phi(Nt):
    
    # rel. Luftfeuchte (Vektor mit Nt floats zwischen 0.63 und 0.96) [%]
    Phi = np.random.uniform(0.63, 0.96, Nt)
    # Phi = np.ones(Nt) * 0.7
    
    return Phi
