# -*- coding: utf-8 -*-
""" Ermittlung der Systemleistung anhand einer Leistungsbilanz an der Oberfläche
    des Heizelements für jeden Zeitschritt

    Q. = (T_g - T_surf) / R_th_tot = Summe Wärmeströme am Heizelement
    Nullstellenproblem:
        F(Q.) = -Q. + Summe Wärmeströme am Heizelement = 0

    Autor: Yannick Apfel
"""
import sympy as sp
import numpy as np
from alpha_konv import alpha_konv


def load(u_inf, A_he, T_b, R_th, T_inf):
    # einfaches System: nur Konvektion an der Oberfläche
    Q = sp.symbols('Q')
    F = sp.Eq(alpha_konv(u_inf) * A_he * (T_b - Q * R_th - T_inf) - Q, 0)

    Q_sol = np.array(sp.solve(F, Q))
    if Q_sol[0] < 0:  # Q. < 0 bei Gravitationswärmerohren nicht möglich
        Q_sol[0] = 0

    return Q_sol
