# -*- coding: utf-8 -*-
""" Wärmeübergangskoeffizient (konvektiv) aufgrund erzwungener Konvektion
    nach [Bentz D. P. 2000]
    alpha = alpha(u_air)

    Autor: Yannick Apfel
"""


def alpha_konv(u_air):  # alpha [W/m²K], u_air [m/s]
    if u_air <= 5:
        alpha_konv = 5.6 + 4 * u_air
    else:
        alpha_konv = 7.2 * u_air ** 0.78

    return alpha_konv
