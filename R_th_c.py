# -*- coding: utf-8 -*-
""" Kontaktwiderstand Bohrloch-Hinterfüllung
    [VDI-Wärmeatlas 2013 - Wärmeübergangskoeffizient Wand-Schüttung]

    Autor: Yannick Apfel
"""
def R_th_c(borefield):

    import math, boreholes
    from scipy.constants import pi, R, Stefan_Boltzmann
    from boreholes import length_field

    # %% 1.) Stoffwerte und Parameter

    phi = 0.8               # Flächenbedeckungsgrad [-]
    lambda_g = 0.025        # Wärmeleitfähigkeit Gas [W/mK]
    d = 1e-3                # Partikeldurchmesser [m]
    delta = 250 * 1e-6      # Oberflächenrauhigkeit Partikel [m]
    C = 2.8                 # materialabhängige Konstante [-]
    M = 0.02896             # molare Masse Gas [kg/mol]
    T = 283                 # Temperatur der Kontaktzone [K]
    c_pg = 1007             # spez. Wärmekapazität Gas [J/kgK]
    epsilon_S = 0.2         # Emissionskoeffizient Schüttung [-]
    epsilon_W = 0.2         # Emissionskoeffizient Wand [-]
    p = 100000              # Gasdruck [Pa]

    # %% 2a.) Hilfsparameter

    # Akkomodationskoeffizient
    gamma = (10 ** (0.6 - (1000 / T + 1) / C) + 1) ** -1

    # freie Weglänge der Gasmoleküle

    l_frei = 2 * (2 - gamma) / gamma * math.sqrt(2 * pi * R * T / M) * \
        lambda_g / (p * (2 * c_pg - R / M))

    # Verhältnisse der Emissionsgrade
    C_WS = Stefan_Boltzmann / (1 / epsilon_W + 1 / epsilon_S - 1)

    # %% 2b.) Wärmeübergangskoeffizient Wand-Schüttung

    # Anteil Wärmeleitung
    alpha_WP = 4 * lambda_g / d * ((1 + 2 * (l_frei + delta) / d) *
                                   math.log(1 + d / (2 * (l_frei + delta)))
                                   - 1)

    # Anteil Strahlung
    alpha_rad = 4 * C_WS * T ** 3

    # Wärmeübergangskoeffizient Wand-Schüttung
    alpha_WS = phi * alpha_WP + alpha_rad

    # %% 3.) thermischer Kontaktwiderstand des Erwärmesondenfelds

    r_b = borefield[0].radius()
    H_field = length_field(borefield)

    R_th_c = (2 * pi * r_b * alpha_WS * H_field) ** -1

    return R_th_c
