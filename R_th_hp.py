# -*- coding: utf-8 -*-
""" Thermischer Widerstand im Thermosiphon (verdampferseitig) mit
    simpler Charakteristik im Q.-deltaT - Diagramm (linearer Anstieg)

    Heatpipe-Charakteristik: Q. = (500W / 1K) * deltaT (ohne Entrainment)

    Autor: Yannick Apfel
"""


def R_th_hp(borefield, hp):

    # %% 1.) Parameter

    N = hp.N                                # Anzahl Heatpipes pro Bohrloch [-]
    N_b = len(borefield)                    # Anzahl Bohrlöcher

    # Steigung der Charakteristik (delta_y / delta_x)
    delta_y = 500       # Leistung pro Heatpipe [W]
    delta_x = 1         # delta_T [K]

    # %% 2.) thermischer Widerstand (Entrainment wird vernachlässigt)

    R_th_hp = delta_x / (delta_y * N * N_b)

    return R_th_hp
