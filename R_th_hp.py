# -*- coding: utf-8 -*-
""" GERDPy - 'R_th_hp.py'
    
    Modul für den thermischen Widerstand im Thermosiphon (Heatpipe)
    
    Thermischer Widerstand im Thermosiphon (verdampferseitig) mit
    simpler Charakteristik im Q.-deltaT - Diagramm (linearer Anstieg)
    
    Vernachlässigung von Leistungsgrenzen (z. B: Entrainment-Limit)

    Heatpipe-Charakteristik: Q. = (500 W * N_Anzahl_Heatpipes / 1 K) * deltaT

    Autor(en): Yannick Apfel
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
