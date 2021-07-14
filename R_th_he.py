# -*- coding: utf-8 -*-
""" Thermischer Widerstand im Heizelement
    ermittelt aus den Ergebnissen der FEM-Simulationen von Yannick Apfel für
    die Verrohrung des Oberflächenheizelements mit Rohrschlaufen und ca. 55 mm
    Abstand zwischen den Rohren im Verteilsystem

    variable Größen:
        - minimaler Oberflächenabstand x_min [mm]
        - Wärmeleitfähigkeit des Betons des Heizelements [W/mK]

    Autor: Yannick Apfel
"""


def R_th_he(he):

    # %% 1.) Parameter

    # Fläche Heizelement [m²]
    A_he = he.A_he

    # minimaler Oberflächenabstand [mm]
    x_min = he.x_min

    # Wärmeleitfähigkeit des Betonelements [W/mk]
    lambda_Bet = he.lambda_Bet

    # %% 2.) thermischer Widerstand

    r_th_he = (5 - (- (0.55 / 5) * x_min + 2 * (lambda_Bet - 2.3) + 3.4)) / 250

    R_th_he = r_th_he / A_he

    return R_th_he
