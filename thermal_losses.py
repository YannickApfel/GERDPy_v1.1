# -*- coding: utf-8 -*-
""" Thermische Verluste (exkl. Oberfläche des Heizelements)

    Q_V_An: thermische Verluste der Anbindung zwischen Erdwärmesonden und Heizelement
        - Modellierung mittels Péclet-Gleichung für Zylinderschalen (Rohr + Isolierung) mit der Gesamtlänge aller Heatpipes
          im Bereich der Anbindung
        - Ermittlung der Rohrinnentemperatur über die Bohrlochtemperatur Theta_b und den therm. Widerstand R_th_g_hp
        - Annahme Theta_inf (Umgebungstemperatur) als Rohraußentemperatur

    Q_V_He: thermische Verluste an der Unterseite des Heizelements
        - "R_th_he_u": Verrohrung (Heatpipe-Innenseite) bis Heizelement-Unterseite (ohne Isolierung)
        - "R_th_he_iso": Isolationsschicht an Unterseite des Heizelements

    Autor: Yannick Apfel
"""
import math
from heating_element_utils import *


def Q_V_An(Theta_R, Theta_inf, lambda_p, lambda_iso, l_An, r_iso, r_pa, r_pi):  # [W]

    return (Theta_R - Theta_inf) * (2 * math.pi * l_An) \
        * (math.log(r_pa / r_pi) / lambda_p + math.log(r_iso / r_pa) / lambda_iso) ** -1


def Q_V_He(he, lambda_iso, Theta_R, Theta_inf):
    
    x_o = he.x_min + 0.5 * he.d_R_a  # Achsabstand Kondensatrohr nach oben [m]
    x_u = he.D - x_o  # Achsabstand Kondensatrohr nach unten [m]

    R_th_he_u = (he.l_R * q_l(x_u, x_o, he.d_R_a, he.d_R_i, he.lambda_B, he.lambda_R, he.s_R, 1, 0, state_u_insul=True)) ** -1
    # x_o und x_u sind vertauscht, da Unterseite betrachtet wird (Vergl. 'R_th_he.py')

    R_th_he_iso = 1 / lambda_iso * he.D_iso / he.A_he

    return (Theta_R - Theta_inf) * (R_th_he_u + R_th_he_iso) ** -1
