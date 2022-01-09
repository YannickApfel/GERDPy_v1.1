# -*- coding: utf-8 -*-
""" Thermische Verluste (exkl. Oberfläche des Heizelements)

    Q_V_An: thermische Verluste der Anbindung zwischen Erdwärmesonden und Heizelement
            - Modellierung mittels Péclet-Gleichung für eine Zylinderschale mit der Gesamtlänge aller Heatpipes
              im Bereich der Anbindung
            - Ermittlung der Rohrinnentemperatur über die Bohrlochtemperatur Theta_b und den therm. Widerstand R_th_g_hp
            - Annahme Theta_inf (Umgebungstemperatur) als Rohraußentemperatur
    

    Autor: Yannick Apfel
"""
import math


def Q_V_An(Theta_b, Q, R_th_ghp, Theta_inf, lambda_p, l_R_An, N, r_pa, r_pi):  # [W]

    return (Theta_b - Q * R_th_ghp - Theta_inf) * (2 * math.pi * lambda_p * N * l_R_An / math.log(r_pa / r_pi))