# -*- coding: utf-8 -*-
""" Thermische Verluste (exkl. Oberfläche des Heizelements)

    Q_V_An: thermische Verluste der Anbindung zwischen Erdwärmesonden und Heizelement
            - Modellierung mittels Péclet-Gleichung für Zylinderschalen (Rohr + Isolierung) mit der Gesamtlänge aller Heatpipes
              im Bereich der Anbindung
            - Ermittlung der Rohrinnentemperatur über die Bohrlochtemperatur Theta_b und den therm. Widerstand R_th_g_hp
            - Annahme Theta_inf (Umgebungstemperatur) als Rohraußentemperatur
    
    Q_V_He: thermische Verluste an der Unterseite des Heizelements

    Autor: Yannick Apfel
"""
import math


def Q_V_An(Theta_R, Theta_inf, lambda_p, lambda_iso, l_An, r_iso, r_pa, r_pi):  # [W]
 
    return (Theta_R - Theta_inf) * (2 * math.pi * l_An) \
        * (math.log(r_pa / r_pi) / lambda_p + math.log(r_iso / r_pa) / lambda_iso) ** -1


def Q_V_He():
    pass