# -*- coding: utf-8 -*-
""" Thermischer Widerstand des Heizelements

    1.) R_th_he_PS: (PS - Prüfstand)

        - thermischer Widerstand des Heizelements des ZAE Klimaprüfstands

        - parametrierte Geradengleichung aus numerischen Simulationen (FEM) des
        Oberflächenheizelements für die Klimakammer des ZAE, basierend auf der
        Masterarbeit [Apfel2020]; Heizleistung: 250 W, delta-T: 5 K

        - Geometrie: Rohrschlaufen mit ca. 55 mm Abstand der Kondensatorrohre

        - variable Größen:
            - Fläche des Heizelements [m²]
            - minimaler Oberflächenabstand x_min [m]
            - Wärmeleitfähigkeit des Betons des Heizelements [W/mK]

    2.) R_th_he_an: (an - analytisch))

        - analytisch ermittelter therm. Widerstand basierend auf der Norm
        VDI 2055-1 [2019] für die thermischen Verluste gedämmter Rohrleitungen
        verlegt in Registeranordnung in ein- oder mehrschichtigen Platten
        (z. B. Fußbodenheizungen) - vereinfacht und angepasst in [Apfel2020]
        
        - variable Größen:
            - minimaler Oberflächenabstand x_min [m]
            - Wärmeleitfähigkeiten Beton und Kondensatrohre [W/mK]
            - Rohrgeometrie des Kondensator (Durchmesser, Rohrabstand,
                                             Rohrlänge) [m]

    Autor: Yannick Apfel
"""
from heating_element_utils import *


def R_th_he_PS(he):  # thermischer Widerstand [K/W]

    r_th_he = (5 - (- (0.55 / 5) * he.x_min * 1000 + 2 * (he.lambda_B - 2.3) +
                    3.4)) / 250  # [Km²/W]

    R_th_he = r_th_he / he.A_he  # [K/W]

    return R_th_he


def R_th_he_an(he):  # thermischer Widerstand [K/W]

    # 1.) zusätzliche Parameter des Heizelements
    x_o = he.x_min + 0.5 * he.d_R_a  # Achsabstand Kondensatrohr nach oben [m]
    x_u = he.D - x_o  # Achsabstand Kondensatrohr nach unten [m]

    # 2.) delta-T zwischen Kondensatrohren und Oberfläche zu 1 K definieren
    Theta_R = 1
    Theta_inf_o = 0

    # 3.) thermischer Widerstand des Heizelements in [K/W] --> R_th = 1 K / (q_l * l_R)
    ''' Unterseite wird als thermisch isoliert betrachtet: state_u_insul=True
    '''
    R_th_he = (he.l_R * q_l(x_o, x_u, he.d_R_a, he.d_R_i, he.lambda_B, he.lambda_R, he.s_R, Theta_R, Theta_inf_o, state_u_insul=True)) ** -1

    return R_th_he
