# -*- coding: utf-8 -*-
""" Thermischer Widerstand des Heizelements

    1.) R_th_he_PS: (PS - Prüfstand))

        - parametrierte Geradengleichung aus numerischen Simulationen (FEM) des
        Oberflächenheizelements für die Klimakammer des ZAE, basierend auf der
        Masterarbeit [Apfel2020]; Heizleistung: 250 W, delta-T: 5 K

        - Geometrie: Rohrschlaufen mit ca. 55 mm Abstand der Kondensatorrohre

        - variable Größen:
            - minimaler Oberflächenabstand x_min [mm]
            - Wärmeleitfähigkeit des Betons des Heizelements [W/mK]

    2.) R_th_he_an: (an - analytisch))

        - analytisch ermittelter therm. Widerstand basierend auf der Norm
        VDI 2055-1 [2019] für die thermischen Verluste gedämmter Rohrleitungen
        verlegt in Registeranordnung in ein- oder mehrschichtigen Platten
        (z. B. Fußbodenheizungen) - vereinfacht und angepasst in [Apfel2020]

    Autor: Yannick Apfel
"""
import math
import numpy as np

import matplotlib.pyplot as plt


def R_th_he_PS(he):  # thermischer Widerstand [K/W]

    r_th_he = (5 - (- (0.55 / 5) * he.x_min * 1000 + 2 * (he.lambda_Bet - 2.3) +
                    3.4)) / 250  # [Km²/W]

    R_th_he = r_th_he / he.A_he  # [K/W]

    return R_th_he


# Ermittlung des sum-Terms der analyt. Lösung aus VDI 2055-1
def sum_fct(kappa_o, kappa_u, s, s_c, x_o, x_u, lambda_B):

    # 1.) Definition der therm. Größen
    beta_o = kappa_o * s / lambda_B
    beta_u = kappa_u * s / lambda_B

    # 2.) Approximation der unendlichen Reihensumme ssum
    ssum_temp = 0
    j = 0
    error = 1e-6

    while 1:  # unendliche Schleife, bis break-Befehl
    
        j += 1

        N_1 = 1 - (beta_u + 2 * math.pi * j) / (beta_u - 2 * math.pi * j) * math.exp(4 * math.pi * j * s_c / s)
        N_2 = 1 - (beta_u - 2 * math.pi * j) / (beta_u + 2 * math.pi * j) * math.exp(-4 * math.pi * j * s_c / s)

        gamma = (beta_o - 2 * math.pi * j) / (beta_o + 2 * math.pi * j) * math.exp(-4 * math.pi * j * (x_u + x_o) / s)

        e_o = ((lambda_B + lambda_B / N_1 - lambda_B / N_2) * (math.exp(-4 * math.pi * j * x_u / s) - gamma)) /\
              (lambda_B * (1 + gamma) + (lambda_B / N_2 - lambda_B / N_1) * (1 - gamma))

        e_u = - (beta_o - 2 * math.pi * j) / (beta_o + 2 * math.pi * j) * math.exp(-4 * math.pi * j * x_o / s) * (1 + e_o)
        
        ssum = ssum_temp + (e_o + e_u) / j
        
        if abs(ssum - ssum_temp) < error:
            break
        
        ssum_temp = ssum
    
    return ssum


# analyt. Lösung aus VDI 2055-1 (Wärmeleistung pro Rohrmeter)
def q_l():  # [W/m]
    
    # 1.) geometrische Größen
    s_c = 0.0  # [m] -> setze 0, da keine Zusatzschichten vorhanden
    x_o = 0.034
    x_u = 1e10  # [m] -> Wärmeverluste einseitig (halbunendlicher Raum)
    d_R_a = 0.008
    d_R_i = 0.006
    
    d_insul_a = d_R_a  # Option: Modellierung Kontaktwiderstand als Isolationsschicht
    
    # 2.) therm. Größen
    delta_T = 4.5
    lambda_B = 2.1
    lambda_R = 14
    lambda_insul = 0.035
    alpha_o = 1e10  # großer WÜK an der Oberfläche, da T-RB.
    alpha_u = 1e-10  # andere Seite ist quasi isoliert
    
    # 3.) Wärmedurchgänge  [W/m²K]
    kappa_o = (1 / alpha_o + s_c / lambda_B) ** -1
    kappa_u = (1 / alpha_u + s_c / lambda_B) ** -1
    kappa_o_ = (kappa_o ** -1 + x_o / lambda_B) ** -1
    kappa_u_ = (kappa_u ** -1 + x_u / lambda_B) ** -1
    
    # 4.) Wärmestrom pro Meter [W/m]
    s_vec = np.arange(d_insul_a, 1, .001)
    ssum = np.zeros(len(s_vec))
    q_l = np.zeros(len(s_vec))
    
    for i in range(len(s_vec)):
        
        ssum[i] = sum_fct(kappa_o, kappa_u, s_vec[i], s_c, x_o, x_u, lambda_B)  # Ermittlung des sum-Terms
        
        q_l[i] = 2 * math.pi * lambda_B * delta_T / (lambda_B / lambda_R * math.log(d_R_a / d_R_i) \
                    + lambda_B / lambda_insul * math.log(d_insul_a / d_R_a) + math.log(s_vec[i] / (math.pi * d_insul_a)) \
                    + (2 * math.pi * lambda_B) / (s_vec[i] * (kappa_o_ + kappa_u_)) + ssum[i])  # [W/m]
                                                
    return q_l, s_vec


def R_th_he_an():  # thermischer Widerstand [K/W]

    # Parameter Heizelement
    # x_o = he.x_min + 0.5 * he.d_R_a  # Achsabstand Kondensatrohr nach oben [m]
    # x_u = he.D - x_o  # Achsabstand Kondensatrohr nach unten [m]
                                                    
    plt.plot(q_l()[1], q_l()[0])
    plt.grid('major')
    
    print(q_l()[0][3])


        
        