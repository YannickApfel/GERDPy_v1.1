# -*- coding: utf-8 -*-
""" Analytische Reihenlösung für den Wärmestrom in einem Rohrregister mit
    äquidistanten Rohren nach VDI 2055-1 (T-Randbedingung)

    q_l - Wärmeleistung pro Meter verbautem Rohr [W/m]

    - verwendet zur Ermittlung des thermischen Widerstands des Heizelements

    Autor: Yannick Apfel
"""
import math


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
def q_l(x_o, x_u, d_R_a, d_R_i, lambda_B, lambda_R, s, Theta_R, Theta_inf_o, state_u_insul):  # [W/m]

    # 1.) zusätzliche geometrische Größen
    s_c = 0.0  # [m] -> setze 0, da keine Zusatzschichten vorhanden
    d_insul_a = d_R_a + 0.0002  # Option: Modellierung Kontaktwiderstand als Luftschicht (1/10 mm)

    # 2.) zusätzliche therm. Größen
    Theta_inf_u = Theta_inf_o
    lambda_insul = 0.0262  # Luft
    alpha_o = 1e10  # WÜK an Oberseite --> unendlich, da T-RB
    alpha_u = alpha_o
    
    if state_u_insul:
        alpha_u = 1e-10  # WÜK an Unterseite: --> 0, isoliert
        x_u = 1e10  # x_u überschreiben mit halbunendlichem Raum

    # 3.) Wärmedurchgänge  [W/m²K]
    kappa_o = (1 / alpha_o + s_c / lambda_B) ** -1
    kappa_u = (1 / alpha_u + s_c / lambda_B) ** -1
    kappa_o_ = (kappa_o ** -1 + x_o / lambda_B) ** -1
    kappa_u_ = (kappa_u ** -1 + x_u / lambda_B) ** -1

    # 4.) Wärmestrom pro Meter [W/m]

    ssum = sum_fct(kappa_o, kappa_u, s, s_c, x_o, x_u, lambda_B)  # Ermittlung des sum-Terms

    q_l = 2 * math.pi * lambda_B * (Theta_R - (Theta_inf_o * kappa_o_ + Theta_inf_u * kappa_u_) / (kappa_o_ + kappa_u_)) \
          / (lambda_B / lambda_R * math.log(d_R_a / d_R_i)
          + lambda_B / lambda_insul * math.log(d_insul_a / d_R_a) + math.log(s / (math.pi * d_insul_a))
          + (2 * math.pi * lambda_B) / (s * (kappa_o_ + kappa_u_)) + ssum)  # [W/m]

    return q_l
