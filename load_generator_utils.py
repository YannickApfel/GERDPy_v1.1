# -*- coding: utf-8 -*-
""" Stoffwerte und Unterfunktionen für 'load_generator.py'

    Legende:
        - Temperaturen:
            - T in Kelvin [K] - für (kalorische) Gleichungen
            - Theta in Grad Celsius [°C] - Input aus dem Wetterdatenfile

    Quellen: [Konrad2009], [Fuchs2021], [ASHRAE2015]

    Autor: Yannick Apfel
"""
import math, sys
import CoolProp.CoolProp as CP


# Stoffwerte (Normbedingungen, wo nicht anders vermerkt)
lambda_l = 0.0262  # Wärmeleitfähigkeit von Luft [W/m*K]
a_l = lambda_l / (1.293 * 1005)  # T-Leitfähigkeit von Luft [m²/s]
theta_l = 13.3e-6  # kin. Vis. von Luft [m²/s]
rho_w = 997  # Dichte von Wasser (bei 25 °C) [kg/m³]
rho_l = 1.33  # Dichte trockener Luft (bei 0 °C, 1 bar) [kg/m³]
h_Ph_sl = 333e3  # Phasenwechselenthalpie Schnee <=> Wasser [J/kg]
h_Ph_lg = 2499e3  # Phasenwechselenthalpie Wasser <=> Dampf [J/kg]
c_p_s = 2.04e3  # spez. Wärmekapazität Eis / Schnee (bei 0 °C) [J/kgK]
c_p_w = 4212  # spez. Wärmekapazität Wasser (bei 0 °C) [J/kgK]
c_p_l = 1005  # spez. Wärmekapazität Luft (bei 0 °C, 1 bar) [J/kgK]
Theta_Schm = 0  # Schmelztemperatur von Eis / Schnee [°C]
H_max = 2  # maximal erlaubte Wasserhöhe auf dem Heizelement [mm]


# korrigierte Windgeschwindigkeit (Wind-Shear) [m/s]
def u_eff(v):
    z_1 = 10  # Höhe, für die Windgeschwindigkeit bekannt ist (Wetterdaten) [m]
    z_n = 1  # Bezugshöhe (liegt für das Problem definitorisch bei 1 m)
    z_0 = 0.005  # Rauhigkeitshöhe

    return v * (math.log10(z_n) - math.log10(z_0)) / \
        (math.log10(z_1) - math.log10(z_0))


# höhenkorrigierter Umgebungsdruck [Pa]
def p_inf(h_NHN):
    return 101325 * (1 - 0.0065 * h_NHN / 288.2) ** 5.265


# Sättigungs-Dampfdruck nach ASHRAE2013 [Pa]
def p_s_ASHRAE(T):  # Input in [K]
    C_1 = -5.6745359e3
    C_2 = 6.3925247e0
    C_3 = -9.6778430e-3
    C_4 = 6.2215701e-7
    C_5 = 2.0747825e-9
    C_6 = -9.4840240e-13
    C_7 = 4.1635019e0
    C_8 = -5.8002206e3
    C_9 = 1.3914993e0
    C_10 = -4.8640239e-2
    C_11 = 4.1764768e-5
    C_12 = -1.4452093e-8
    C_13 = 6.5459673e0

    if ((T - 273.15) > -100) and ((T - 273.15) < 0):  # "over ice"
        return math.exp(C_1 / T + C_2 + C_3 * T + C_4 * T ** 2 + C_5 * T ** 3
                        + C_6 * T ** 4 + C_7 * math.log(T))
    elif ((T - 273.15) >= 0) and ((T - 273.15) <= 200):  # "over liquid water"
        return math.exp(C_8 / T + C_9 + C_10 * T + C_11 * T ** 2
                        + C_12 * T ** 3 + C_13 * math.log(T))
    else:
        print('Interner Fehler: erlaubter T-Bereich für \
              Sättigungs-Dampfdruckformel nach ASHRAE2013 \
                  unter-/überschritten!')
        sys.exit()


# Wärmeübergangskoeffizient [W/m²K]
# nach [Bentz D. P. 2000] (nur erzwungene Konvektion)
# alpha = alpha(u_air)
def alpha_kon_Bentz(u):
    if u <= 5:
        alpha = 5.6 + 4 * u
    else:
        alpha = 7.2 * u ** 0.78

    return alpha


# Wassermengen-Bilanz an Heizelement-Oberfläche (für Verdunstung)
def m_Restwasser(m_Rw_0, RR, A_he, Q_eva):
    m_Rw_1 = m_Rw_0 + (RR * rho_w * A_he) / 1000 - (Q_eva / h_Ph_lg) * 3600

    if (m_Rw_1 / (rho_w * A_he)) > (H_max / 1000):  # Maximale Wasserhöhe
        m_Rw_1 = (H_max / 1000) * rho_w * A_he

    if m_Rw_1 < 0:  # Restwassermenge kann nicht geringer werden als 0
        m_Rw_1 = 0

    return m_Rw_1


# Schneemengen-Bilanz an Heizelement-Oberfläche
def m_Restschnee(m_Rs_0, S_w, A_he, Q_lat, sb_active):
    if (sb_active == 1):
        m_Rs_1 = m_Rs_0 + (S_w * rho_w * A_he) / 1000 - (Q_lat / h_Ph_sl) \
            * 3600
    else:
        m_Rs_1 = 0

    if m_Rs_1 < 0:  # Restschneemenge kann nicht geringer werden als 0
        m_Rs_1 = 0

    return m_Rs_1


# Emissionskoeffizient des Heizelements [-]
def epsilon_surf(material):
    if material == 'Beton':
        return 0.94
    else:
        return 0.94


# mittlere Strahlungstemperatur der Umgebung [K]
def T_MS(S_w, Theta_inf, B, Phi):
    if S_w > 0:  # entspricht bei Schneefall der Umgebungstemperatur
        return (Theta_inf + 273.15)
    else:  # ohne Schneefall: Funktion von Umgebungstemp. und rel. Luftfeuchte
        T_H = (Theta_inf + 273.15) - (1.1058e3 - 7.562 * (Theta_inf + 273.15) +
                                      1.333e-2 * (Theta_inf + 273.15) ** 2
                                      - 31.292 * Phi + 14.58 * Phi ** 2)
        T_W = (Theta_inf + 273.15) - 19.2

        if T_H > T_W:
            T_H = T_W

        T_MS = (T_W ** 4 * B + T_H ** 4 * (1 - B)) ** 0.25

        return T_MS


# binärer Diffusionskoeffizient [-]
def delta(Theta_inf, h_NHN):

    return (2.252 / p_inf(h_NHN)) * ((Theta_inf + 273.15) / 273.15) ** 1.81


# Stoffübergangskoeffizient [m/s]
def beta_c(Theta_inf, u, h_NHN):
    Pr = theta_l / a_l  # Prandtl-Zahl für Luft
    Sc = theta_l / delta(Theta_inf, h_NHN)  # Schmidt-Zahl
    alpha = alpha_kon_Bentz(u)  # Verwendung der einfachen Korrelation f alpha

    beta_c = (Pr / Sc) ** (2/3) * alpha / (rho_l * c_p_l)

    return beta_c


# Wasserdampfbeladung der gesättigten Luft bei Theta_inf [kg Dampf / kg Luft]
def X_D_inf(Theta_inf, Phi, h_NHN):
    ''' Sättigungsdampfdruck in der Umgebung bei Taupunkttemperatur:
        p_D = p_s_ASHRAE(T_tau(Theta_inf, Phi))
    '''
    T_tau = CP.HAPropsSI('DewPoint', 'T', (Theta_inf + 273.15),
                         'P', 101325, 'R', Phi)  # Input in [K]
    p_D = p_s_ASHRAE(T_tau)  # Input in [K]

    return 0.622 * p_D / (p_inf(h_NHN) - p_D)


# Wasserdampfbeladung der gesättigten Luft bei Theta_surf [kg Dampf / kg Luft]
def X_D_sat_surf(Theta_surf, h_NHN):
    ''' Sättigungsdampfdruck an der Heizelementoberfläche bei Theta_surf:
        p_D = p_s_ASHRAE(Theta_surf)
    '''
    p_D = p_s_ASHRAE(Theta_surf + 273.15)  # Input in [K]

    return 0.622 * p_D / (p_inf(h_NHN) - p_D)
