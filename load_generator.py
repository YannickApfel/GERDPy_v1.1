# -*- coding: utf-8 -*-
""" Ermittlung der Systemleistung anhand einer Leistungsbilanz an der Oberfläche
    des Heizelements für jeden Zeitschritt

    Q. = (T_g - Theta_surf) / R_th_tot = Summe Wärmeströme am Heizelement
    Nullstellenproblem:
        F(Q.) = Summe Wärmeströme am Heizelement - Q. = 0
        
    Legende:
        - Temperaturen:
            - T in Kelvin [K] - für (kalorische) Gleichungen
            - Theta in Grad Celsius [°C] - Input aus dem Wetterdatenfile

    Autor: Yannick Apfel
"""
import math, sys
import sympy as sp
import numpy as np
import CoolProp.CoolProp as CP  # verwendete Unterfunktionen: PropsSI, HAPropsSI (humid air)
from scipy.constants import sigma


# Definitionen

# a) Stoffwerte
lambda_l = 0.0262  # Wärmeleitfähigkeit von Luft (Normbedingungen) [W/m*K]
a_l = lambda_l / (1.293 * 1005)  # T-Leitfähigkeit von Luft (Normbedingungen) [m²/s]
beta = lambda Theta_inf: 1 / (Theta_inf + 273.15)  # isobarer therm. Ausdehnungskoeffizient f. ideale Gase [1/K]
theta_l = 13.3e-6   # kin. Vis. von Luft [m²/s]
rho_w = 997  # Dichte von Wasser (bei 25 °C) [kg/m³]
rho_l = 1.276  # Dichte trockener Luft (bei 0 °C, 1 bar) [kg/m³]
h_Ph_sl = 333e3  # Phasenwechselenthalpie Schnee <=> Wasser [J/kg]
h_Ph_lg = 2256.6e3  # Phasenwechselenthalpie Wasser <=> Dampf [J/kg]
c_p_s = 2.04e3  # spez. Wärmekapazität Eis / Schnee (bei 0 °C) [J/kgK]
c_p_w = 4212  # spez. Wärmekapazität Wasser (bei 0 °C) [J/kgK]
c_p_l = 1006  # spez. Wärmekapazität Luft (bei 0 °C, 1 bar) [J/kgK]
Theta_Schm = 0  # Schmelztemperatur von Eis / Schnee [°C]
H_max = 1  # maximal erlaubte Wasserhöhe auf dem Heizelement [mm]


# korrigierte Windgeschwindigkeit (Wind-Shear) [m/s]
def u_eff(v):
    z_1 = 10  # Höhe, für die die Windgeschwindigkeit bekannt ist (Wetterdaten) [m]
    z_n = 1  # Bezugshöhe (liegt für das Problem definitorisch bei 1 m)
    z_0 = 0.005  # Rauhigkeitshöhe
    
    return v * (math.log10(z_n) - math.log10(z_0)) / (math.log10(z_1) - math.log10(z_0))


# höhenkorrigierter Umgebungsdruck [Pa]
def p_inf(h_NHN):
    return 101325 * (1 - 0.0065 * h_NHN / 288.2) ** 5.265

# Sättigungs-Dampfdruck nach ASHRAE2013 [Pa]
def p_s(T):  # Input in [K]
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
        return math.exp(C_1 / T + C_2 + C_3 * T + C_4 * T ** 2 + C_5 * T ** 3 + C_6 * T ** 4 + C_7 * math.log(T))
    elif ((T - 273.15) >= 0) and ((T - 273.15) <= 200):  # "over liquid water"
        return math.exp(C_8 / T + C_9 + C_10 * T + C_11 * T ** 2 + C_12 * T ** 3 + C_13 * math.log(T))
    else:
        print('Interner Fehler: erlaubter T-Bereich für Sättigungs-Dampfdruckformel nach ASHRAE2013 unter-/überschritten!')
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
def m_Restwasser(m_Rw, RR, Q_eva, A_he):
    if m_Rw < 0:  # Restwassermenge kann nicht geringer werden als 0
        m_Rw_sol = 0
    else:
        m_Rw_sol = m_Rw + (RR * rho_w * 1) / 1000 - (Q_eva / h_Ph_lg) * 3600
    
    if (m_Rw_sol / (rho_w * A_he)) > (H_max / 1000):  # H_max deckelt die max. mögl. Wasserhöhe
        m_Rw_sol = (H_max / 1000) * rho_w * A_he
    
    return m_Rw_sol


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
    else:  # ohne Schneefall: als Funtion von Umgebungstemp. und rel. Luftfeuchte
        T_H = (Theta_inf + 273.15) - (1.1058e3 - 7.562 * (Theta_inf + 273.15) +
                                  1.333e-2 * (Theta_inf + 273.15) ** 2 - 31.292 *
                                  Phi + 14.58 * Phi ** 2)
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
    # Sättigungsdampfdruck in der Umgebung bei Taupunkttemperatur: p_D = p_s(T_tau(Theta_inf, Phi))
    T_tau = CP.HAPropsSI('DewPoint', 'T', (Theta_inf + 273.15), 'P', 101325, 'R', Phi)  # Input in [K]
    p_D = p_s(T_tau)  # Input in [K]
    
    return 0.622 * p_D / (p_inf(h_NHN) - p_D)


# Wasserdampfbeladung der gesättigten Luft bei Theta_surf [kg Dampf / kg Luft]
def X_D_sat_surf(Theta_surf, h_NHN):
    # Sättigungsdampfdruck an der Heizelementoberfläche bei Theta_surf: p_D = p_s(Theta_surf)
    p_D = p_s(Theta_surf + 273.15)  # Input in [K]

    return 0.622 * p_D / (p_inf(h_NHN) - p_D)


# Definition & Bilanzierung der Einzellasten
def load(h_NHN, v, Theta_inf, S_w, A_he, Theta_b_0, R_th, Theta_surf_0, B, Phi, RR, m_Rw_0):  # Theta_x_0: Temp. des vorhergehenden Zeitschritts
                                                                   # Input-Temperaturen in [°C]
        
    # 0.) Preprocessing
    u_inf = u_eff(v)  # Reduzierte Windgeschwindigkeit (logarithmisches Windprofil)
    Q = sp.symbols('Q')  # thermische Leistung Q als Variable definieren
    
    if Theta_surf_0 < 0 and S_w > 0:  # snow-free area ratio dynamisch
        R_f = 0.25  # Schneedecke bildet sich, Schnee nicht ganz "ideal isolierend"
    else:
        R_f = 1  # Oberfläche ist schnee-frei (alle Verlustmechanismen treffen maximal ein)
    
    # 1.) Teil-Wärmeströme aktivieren oder deaktivieren (für unit-testing))
    lat = False
    sen = False
    con = False
    rad = False
    eva = True

    # Q_latent
    if lat is True:
        Q_lat = rho_w * S_w * h_Ph_sl * (3.6e6)**-1 * A_he
    else:
        Q_lat = 0

    # Q_sensibel
    if sen is True:
        Q_sen = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_b_0 - Q * R_th - Theta_Schm)) * (3.6e6)**-1 * A_he
    else:
        Q_sen = 0

    # Q_Konvektion
    if con is True:
        Q_con = alpha_kon_Bentz(u_inf) * (Theta_b_0 - Q * R_th - Theta_inf) * A_he
    else:
        Q_con = 0

    # Q_Strahlung
    if rad is True:  # Theta_surf_0 statt Theta_b_0 - Q * R_th, da sonst 4 Nullstellen
        Q_rad = sigma * epsilon_surf('Beton') * ((Theta_surf_0 + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he
    else:
        Q_rad = 0

    # Q_Verdunstung
    ''' Voraussetzungen: 
            - Theta_surf >= 0 °C
            - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
    '''    
    if (eva is True and Theta_surf_0 >= 0 and m_Rw_0 > 0):
        Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf_0, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he
    else:
        Q_eva = 0
        
    if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
        Q_eva = 0

    # Ermittlung Restwassermenge
    m_Rw_1 = m_Restwasser(m_Rw_0, RR, Q_eva, A_he)

    # 2.) stationäre Leistungbilanz am Heizelement (Kopplung Oberfläche mit Erdboden)
    F_Q = sp.Eq(Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva) - Q, 0)

    # 3.) Auflösung der Leistungsbilanz nach Q
    Q_sol = float(np.array(sp.solve(F_Q, Q)))
    
    # 4.) Neuermittlung der Oberflächentemperatur (reduzierte Leistungsbilanz), falls Leistung negativ (Umgebung erwärmt Heizelement)
    net_neg = False  # Variable zur Übergabe an _main
    Theta_surf_sol = 0  # Initialisierung der neuen Oberflächentemperatur
    
    if Q_sol < 0:
        net_neg = True  # im Falle eines negativen Wärmestroms wahr
        Theta_surf_ = sp.symbols('Theta_surf_')  # Oberflächentemperatur als Variable definieren

        # Q_latent
        Q_lat_red = Q_lat

        # Q_sensibel
        if sen is True:
            Q_sen_red = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_surf_ - Theta_Schm)) * (3.6e6)**-1 * A_he
        else:
            Q_sen_red = 0

        # Q_Konvektion
        if con is True:
            Q_con_red = alpha_kon_Bentz(u_inf) * (Theta_surf_ - Theta_inf) * A_he
        else:
            Q_con_red = 0

        # Q_Strahlung
        if rad is True:  # Theta_surf_0 statt Theta_surf_, da sonst 4 Nullstellen
            Q_rad_red = sigma * epsilon_surf('Beton') * ((Theta_surf_0 + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he
        else:
            Q_rad_red = 0

        # Q_Verdunstung
        Q_eva_red = Q_eva

        # reduzierte Leistungsbilanz
        F_Q_red = sp.Eq(Q_lat_red + Q_sen_red + R_f * (Q_con_red + Q_rad_red + Q_eva_red), 0)
        
        Theta_surf_sol = float(np.array(sp.solve(F_Q_red, Theta_surf_)))
    
    # 4.) return an _main
    if Q_sol < 0:  # Q. < 0 bei Gravitationswärmerohren nicht möglich
        Q_sol = 0

    return Q_sol, net_neg, Theta_surf_sol, m_Rw_1
