# -*- coding: utf-8 -*-
""" Ermittlung der Systemleistung anhand einer Leistungsbilanz an der Oberfläche
    des Heizelements für jeden Zeitschritt

    Q. = (T_g - T_surf) / R_th_tot = Summe Wärmeströme am Heizelement
    Nullstellenproblem:
        F(Q.) = Summe Wärmeströme am Heizelement - Q. = 0

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
a_l = lambda_l / (1.293 * 1005)  # T-Leitfähigkeit von Luft (Normbedingungen) [m2/s]
beta = lambda T_inf: 1 / (T_inf + 273.15)  # isobarer therm. Ausdehnungskoeffizient f. ideale Gase [1/K]
theta_l = 13.3e-6   # kin. Vis. von Luft [m²/s]
rho_w = 997  # Dichte von Wasser (bei 25 °C) [kg/m³]
rho_l = 1.276  # Dichte trockener Luft (bei 0 °C, 1 bar) [kg/m³]
h_Ph_sl = 333e3  # Phasenwechselenthalpie Schnee <=> Wasser [J/kg]
h_Ph_lg = 2256.6e3  # Phasenwechselenthalpie Wasser <=> Dampf [J/kg]
c_p_s = 2.04e3  # spez. Wärmekapazität Eis / Schnee (bei 0 °C) [J/kgK]
c_p_w = 4212  # spez. Wärmekapazität Wasser (bei 0 °C) [J/kgK]
c_p_l = 1006  # spez. Wärmekapazität Luft (bei 0 °C, 1 bar) [J/kgK]
T_Schm = 0  # Schmelztemperatur von Eis / Schnee [°C]

# b) Snow-free area ratio
R_f = 1


# korrigierte Windgeschwindigkeit (Wind-Shear)
def wind_speed_eff():
    pass


# höhenkorrigierter Umgebungsdruck
def p_inf(h_NHN):
    return 101325 * (1 - 0.0065 * h_NHN / 288.2) ** 5.265


# Wärmeübergangskoeffizient
# nach [Bentz D. P. 2000] (nur erzwungene Konvektion)
# alpha = alpha(u_air)
def alpha_kon_Bentz(u):  # alpha [W/m²K]

    if u <= 5:
        alpha = 5.6 + 4 * u
    else:
        alpha = 7.2 * u ** 0.78

    return alpha


# Wärmeübergangskoeffizient
# nach [VDI2013] (erwungene und freie Konvektion)
# alpha = alpha(Pr, Gr, Re)
def alpha_kon_VDI(A_he, u, T_inf, T_surf):  # alpha [W/m²K]

    # 1.) elementare charakteristische Größen des Strömungsproblems: L_c, Pr, Re, Gr
    L_c_er = math.sqrt(A_he)  # charakteristisches Längenmaß erzwungene Konvektion [m] (quadratisches He)
    L_c_fr = A_he / (4 * math.sqrt(A_he))  # charakteristisches Längenmaß freie Konvektion [m] (quadratisches He)

    Pr = theta_l / a_l  # Prandtl-Zahl für Luft
    Re = u * L_c_er / theta_l  # Reynolds-Zahl (Turbulenzmaß für erzwungene Konvektion)

    Gr = L_c_fr ** 3 * 9.81 * beta(T_inf) * (T_surf - T_inf) / (theta_l ** 2)  # Grashof-Zahl (Turbulenzmaß für freie Konvektion)

    # 2.) Verhältnis (Gr:Re²) als Kriterium für die Fallunterscheidung nach: freier Konvektion, Mischkonvektion, erzwungener Konvektion
    if (Gr / (Re ** 2)) > 10:  # nur freie Konvektion: Vernachlässigung erzwungene Konvektion
        kon_fr = 1
        kon_er = 0
    elif (Gr / (Re ** 2)) < 0.1:  # nur erzwungene Konvektion: Vernachlässigung freie Konvektion
        kon_fr = 0
        kon_er = 1
    else:  # Mischkonvektion
        kon_fr = 1
        kon_er = 1

    if (kon_fr == 1):  # nur freie Konvektion: Vernachlässigung erzwungene Konvektion
        L_c = L_c_fr
        Ra = Gr * Pr  # Rayleigh-Zahl
        f1_Pr = (1 + (0.322 / Pr) ** (11 / 20)) ** (-20/11)  # Prandtl-Korrekturfaktor

        if (Ra * f1_Pr) <= 7e4:
            Nu_fr = 0.766 * (Ra * f1_Pr) ** (1/5)  # laminare Nusselt-Zahl (freie Konvektion)
        else:
            Nu_fr = 0.15 * (Ra * f1_Pr) ** (1/3)  # turbulente Nusselt-Zahl (freie Konvektion)

        Nu = Nu_fr

    if (kon_er == 1):  # nur erzwungene Konvektion: Vernachlässigung freie Konvektion
        L_c = L_c_er

        Nu_er_lam = 0.664 * math.sqrt(Re) * Pr ** (1/3)  # laminare Nusselt-Zahl (erzwungene Konvektion)
        Nu_er_tur = 0.037 * Re ** 0.8 * Pr / (1 + 2.443 * Re ** -0.1 * (Pr ** (2/3) - 1))  # turb. Nusselt-Zahl (erzwungene Konvektion) 

        if ((50 < Re) and (Re < 10e7)) and ((0.5 < Pr) and (Pr < 2000)):  # Probe auf gültigen Bereich für Strömungsmilieu
            Nu_er_0 = math.sqrt(Nu_er_lam ** 2 + Nu_er_tur ** 2)  # mittlere Nusselt-Zahl für Mischbedingungen
        else:
            print("Prandtl- und/oder Reynoldszahl nicht im gültigen Bereich")
            sys.exit()

        Nu_er = ((T_inf + 273.15) / (T_surf + 273.15)) ** 0.12 * Nu_er_0  # mittlere Nu - korrigiert für T-abhängige Stoffwerte
        Nu = Nu_er

    if (kon_fr == 1 and kon_er == 1):  # Mischkonvektion
        L_c = L_c_er
        L_c_mix = L_c_er  # Definition der char. Länge entspricht der der erzw. Konvektion (beide Konvektionsarten)
        Gr_mix = L_c_mix ** 3 * 9.81 * beta(T_inf) * (T_surf - T_inf) / (theta_l ** 2)

        if (0 < Pr) and (Pr < 100):  # Probe auf gültigen Bereich für Strömungsmilieu
            Nu_fr_mix = ((Pr / 5) ** (1/5)) * (Pr ** 0.5) / (0.25 + 1.6 * Pr ** 0.5) * Gr_mix ** (1/5)
        else:
            print("Prandtl- und/oder Reynoldszahl nicht im gültigen Bereich")
            sys.exit()

        Nu_mix = (Nu_er ** (7/2) + Nu_fr_mix ** (7/2)) ** (2/7)
        Nu = Nu_mix

    # 3.) Ermittlung des WÜK [W/m²K]
    alpha = (Nu * lambda_l) / L_c

    return alpha


# Emissionskoeffizient des Heizelements
def epsilon_surf(material):
    if material == 'Beton':
        return 0.94
    else:
        return 0.94
    

# mittlere Strahlungstemperatur der Umgebung [K]
def T_MS(S_w, T_inf, B, Phi):
    if S_w > 0:  # entspricht bei Schneefall der Umgebungstemperatur
        return (T_inf + 273.15)
    else:  # ohne Schneefall: als Funtion von Umgebungstemp. und rel. Luftfeuchte
        T_H = (T_inf + 273.15) - (1.1058e3 - 7.562 * (T_inf + 273.15) +
                                  1.333e-2 * (T_inf + 273.15) ** 2 - 31.292 *
                                  Phi + 14.58 * Phi ** 2)
        T_W = (T_inf + 273.15) - 19.2

        if T_H > T_W:
            T_H = T_W

        T_MS = (T_W ** 4 * B + T_H ** 4 * (1 - B)) ** 0.25

        return T_MS
    
    
# binärer Diffusionskoeffizient
def delta(T_inf, h_NHN):

    return (2.252 / p_inf(h_NHN)) * ((T_inf + 273.15) / 273.15) ** 1.81


# Stoffübergangskoeffizient [m/s]
def beta_c(T_inf, u, h_NHN):
    Pr = theta_l / a_l  # Prandtl-Zahl für Luft
    Sc = theta_l / delta(T_inf, h_NHN)  # Schmidt-Zahl
    alpha = alpha_kon_Bentz(u)  # Verwendung der einfachen Korrelation f alpha

    beta_c = (Pr / Sc) ** (2/3) * alpha / (rho_l * c_p_l)
    
    return beta_c


# Wasserdampfbeladung der gesättigten Luft bei T_inf [kg Dampf / kg Luft]
def X_D_inf(T_inf, Phi, h_NHN):
    # Sättigungsdampfdruck in der Umgebung bei Taupunkttemperatur: p_D = p_s(T_tau(T_inf, Phi))
        # der zul. Wertebereich von PropsSI aus CoolProp nicht ausreichend....
        # ....Verwendung einer Korrelation für den Sättigungsdruck nach [Konrad2009, S. 79]
    T_tau = CP.HAPropsSI('DewPoint', 'T', (T_inf + 273.15), 'P', 101325, 'R', Phi)  # Input in [K]
    p_D = 6.108 * math.exp((17.081 * (T_tau - 273.15)) / (234.175 + (T_tau - 273.15))) * 100  # Input in [°C]
    
    return 0.622 * p_D / (p_inf(h_NHN) - p_D)



# Wasserdampfbeladung der gesättigten Luft bei T_surf [kg Dampf / kg Luft]
def X_D_sat_surf(T_surf, h_NHN):
    # Sättigungsdampfdruck an der Heizelementoberfläche bei T_surf: p_D = p_s(T_surf)
        # der zul. Wertebereich von PropsSI aus CoolProp nicht ausreichend....
        # ....Verwendung einer Korrelation für den Sättigungsdruck nach [Konrad2009, S. 79]
    p_D = 6.108 * math.exp((17.081 * T_surf) / (234.175 + T_surf)) * 100  # Input in [°C]

    return 0.622 * p_D / (p_inf(h_NHN) - p_D)


# Definition & Bilanzierung der Einzellasten
def load(h_NHN, u_inf, T_inf, S_w, A_he, T_b_0, R_th, T_surf_0, B, Phi):  # T_x_0: Temp. des vorhergehenden Zeitschritts
                                                                   # Input-Temperaturen in [°C]
                                                                   
    # 0.)
    Q = sp.symbols('Q')  # thermische Leistung Q als Variable definieren
    
    # 1.) Teil-Wärmeströme aktivieren oder deaktivieren
    lat = True
    sen = True
    kon = True
    sstr = True
    eva = True

    # Q_latent
    if lat is True:
        Q_lat = rho_w * S_w * h_Ph_sl * (3.6e6)**-1 * A_he
    else:
        Q_lat = 0

    # Q_sensibel
    if sen is True:
        Q_sen = rho_w * S_w * (c_p_s * (T_Schm - T_inf) + c_p_w * (T_b_0 - Q * R_th - T_Schm)) * (3.6e6)**-1 * A_he
    else:
        Q_sen = 0

    # Q_Konvektion
    if kon is True:
        Q_kon = alpha_kon_VDI(A_he, u_inf, T_inf, T_surf_0) * (T_b_0 - Q * R_th - T_inf) * A_he
    else:
        Q_kon = 0

    # Q_Strahlung
    if sstr is True:  # T_surf_0 statt T_b_0 - Q * R_th, da sonst 4 Nullstellen
        Q_str = sigma * epsilon_surf('Beton') * ((T_surf_0 + 273.15) ** 4 - T_MS(S_w, T_inf, B, Phi) ** 4) * A_he
    else:
        Q_str = 0

    # Q_Verdunstung
    if eva is True:  # T_surf_0 statt T_b_0 - Q * R_th
        if T_surf_0 > 0:  # bei Oberflächentemp. <= 0 °C keine Verdunstung
            Q_eva = rho_l * beta_c(T_inf, u_inf, h_NHN) * (X_D_sat_surf(T_surf_0, h_NHN) - X_D_inf(T_inf, Phi, h_NHN)) * h_Ph_lg * A_he
        else:
            Q_eva = 0
    else:
        Q_eva = 0

    # 2.) stationäre Leistungbilanz am Heizelement (Kopplung Oberfläche mit Erdboden)
    F_Q = sp.Eq(Q_lat + Q_sen + R_f * (Q_kon + Q_str + Q_eva) - Q, 0)

    # 3.) Auflösung der Leistungsbilanz nach Q
    Q_sol = float(np.array(sp.solve(F_Q, Q)))
    net_neg = False  # Variable zur Übergabe an _main
    T_surf_sol = 0  # Initialisierung der neuen Oberflächentemperatur
    
    if Q_sol < 0:  # Neuermittlung der Oberflächentemperatur (reduzierte Leistungsbilanz), falls Leistung negativ (Umgebung erwärmt Heizelement)
        net_neg = True  # im Falle eines negativen Wärmestroms wahr
        T_surf_ = sp.symbols('T_surf_')  # Oberflächentemperatur als Variable definieren

        # Q_latent
        if lat is True:
            Q_lat_red = Q_lat
        else:
            Q_lat_red = 0

        # Q_sensibel
        if sen is True:
            Q_sen_red = rho_w * S_w * (c_p_s * (T_Schm - T_inf) + c_p_w * (T_surf_ - T_Schm)) * (3.6e6)**-1 * A_he
        else:
            Q_sen_red = 0

        # Q_Konvektion
        if kon is True:
            Q_kon_red = alpha_kon_Bentz(u_inf) * (T_surf_ - T_inf) * A_he
        else:
            Q_kon_red = 0

        # Q_Strahlung
        if sstr is True:  # T_surf_0 statt T_surf_, da sonst 4 Nullstellen
            Q_str_red = sigma * epsilon_surf('Beton') * ((T_surf_0 + 273.15) ** 4 - T_MS(S_w, T_inf, B, Phi) ** 4) * A_he
        else:
            Q_str_red = 0

        # Q_Verdunstung
        if eva is True:
            if T_surf_0 > 0:  # bei Oberflächentemp. <= 0 °C keine Verdunstung
                Q_eva_red = rho_l * beta_c(T_inf, u_inf, h_NHN) * (X_D_sat_surf(T_surf_0, h_NHN) - X_D_inf(T_inf, Phi, h_NHN)) * h_Ph_lg * A_he
            else:
                Q_eva_red = 0
        else:
            Q_eva_red = 0

        F_Q_red = sp.Eq(Q_lat_red + Q_sen_red + R_f * (Q_kon_red + Q_str_red + Q_eva_red), 0)
        
        T_surf_sol = float(np.array(sp.solve(F_Q_red, T_surf_)))
    
    # 4.) return an _main
    if Q_sol < 0:  # Q. < 0 bei Gravitationswärmerohren nicht möglich
        Q_sol = 0

    return Q_sol, net_neg, T_surf_sol
