# -*- coding: utf-8 -*-
""" Ermittlung der Systemleistung anhand einer stationären Leistungsbilanz an der Oberfläche
    des Heizelements für jeden Zeitschritt

    Q. = (Theta_b - Theta_surf) / R_th_tot 
       = Q._lat + Q._sen + R_f * (Q._con + Q._rad + Q._eva)
    
    Leistungsbilanzen: (ohne Betrachtung der Verluste)
        - F_Q = Q._lat + Q._sen + R_f * (Q._con + Q._rad + Q._eva) - Q.
        - F_T = Q._lat + Q._sen + R_f * (Q._con + Q._rad + Q._eva)

    Fallunterscheidung:
        - Theta_b >= Theta_surf: Lösung der vollen Leistungsbilanz 
            F_Q = 0, Q. >= 0 (positiver Wärmeentzug aus dem Boden)

            {Erdboden + Oberfläche + Umgebung} -> Auflösung nach Q.

        - Theta_b < Theta_surf: Lösung der reduzierten Leistungsbilanz 
            F_T = 0, Q. := 0 (kein Wärmeentzug aus dem Boden)

            {Oberfläche + Umgebung} -> Auflösung nach Theta_surf (Oberflächentemperatur)

    Definition der Einzellasten:

        lat - latent
        sen - sensibel
        con - konvektiv
        rad - Strahlung
        eva - Verdunstung

    &

    Lösung der Leistungsbilanz - Verfahren: iterative Nullstellensuche

    Legende:
        - Temperaturen:
            - T in Kelvin [K] - für (kalorische) Gleichungen
            - Theta in Grad Celsius [°C] - Input aus dem Wetterdatenfile

    Algorithmus basierend in Teilen auf [Konrad2009]
    
    Anmerkungen zu Variablen:
        - calc_T_surf:
            - "False": Leistungsentzug aus dem Boden Q >= 0 (positiv), Oberflächentemp. Theta_surf wird in "main.py" ermittelt
            (Simulationsmodi 2 und 4)
            - "True": kein Leistungsentzug aus dem Boden, Oberflächentemp. Theta_surf wird in "load_generator.py" ermittelt
        - sb_active:
            - "True": Schneeschichtbilanzierung aktiv (es kann sich eine Schneeschicht bilden)
            - "False": Schneeschichtbilanzierung inaktiv (die Schneelast wird in jedem Zeitschritt abgeschmolzen, keine Bildung einer Schneedecke)

    Autor: Yannick Apfel
"""
import math
from scipy.constants import sigma

# Import der physikalischen Modellgleichungen
from load_generator_utils import *
from heating_element_utils import *


# Q_Konvektion = fkt(Q.) - Input für Leistungsbilanz F_Q = 0
def Q_con_Q(Q, con, u_inf, Theta_b_0, R_th, Theta_inf, A_he):  # [W]
    Q_con = 0
    if con:
        Q_con = alpha_kon_Bentz(u_inf) * (Theta_b_0 - Q * R_th - Theta_inf) * A_he

    return Q_con


# Q_Konvektion = fkt(Theta_surf) - Input für vereinfachte Leistungsbilanz F_T = 0
def Q_con_T(Theta_surf, con, u_inf, Theta_inf, A_he):  # [W]
    Q_con = 0
    if con:
        Q_con = alpha_kon_Bentz(u_inf) * (Theta_surf - Theta_inf) * A_he

    return Q_con


# Q_Strahlung = fkt(Q.) - Input für Leistungsbilanz F_Q = 0
def Q_rad_Q(Q, rad, Theta_b_0, R_th, S_w, Theta_inf, B, Phi, A_he):  # [W]
    Q_rad = 0
    if rad:
        Q_rad = sigma * epsilon_surf('Beton') * ((Theta_b_0 - Q * R_th + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

    return Q_rad


# Q_Strahlung = fkt(Theta_surf) - Input für vereinfachte Leistungsbilanz F_T = 0
def Q_rad_T(Theta_surf, rad, S_w, Theta_inf, B, Phi, A_he):  # [W]
    Q_rad = 0
    if rad:
        Q_rad = sigma * epsilon_surf('Beton') * ((Theta_surf + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

    return Q_rad


# Q_Verdunstung = fkt(Q.) - Input für Leistungsbilanz F_Q = 0
def Q_eva_Q(Q, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Theta_b_0, R_th, Phi, A_he):  # [W]
    ''' Voraussetzungen:
        - Theta_surf >= 0 °C
        - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
    '''
    Q_eva = 0
    if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
        Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_b_0 - Q * R_th, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he

    if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv!
        Q_eva = 0

    return Q_eva


# Q_Verdunstung = fkt(Theta_surf) - Input für vereinfachte Leistungsbilanz F_T = 0
def Q_eva_T(Theta_surf, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Phi, A_he):  # [W]
    ''' Voraussetzungen:
        - Theta_surf >= 0 °C
        - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
    '''
    Q_eva = 0
    if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
        Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he

    if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv!
        Q_eva = 0

    return Q_eva


# Q_sensibel = fkt(Q.) - Input für Leistungsbilanz F_Q = 0
def Q_sen_Q(Q, sen, S_w, Theta_inf, Theta_b_0, R_th, A_he):  # [W]
    Q_sen = 0
    if sen:
        Q_sen = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_b_0 - Q * R_th - Theta_Schm)) * (3.6e6)**-1 * A_he

    return Q_sen


# Q_sensibel = fkt(Theta_surf) - Input für vereinfachte Leistungsbilanz F_T = 0
def Q_sen_T(Theta_surf, sen, S_w, Theta_inf, A_he):  # [W]
    Q_sen = 0
    if sen:
        Q_sen = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_surf - Theta_Schm)) * (3.6e6)**-1 * A_he

    return Q_sen


# Q_latent - identisch für beide Leistungsbilanzen F_Q = 0 sowie F_T = 0
def Q_lat(lat, S_w, A_he):  # [W]
    Q_lat = 0
    if lat:
        Q_lat = rho_w * S_w * h_Ph_sl * (3.6e6)**-1 * A_he

    return Q_lat


""" Thermische Verluste (exkl. Oberfläche des Heizelements) Q_V = Q_V_An + Q_V_He

    Q_V_An: thermische Verluste der Anbindung zwischen Erdwärmesonden und Heizelement
        - Modellierung mittels Péclet-Gleichung für Zylinderschalen (Rohr + Isolierung) mit der Gesamtlänge aller Heatpipes
          im Bereich der Anbindung
        - Ermittlung der Rohrinnentemperatur über die Bohrlochtemperatur Theta_b und den therm. Widerstand R_th_g_hp
        - Annahme Theta_inf (Umgebungstemperatur) als Rohraußentemperatur

    Q_V_He: thermische Verluste an der Unterseite des Heizelements
        - "R_th_he_u": Verrohrung (Heatpipe-Innenseite) bis Heizelement-Unterseite (ohne Isolierung)
        - "R_th_he_iso": Isolationsschicht an Unterseite des Heizelements
"""
def Q_V(Theta_R, Theta_inf, lambda_p, lambda_iso, l_R_An, r_iso, r_pa, r_pi, he):
    
    def Q_V_An(Theta_R, Theta_inf, lambda_p, lambda_iso, l_R_An, r_iso, r_pa, r_pi):  # [W]

        Q_V_an = (Theta_R - Theta_inf) * (2 * math.pi * l_R_An) \
                * (math.log(r_pa / r_pi) / lambda_p + math.log(r_iso / r_pa) / lambda_iso) ** -1
        if Q_V_an < 0:  # Q. < 0 bei Gravitationswärmerohren nicht möglich
            Q_V_an = 0
            
        return Q_V_an
    
    def Q_V_He(he, lambda_iso, Theta_R, Theta_inf):  # [W]
        
        # therm. Widerstand Innenseite Verrohrung bis Unterseite Heizelement (ohne Isolierung)
        R_th_he_u = (he.l_R * q_l(he.D - (he.x_min + 0.5 * he.d_R_a), (he.x_min + 0.5 * he.d_R_a), he.d_R_a, he.d_R_i, he.lambda_B, he.lambda_R, he.s_R, 1, 0, state_u_insul=True)) ** -1
        
        # thermischer Widerstand der Isolierung
        R_th_he_iso = 1 / lambda_iso * he.D_iso / he.A_he
        
        Q_V_he = (Theta_R - Theta_inf) * (R_th_he_u + R_th_he_iso) ** -1
        if Q_V_he < 0:  # Q. < 0 bei Gravitationswärmerohren nicht möglich
            Q_V_he = 0
    
        return Q_V_he
    
    Q_V = Q_V_An(Theta_R, Theta_inf, lambda_p, lambda_iso, l_R_An, r_iso, r_pa, r_pi) + Q_V_He(he, lambda_iso, Theta_R, Theta_inf)
    
    return Q_V


# Leistungsbilanz F_Q {Erdboden, Oberfläche, Umgebung}, Fall Q >= 0
def F_Q(R_f, lat, S_w, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he):

    F_Q = Q_lat(lat, S_w, A_he) \
        + Q_sen_Q(Q, sen, S_w, Theta_inf, Theta_b_0, R_th, A_he) \
        + R_f \
        * (Q_con_Q(Q, con, u_inf, Theta_b_0, R_th, Theta_inf, A_he)
        + Q_rad_Q(Q, rad, Theta_b_0, R_th, S_w, Theta_inf, B, Phi, A_he)
        + Q_eva_Q(Q, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Theta_b_0, R_th, Phi, A_he)) \
        - Q

    return F_Q


# Leistungsbilanz F_T {Oberfläche, Umgebung}, Fall Q < 0
def F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he):
    F_T = Q_lat(lat, S_w, A_he) \
        + Q_sen_T(Theta_surf, sen, S_w, Theta_inf, A_he) \
        + R_f \
        * (Q_con_T(Theta_surf, con, u_inf, Theta_inf, A_he)
        + Q_rad_T(Theta_surf, rad, S_w, Theta_inf, B, Phi, A_he)
        + Q_eva_T(Theta_surf, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Phi, A_he))

    return F_T


# Solver für die Lösung der Leistungsbilanz F_Q := 0 nach Q
def solve_F_Q(R_f, con, rad, eva, sen, lat, S_w, Theta_inf, Theta_b_0, R_th, u_inf, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he):
    step_refine = 0  # Hilfsvariable zur Verfeinerung der Schrittweite
    step = 30  # doppelter Startwert als Iterationsschrittweite für Q
    res = 0.001  # zulässiges Residuum für F_Q (Restfehler)

    Q = 0  # Startwert für Q

    error = abs(F_Q(R_f, lat, S_w, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he))

    # Ermittlung von Q für F_Q = 0 (stationäres System)
    while error > res:
        if F_Q(R_f, lat, S_w, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) > 0:
            step_refine += 1
            while F_Q(R_f, lat, S_w, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) > 0:
                Q += (step / (2 * step_refine))  # Halbierung der Schrittweite für eine weitere Überschreitung/Unter- des Zielwerts
        elif F_Q(R_f, lat, S_w, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) < 0:
            step_refine += 1
            while F_Q(R_f, lat, S_w, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) < 0:
                Q -= (step / (2 * step_refine))  # Halbierung der Schrittweite für eine weitere Überschreitung/Unter- des Zielwerts

        error = abs(F_Q(R_f, lat, S_w, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he))

    # Einsetzen des ermittelten Q in die für Schneeschmelze und Verdunstung relevanten Terme
    Q_lat_sol = Q_lat(lat, S_w, A_he)
    Q_sen_sol = Q_sen_Q(Q, sen, S_w, Theta_inf, Theta_b_0, R_th, A_he)
    Q_eva_sol = Q_eva_Q(Q, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Theta_b_0, R_th, Phi, A_he)
    Q_load = Q

    return Q_load, Q_lat_sol, Q_sen_sol, Q_eva_sol


# Solver für die Lösung der Leistungsbilanz F_T := 0 nach Theta_surf
def solve_F_T(R_f, con, rad, eva, sen, lat, S_w, Theta_inf, u_inf, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he):
    step_refine = 0  # Hilfsvariable zur Verfeinerung der Schrittweite
    step = 10  # doppelter Startwert als Iterationsschrittweite für T
    res = 0.01  # zulässiges Residuum für F_T (Restfehler)

    Theta_surf = 0  # Startwert für Q

    error = abs(F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he))

    # Ermittlung von Theta_surf für F_T = 0 (stationäres System)
    while error > res:
        if F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) > 0:
            step_refine += 1
            while F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) > 0:
                Theta_surf -= (step / (2 * step_refine))  # Halbierung der Schrittweite für eine weitere Überschreitung/Unter- des Zielwerts
        elif F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) < 0:
            step_refine += 1
            while F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) < 0:
                Theta_surf += (step / (2 * step_refine))  # Halbierung der Schrittweite für eine weitere Überschreitung/Unter- des Zielwerts

        error = abs(F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he))

    # Einsetzen des ermittelten Theta_surf in die für Schneeschmelze und Verdunstung relevanten Terme
    Q_lat_sol = Q_lat(lat, S_w, A_he)
    Q_sen_sol = Q_sen_T(Theta_surf, sen, S_w, Theta_inf, A_he)
    Q_eva_sol = Q_eva_T(Theta_surf, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Phi, A_he)
    Theta_surf_sol = Theta_surf

    return Theta_surf_sol, Q_lat_sol, Q_sen_sol, Q_eva_sol


def load(h_NHN, v, Theta_inf, S_w, he, Theta_b_0, R_th, R_th_ghp, Theta_surf_0, B, Phi, RR, m_Rw_0, m_Rs_0, start_sb, 
         l_R_An, lambda_p, lambda_iso, r_iso, r_pa, r_pi):
    ''' Hauptfunktion zur Ermittlung der Entzugsleistung pro Zeitschritt, Aufteilung in:
        - Q_load: aus Oberflächenbilanzen ermittelte thermische Leistung der Oberfläche
        - Q_N: Anteil von Q_load, der tatsächlich zur Schneeschmelze benutzt wird - also Q_load abzügl. Oberflächenverlusten Q_con, Q_rad, Q_eva
        - Q_V: Verlustleistungen (gehen nicht in Leistungsbilanzen mit ein)
            - Anbindung (An) zwischen Bohrloch und Heizelement (Heatpipes)
            - Unterseite Heizelement (He)
    '''
    # Theta_x_0: Temp. des vorhergehenden Zeitschritts

    # 0.) Preprocessing
    u_inf = u_eff(v)  # Reduzierte Windgeschwindigkeit (logarithmisches Windprofil)

    # Hilfsvariablen
    calc_T_surf = False  # "True" falls Oberlächentemp. bereits in diesem Modul ermittelt wird
    Theta_surf_sol = None

    ''' Simulationsmodi 1-3: Schnee wird verzögert abgeschmolzen (Bildung einer Schneedecke ist möglich)
        Simulationsmodi 4-5: Schnee wird instantan (=innerhalb des Zeitschritts) abgeschmolzen
    '''
    # Simulationsmodus ermitteln:
    if (m_Rs_0 > 0 or start_sb is True):  # Schneedecke bildet sich
        sb_active = 1  # snow-balancing aktivieren
    else:  # Oberfläche ist schnee-frei
        sb_active = 0

    # 1.) Teil-Wärmeströme
    con = True  # aktivieren oder deaktivieren (für unit-testing)
    rad = True
    eva = True
    sen = True
    lat = True

    # 2.) Ermittlung Entzugsleistung Q_load und Oberflächentemperatur Theta_surf_sol
    ''' Simulationsmodi 1-3'''
    if (sb_active == 1):  # "Schnee wird verzögert abgeschmolzen" (Bildung einer Schneedecke)

        ''' Erdboden, Heizelement-Oberfläche und Umgebung bilden ein stationäres System, wobei
            sich eine Schneedecke bildet, sodass T_surf = T_Schm (Schmelzwasser).
            - Q._0 = (T_b - T_surf) / R_th definiert zur Verfügung stehende Leistung (delta-T)
            - Q._R = Q._0 - (Q._con + Q._rad + Q._eva) ergibt die zur Schneeschmelze (= Q._lat + Q._sen) vor-
            handene restliche Leistung

            - Simulationsmodus 1: Q._0 < 0, Lösung von F_T nach Theta_surf
            (kein Wärmeentzug aus dem Boden möglich, da kein delta-T vorhanden - sibirische Verhältnisse)
            
            - Simulationsmodus 2: Q._R < 0, Lösung von F_Q nach Q.
            (keine Restleistung zur Schneeschmelze vorhanden, Leistung geht für Konvektions- und Strahlungsverluste drauf)
            
            - Simulationsmodus 3: Q._R > 0, setze Theta_surf := Theta_Schmelz (Oberfläche mit Wasser benetzt)
            (Q._R wird zur Schneeschmelze verwendet)
            
        '''

        # 2.1) Pre-Processing
        R_f = 0.2  # free-area ratio
        Theta_surf_0 = Theta_Schm  # Fixieren der Oberflächentemperatur

        # 2.2) verfügbare Entzugsleistung
        Q_0 = (Theta_b_0 - Theta_surf_0) * R_th ** -1

        # 2.3) Fallunterscheidung Entzugsleistung Q_0
        ''' Simulationsmodus 1'''
        if Q_0 < 0:  # keine nutzbare T-Differenz im Boden vorhanden (sibirische Verhältnisse)

            sim_mod = 1  # Simulationsmodus aufzeichnen

            calc_T_surf = True

            # Q_sensibel, Q_latent, Q_Verdunstung = 0
            sen, lat, eva = 0, 0, 0  # Energie zur Schneeschmelze kommt definitorisch aus dem Boden, nicht der Umgebung

            # 2.4) iterative Lösung der stationären Leistungsbilanz F_T = 0 am Heizelement (Oberfläche + Umgebung) nach T
            Theta_surf_sol, Q_lat, Q_sen, Q_eva = solve_F_T(R_f, con, rad, eva, sen, lat, S_w, Theta_inf, u_inf, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, he.A_he)

            Q_load = -1  # keine Entzugsleistung aus dem Boden

        else:  # Q_0 >= 0, nutzbare T-Differenz im Boden vorhanden (Regelfall)
            ''' Simulationsmodi 2 & 3'''
            # 2.4) Oberflächenverluste (explizit für Theta_surf_0 = Theta_Schm formuliert)

            # Q_Konvektion
            Q_con = Q_con_T(Theta_surf_0, con, u_inf, Theta_inf, he.A_he)

            # Q_Strahlung
            Q_rad = Q_rad_T(Theta_surf_0, rad, S_w, Theta_inf, B, Phi, he.A_he)

            # Q_Verdunstung
            Q_eva = 0

            # 2.5) Ermittlung vorhandene Restleistung (f. Schnee- und Eisschmelze)
            Q_R = Q_0 - R_f * (Q_con + Q_rad + Q_eva)

            # 2.6) Fallunterscheidung Restleistung Q_R
            ''' Simulationsmodus 2'''
            if Q_R < 0:  # keine Restleistung für Schneeschmelze vorhanden

                sim_mod = 2  # Simulationsmodus aufzeichnen

                # Q_sensibel, Q_latent, Q_Verdunstung = 0
                sen, lat, eva = 0, 0, 0  # Energie zur Schneeschmelze kommt definitorisch aus dem Boden, nicht der Umgebung

                # 2.7) iterative Lösung der stationären Leistungsbilanz F_Q = 0 am Heizelement (Erdboden + Oberfläche + Umgebung))
                Q_load, Q_lat, Q_sen, Q_eva = solve_F_Q(R_f, con, rad, eva, sen, lat, S_w, Theta_inf, Theta_b_0, R_th, u_inf, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, he.A_he)

            else:  # Q_R >= 0, Restleistung für Schneeschmelze vorhanden
                ''' Simulationsmodus 3'''
                sim_mod = 3  # Simulationsmodus aufzeichnen

                calc_T_surf = True

                # 2.7) Volumenstrom der Schneeschmelze
                V_s = Q_R / (rho_w * (h_Ph_sl + c_p_s * (Theta_Schm - Theta_inf)))
                if V_s < 0:  # Schmelz-Volumenstrom ist definitorisch positiv
                    V_s = 0

                # 2.8) Schmelzterme (explizit)

                # Q_sensibel
                Q_sen = 0
                if sen:
                    Q_sen = rho_w * c_p_s * (Theta_Schm - Theta_inf) * V_s

                # Q_latent
                Q_lat = 0
                if lat:
                    Q_lat = rho_w * h_Ph_sl * V_s

                Theta_surf_sol = Theta_Schm
                Q_load = Q_0

                ''' Simulationsmodi 4 & 5'''
                ''' Erdboden, Heizelement-Oberfläche und Umgebung bilden ein stationäres System, wobei sich
                        die Oberflächentemp. Theta_surf entsprechend der Oberflächenlasten ergibt.
                        Die Leistungsbilanz F_Q = Q._lat + Q._sen + R_f(Q._con + Q._rad + Q._eva) - Q. wird für
                        - Fall Q. >= 0 (Leistungsentzug aus dem Boden) nach Q. aufgelöst - Simulationsmodus 4
                        - Fall Q. < 0 : Lösung der vereinfachten Leistungsbilanz F_T = F_Q(Q.=0) nach Theta_surf
                            (kein Leistungsentzug aus dem Boden, Umgebung erwärmt Heizelement) - Simulationsmodus 5
                '''
    else:  # Schnee wird instantan abgeschmolzen (keine Bildung einer Schneedecke)
        ''' Simulationsmodus 4'''
        sim_mod = 4  # Simulationsmodus aufzeichnen

        # 2.1) Pre-Processing
        R_f = 1  # free-area ratio

        # 2.2) iterative Lösung der stationären Leistungsbilanz F_Q = 0 am Heizelement (Erdboden + Oberfläche + Umgebung) nach Q.
        Q_load, Q_lat, Q_sen, Q_eva = solve_F_Q(R_f, con, rad, eva, sen, lat, S_w, Theta_inf, Theta_b_0, R_th, u_inf, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, he.A_he)

        # 2.3) Fall Q. < 0
        ''' Simulationsmodus 5'''
        if Q_load < 0:  # kein Wärmeentzug aus Erdboden

            sim_mod = 5  # Simulationsmodus aufzeichnen

            calc_T_surf = True

            # Q_sensibel, Q_latent = 0
            sen, lat = 0, 0  # Energie zur Schneeschmelze kommt definitorisch aus dem Boden, nicht der Umgebung

            # 2.4) iterative Lösung der stationären Leistungsbilanz F_T = 0 am Heizelement (Oberfläche + Umgebung) nach T
            Theta_surf_sol, Q_lat, Q_sen, Q_eva = solve_F_T(R_f, con, rad, eva, sen, lat, S_w, Theta_inf, u_inf, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, he.A_he)

    # 3.) Wasser- und Schneehöhenbilanz auf Heizelement

    # Ermittlung Restwassermenge
    m_w_1 = m_Restwasser(m_Rw_0, RR, he.A_he, Q_eva)

    # Ermittlung Restschneemenge
    m_s_1 = m_Restschnee(m_Rs_0, S_w, he.A_he, Q_lat, sb_active)

    # 4.) Auswertung der Entzugsleistung Q_load, Nutzleistung Q_N und Verlustleistung Q_V (Anbindung und Unterseite Heizelement) [W]

    # 4.1) Q_load [W]
    if Q_load < 0:  # Q. < 0 bei Gravitationswärmerohren nicht möglich
        Q_load = 0

    # 4.2) Q_N [W]
    Q_N = Q_lat + Q_sen

    # 4.3) Q_V [W]
    Q_V_sol = Q_V(Theta_b_0 - Q_load * R_th_ghp, Theta_inf, lambda_p, lambda_iso, l_R_An, r_iso, r_pa, r_pi, he)

    return Q_load, Q_N, Q_V_sol, calc_T_surf, Theta_surf_sol, m_w_1, m_s_1, sb_active, sim_mod
