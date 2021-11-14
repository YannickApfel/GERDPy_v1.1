# -*- coding: utf-8 -*-
""" Ermittlung der Systemleistung anhand einer Leistungsbilanz an der Oberfläche
    des Heizelements für jeden Zeitschritt

    Q. = (T_g - Theta_surf) / R_th_tot = Summe Wärmeströme am Heizelement

    Legende:
        - Temperaturen:
            - T in Kelvin [K] - für (kalorische) Gleichungen
            - Theta in Grad Celsius [°C] - Input aus dem Wetterdatenfile

    Autor: Yannick Apfel
"""
import sympy as sp
import numpy as np
from load_generator_utils import *


# Definition & Bilanzierung der Einzellasten
def load(h_NHN, v, Theta_inf, S_w, A_he, Theta_b_0, R_th, Theta_surf_0, B, Phi, RR, m_Rw_0, m_Rs_0, ssb):  
    # Theta_x_0: Temp. des vorhergehenden Zeitschritts

    # 0.) Preprocessing
    u_inf = u_eff(v)  # Reduzierte Windgeschwindigkeit (logarithmisches Windprofil)

    # Hilfsvariablen
    calc_T_surf = False  # "True" falls Oberlächentemp. bereits in diesem Modul ermittelt wird
    Theta_surf_sol = None

    ''' Simulationsmodi:
        - "Schnee wird verzögert abgeschmolzen" (Bildung einer Schneedecke)
        - Schnee wird instantan (=innerhalb des Zeitschritts) abgeschmolzen
    '''
    # Simulationsmodus ermitteln:
    if (m_Rs_0 > 0 or ssb is True):  # Schneedecke bildet sich
        sb_active = 1  # snow-balancing aktivieren
    else:  # Oberfläche ist schnee-frei
        sb_active = 0

    # 1.) Teil-Wärmeströme
    con = True  # aktivieren oder deaktivieren (für unit-testing)
    rad = True
    eva = True
    sen = True
    lat = True

    # 2.) Ermittlung Entzugsleistung Q_sol und Oberflächentemperatur T_surf_sol
    if (sb_active == 1):  # "Schnee wird verzögert abgeschmolzen" (Bildung einer Schneedecke)

        ''' Erdboden, Heizelement-Oberfläche und Umgebung bilden ein stationäres System, wobei
            sich eine Schneedecke bildet, sodass T_surf = T_Schm (Schmelzwasser).
            - Q._0 = (T_b - T_surf) / R_th definiert zur Verfügung stehende Leistung (delta-T)
            - Q._R = Q._0 - (Q._con + Q._rad + Q._eva) ergibt die zur Schneeschmelze (= Q._lat + Q._sen) vor-
            handene restliche Leistung

            Die Leistungsbilanz Q. = Q._lat + Q._sen + R_f(Q._con + Q._rad + Q._eva) wird nach
             - Q. -> Fall Q. >= 0 (Leistungsentzug aus dem Boden) oder
             - Theta_surf -> Fall Q. < 0 (kein Leistungsentzug aus dem Boden)
            aufgelöst, falls Q._R < 0 (keine Restleistung f. Schneeschmelze vorhanden) bzw.
            Q._0 < 0 (kein Wärmeentzug aus dem Boden möglich - kein delta-T vorhanden)
        '''

        # 2.1) Pre-Processing
        R_f = 0.04  # free-area ratio
        Theta_surf_0 = Theta_Schm  # Fixieren der Oberflächentemperatur

        # 2.2) verfügbare Entzugsleistung
        Q_0 = (Theta_b_0 - Theta_surf_0) * R_th ** -1

        # 2.3) Fallunterscheidung Entzugsleistung Q_0
        if Q_0 < 0:  # keine nutzbare T-Differenz im Boden vorhanden (sibirische Verhältnisse)

            calc_T_surf = True

            # 2.4) Verlustleistungsterme (parametriert) - Parameter Theta_surf_
            Theta_surf_ = sp.symbols('Theta_surf_')  # Oberflächentemperatur als Parameter definieren

            # Q_Konvektion
            Q_con = 0
            if con:
                Q_con = alpha_kon_Bentz(u_inf) * (Theta_surf_ - Theta_inf) * A_he

            # Q_Strahlung
            Q_rad = 0
            if rad:  # Theta_surf_0 statt Theta_surf_, da sonst 4 Nullstellen
                Q_rad = sigma * epsilon_surf('Beton') * ((Theta_surf_0 + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

            # Q_Verdunstung
            ''' Voraussetzungen: 
                    - Theta_surf >= 0 °C
                    - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
            '''
            Q_eva = 0
            if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
                Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf_0, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he
            if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
                Q_eva = 0  

            # Q_sen, Q_lat = 0
            Q_sen, Q_lat = 0, 0

            # stationäre Leistungbilanz (Oberfläche + Umgebung) & Auflösung nach Theta_surf_
            F_T = sp.Eq(Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva), 0)
            Theta_surf_sol = float(np.array(sp.solve(F_T, Theta_surf_)))

            Q_sol = -1  # keine Entzugsleistung aus dem Boden

        else:  # Q_0 >= 0 nutzbare T-Differenz im Boden vorhanden (Regelfall)

            # 2.4) Verlustleistungsterme (explizit)

            # Q_Konvektion
            Q_con = 0
            if con:
                Q_con = alpha_kon_Bentz(u_inf) * (Theta_surf_0 - Theta_inf) * A_he

            # Q_Strahlung
            Q_rad = 0
            if rad:
                Q_rad = sigma * epsilon_surf('Beton') * ((Theta_surf_0 + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

            # Q_Verdunstung
            ''' Voraussetzungen: 
                    - Theta_surf >= 0 °C
                    - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
            '''
            Q_eva = 0
            if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
                Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf_0, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he

            if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
                Q_eva = 0

            # 2.5) Ermittlung vorhandene Restleistung (f. Schnee- und Eisschmelze)
            Q_R = Q_0 - R_f * (Q_con + Q_rad + Q_eva)

            # 2.6) Fallunterscheidung Restleistung Q_R
            if Q_R < 0:  # keine Restleistung für Schneeschmelze vorhanden

                calc_T_surf = True

                # 2.7) Verlustleistungsterme (parametriert) - Parameter Q
                Q = sp.symbols('Q')  # thermische Leistung Q als Parameter definieren

                # Q_Konvektion
                Q_con = 0
                if con:
                    Q_con = alpha_kon_Bentz(u_inf) * (Theta_b_0 - Q * R_th - Theta_inf) * A_he

                # Q_Strahlung
                Q_rad = 0
                if rad:  # Theta_surf_0 statt Theta_surf_, da sonst 4 Nullstellen
                    Q_rad = sigma * epsilon_surf('Beton') * ((Theta_surf_0 + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

                # Q_eva (vererbt)

                # Q_sensibel, Q_latent = 0
                Q_sen, Q_lat = 0, 0

                # 2.8) stationäre Leistungbilanz (Erdboden + Oberfläche + Umgebung) & Auflösung nach Q
                F_Q = sp.Eq(Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva) - Q, 0)
                Q_sol = float(np.array(sp.solve(F_Q, Q)))

            else:  # Restleistung für Schneeschmelze vorhanden

                calc_T_surf = True

                # 2.7) Volumenstrom der Schneeschmelze
                V_s = Q_R / (rho_w * (h_Ph_sl + c_p_s * (Theta_Schm - Theta_inf)))

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
                Q_sol = Q_0
    
    else:  # Schnee wird instantan abgeschmolzen (keine Bildung einer Schneedecke)

        ''' Erdboden, Heizelement-Oberfläche und Umgebung bilden ein stationäres System, wobei sich
            die Oberflächentemp. T_surf entsprechend der Oberflächenlasten ergibt.
            Die Leistungsbilanz Q. = Q._lat + Q._sen + R_f(Q._con + Q._rad + Q._eva) wird nach
             - Q. -> Fall Q. >= 0 (Leistungsentzug aus dem Boden) oder
             - Theta_surf -> Fall Q. < 0 (kein Leistungsentzug aus dem Boden, Umgebung erwärmt Heizelement) 
            aufgelöst.
        '''

        # 2.1) Pre-Processing
        R_f = 1  # free-area ratio

        # 2.2) Schmelz- und Verlustleistungsterme (parametriert) - Parameter Q
        Q = sp.symbols('Q')  # thermische Leistung Q als Parameter definieren

        # Q_Konvektion
        Q_con = 0
        if con:
            Q_con = alpha_kon_Bentz(u_inf) * (Theta_b_0 - Q * R_th - Theta_inf) * A_he

        # Q_Strahlung
        Q_rad = 0
        if rad:  # Theta_surf_0 statt Theta_surf_, da sonst 4 Nullstellen
            Q_rad = sigma * epsilon_surf('Beton') * ((Theta_surf_0 + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

        # Q_Verdunstung
        ''' Voraussetzungen: 
                - Theta_surf >= 0 °C
                - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
        '''
        Q_eva = 0
        if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
            Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf_0, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he
        if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
            Q_eva = 0  

        # Q_sensibel
        Q_sen = 0
        if sen:
            Q_sen = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_b_0 - Q * R_th - Theta_Schm)) * (3.6e6)**-1 * A_he

        # Q_latent
        Q_lat = 0
        if lat:
            Q_lat = rho_w * S_w * h_Ph_sl * (3.6e6)**-1 * A_he

        # 2.3) stationäre Leistungbilanz (Erdboden + Oberfläche + Umgebung) & Auflösung nach Q
        F_Q = sp.Eq(Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva) - Q, 0)
        Q_sol = float(np.array(sp.solve(F_Q, Q)))

        # 2.4) Fall Q. < 0
        if Q_sol < 0:  # kein Wärmeentzug aus Erdboden
            calc_T_surf = True

            # Schmelz- und Verlustleistungsterme (parametriert) - Parameter Theta_surf_
            Theta_surf_ = sp.symbols('Theta_surf_')  # Oberflächentemperatur als Parameter definieren

            # Q_Konvektion
            Q_con = 0
            if con:
                Q_con = alpha_kon_Bentz(u_inf) * (Theta_surf_ - Theta_inf) * A_he

            # Q_Strahlung (vererbt)

            # Q_Verdunstung
            ''' Voraussetzungen: 
                    - Theta_surf >= 0 °C
                    - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
            '''
            Q_eva = 0
            if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
                Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf_0, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he
            if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
                Q_eva = 0  

            # Q_sensibel
            Q_sen = 0
            if sen:
                Q_sen= rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_surf_ - Theta_Schm)) * (3.6e6)**-1 * A_he

            # Q_latent (vererbt)

            # stationäre Leistungbilanz (Oberfläche + Umgebung) & Auflösung nach Theta_surf_
            F_T = sp.Eq(Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva), 0)
            Theta_surf_sol = float(np.array(sp.solve(F_T, Theta_surf_)))

    # 3.) Wasser- und Schneehöhenbilanz auf Heizelement

    # Ermittlung Restwassermenge
    m_Rw_1 = m_Restwasser(m_Rw_0, RR, A_he, Q_eva)

    # Ermittlung Restschneemenge
    m_Rs_1 = m_Restschnee(m_Rs_0, S_w, A_he, Q_lat, sb_active)

    # 4.) return an _main
    if Q_sol < 0:  # Q. < 0 bei Gravitationswärmerohren nicht möglich
        Q_sol = 0

    return Q_sol, calc_T_surf, Theta_surf_sol, m_Rw_1, m_Rs_1, sb_active
