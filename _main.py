# -*- coding: utf-8 -*-
""" GERDPy - Main-File
    Steuerungsfile des Auslegungstools für GERDI

    Autor: Yannick Apfel

    todos:
        - GUI für End-User
        - Kontrollstruktur für Wetterdaten (Errorhandling beim Einlesen)
        
    Legende:
        - Temperaturen:
            - T in Kelvin [K] - für (kalorische) Gleichungen
            - Theta in Grad Celsius [°C] - Input aus dem Wetterdatenfile
"""
import sys
import matplotlib.pyplot as plt
import time as tim
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from scipy.constants import pi

# import GERDPy-modules
import boreholes, heatpipes, heating_element, gfunction, load_aggregation
from load_generator import *
from weather_data import get_weather_data
from geometrycheck import check_geometry
from R_th_tot import R_th_tot


def main():
    # -------------------------------------------------------------------------
    # 1.) Parametrierung der Simulation (Geometrien, Stoffwerte, etc.)
    # -------------------------------------------------------------------------

    # 1.0) Standort
    h_NHN = 520                # Höhe über Normal-Null des Standorts

    # 1.1) Erdboden
    a = 1.0e-6                  # Temperaturleitfähigkeit [m2/s]
    lambda_g = 2.0              # Wärmeleitfähigkeit [W/mK]
    Theta_g = 10.0                  # ungestörte Bodentemperatur [°C]

    # 1.2) Erdwärmesondenfeld

    # Geometrie-Import (.txt-Datei)
    boreField = boreholes.field_from_file('./data/custom_field_5.txt')

    # Sondenmeter (gesamt)
    H_field = boreholes.length_field(boreField)

    # Layout-Plot des Erdwärmesondenfelds
    # boreholes.visualize_field(boreField)

    # 1.3) Bohrloch

    # Geometrie
    N = 6                        # Anzahl Heatpipes pro Bohrloch [-]
    r_b = boreField[0].r_b  # Radius der Erdwärmesondenbohrung [m]
    r_w = 0.12                # Radius der Wärmerohr-Mittelpunkte [m]
    r_pa = 0.016                 # Außenradius der Isolationsschicht [m]
    r_iso = 0.016                # Innenradius der Isolationsschicht [m]
    r_pi = 0.015                # Innenradius des Wärmerohrs [m]

    # Wärmeleitfähigkeiten
    lambda_b = 2                # lambda der Hinterfüllung [W/mK]
    lambda_iso = 0.3            # lambda der Isolationsschicht [W/mK]
    lambda_p = 14               # lambda der Heatpipe [W/mK]

    # Geometrie-Erstellung
    hp = heatpipes.Heatpipes(N, r_b, r_w, r_pa, r_iso, r_pi, lambda_b,
                             lambda_iso, lambda_p)
    # Layout-Plot der Wärmerohrkonfiguration
    # hp.visualize_hp_config()

    # 1.4) Heizelement

    # Fläche Heizelement [m2]
    A_he = 50

    # minimaler Oberflächenabstand [mm]
    x_min = 15

    # Wärmeleitfähigkeit des Betonelements [W/mk]
    lambda_Bet = 2.3

    # Geometrie-Erstellung
    he = heating_element.HeatingElement(A_he, x_min, lambda_Bet)

    # 2.) Simulation

    # Simulationsparameter
    dt = 3600.                           # Zeitschrittweite [s]
    tmax = 1 * 12 * (8760./12) * 3600.    # Gesamt-Simulationsdauer [s]
    Nt = int(np.ceil(tmax/dt))           # Anzahl Zeitschritte [-]

    # -------------------------------------------------------------------------
    # 2.) Überprüfung der geometrischen Verträglichkeit (Sonden & Heatpipes)
    # -------------------------------------------------------------------------

    if check_geometry(boreField, hp):
        print(50*'-')
        print('Geometry-Check: not OK! - Simulation aborted')
        print(50*'-')
        sys.exit()
    else:
        print(50*'-')
        print('Geometry-Check: OK!')
        print(50*'-')

    # -------------------------------------------------------------------------
    # 3.) Ermittlung thermischer Widerstand Bohrlochrand bis Oberfläche
    # -------------------------------------------------------------------------

    R_th = R_th_tot(lambda_g, boreField, hp, he)

    # -------------------------------------------------------------------------
    # 4.) Ermittlung der G-Function (Bodenmodell)
    # -------------------------------------------------------------------------

    # Initialisierung Zeitstempel (Simulationsdauer)
    tic = tim.time()

    # Aufsetzen der Simulationsumgebung mit 'load_aggregation.py'
    LoadAgg = load_aggregation.ClaessonJaved(dt, tmax)
    time_req = LoadAgg.get_times_for_simulation()

    # Berechnung der G-Function mit 'gfunction.py'
    gFunc = gfunction.uniform_temperature(boreField, time_req, a,
                                          nSegments=12)

    # Initialisierung der Simulation mit 'load_aggregation.py'
    LoadAgg.initialize(gFunc/(2*pi*lambda_g))

    # -------------------------------------------------------------------------
    # 5.) Import / Generierung der Wetterdaten
    # -------------------------------------------------------------------------

    # Import Wetterdaten aus weather_data.py
    u_inf, Theta_inf, S_w, B, Phi, RR = get_weather_data(Nt)
    ''' u_inf - Windgeschwindigkeit [m/s]
        Theta_inf - Umgebungstemperatur [°C]
        S_w - Schneefallrate (Wasserequivalent) [mm/s]
        B - Bewölkungsgrad [octal units / 8]
        Phi - rel. Luftfeuchte [%]
    '''

    # -------------------------------------------------------------------------
    # 6.) Iterationsschleife (Simulation mit Nt Zeitschritten der Länge dt)
    # -------------------------------------------------------------------------

    time = 0.
    i = -1
    start_sb = False  # sb - snow balancing

    # Initialisierung Temperaturen
    ''' Vektoren:
            Q[i] - Entzugsleistung f. Zeitschritt i 
            Theta_surf[i] - Oberflächentemp. f. Zeitschritt i
    '''
    Theta_b = np.zeros(Nt)      # Bohrlochrand
    Theta_surf = np.zeros(Nt)   # Oberfläche Heizelement
    
    # Initialisierung Vektor für Restwassermenge [mm]
    m_Rw = np.zeros(Nt)
    
    # Initialisierung Vektor für Restschneemenge [mm H2O]
    m_Rs = np.zeros(Nt)

    # Initialisierung Entnahmeleistung [W]
    Q = np.zeros(Nt)

    print('Simulating...')

    # Hilfsgrößen
    start_sb_counter = np.zeros(Nt)
    sb_active = np.zeros(Nt)
    sim_mod = np.zeros(Nt)

    while time < tmax:  # Iterationsschleife (ein Durchlauf pro Zeitschritt)
        
        # Zeitschritt um 1 inkrementieren
        if start_sb == False:
            time += dt
            i += 1

        LoadAgg.next_time_step(time)

        # Ermittlung der Entzugsleistung im 1. Zeitschritt
        if i == 0:  # Annahme Theta_b = Theta_surf = Theta_g, Heizelementoberfläche trocken und Schnee-frei
            Q[i], calc_T_surf, Theta_surf[i], m_Rw[i], m_Rs[i], sb_active[i], sim_mod[i] = load(h_NHN, u_inf[i], Theta_inf[i], S_w[i], A_he, Theta_g, 
                                                         R_th, Theta_g, B[i], Phi[i], RR[i], 0, 0, start_sb)

        if i > 0:  # alle weiteren Zeitschritte (ermittelte Bodentemperatur)
            Q[i], calc_T_surf, Theta_surf[i], m_Rw[i], m_Rs[i], sb_active[i], sim_mod[i] = load(h_NHN, u_inf[i], Theta_inf[i], S_w[i], A_he, Theta_b[i-1], 
                                                         R_th, Theta_surf[i-1], B[i], Phi[i], RR[i], m_Rw[i-1], m_Rs[i-1], start_sb)
            
        start_sb = False  # Variable Start-Schneebilanzierung zurücksetzen

        # Aufprägung der ermittelten Entzugsleistung mit 'load_aggregation.py'
        LoadAgg.set_current_load(Q[i]/H_field)

        # Temperatur am Bohrlochrand
        deltaTheta_b = LoadAgg.temporal_superposition()
        Theta_b[i] = Theta_g - deltaTheta_b

        # Temperatur an der Oberfläche des Heizelements
        if calc_T_surf == False:
            Theta_surf[i] = Theta_b[i] - Q[i] * R_th

        # Schneebilanzierung starten
        ''' Zeitschritt wird einmalig wiederholt, falls sich eine Schneeschicht beginnt zu bilden:
            - Theta_surf[i] < 0 UND
            - S_w[i] > 0 UND
            - m_Rs[i]==0 (keine Restschneemenge vorhanden)
            
            => Schnee bleibt liegen: Schneebilanz an Oberfläche
        '''
        if (Theta_surf[i] < 0 and S_w[i] > 0 and m_Rs[i] == 0):
            start_sb = True
            start_sb_counter[i] = 1
            
        print(f'Iter: {i}, Sim-Mode: {sim_mod[i]}')
            
    # Zeitstempel (Simulationsdauer)
    toc = tim.time()
    print('Total simulation time: {} sec'.format(toc - tic))
            
    # -------------------------------------------------------------------------
    # 7.) zeitlich gemittelte Leistung & Gesamtenergie
    # -------------------------------------------------------------------------
            
    # Leistung - über Zeitintervall gemittelt
    interv = 72  # Zeitintervall der gemittelten Leistung
    
    Q_m = np.zeros(Nt)  # [W]
    for i in range(0, Nt, interv):  # gemittelte Leistungen
        interval = [x for x in Q[i:(i+interv)]]
        Q_m[i:(i+interv)] = np.mean(interval)
        
    E_tot = (np.sum(Q) / len(Q)) * len(Q) / 1e6  # [MWh]
    print(f'Dem Boden wurden {E_tot} MWh entnommen')
    
    # -------------------------------------------------------------------------
    # 8.) Plots
    # -------------------------------------------------------------------------
    
    hours = np.array([(j+1)*dt/3600. for j in range(Nt)])
    
    # -------------------------------------------------------------------------
    # 8.1) Figure 1
    # -------------------------------------------------------------------------
    
    plt.rc('figure')
    fig1 = plt.figure()

    font = {'weight': 'bold', 'size': 22}
    plt.rc('font', **font)

    # Lastprofil (thermische Leistung Q. über die Simulationsdauer)
    ax1 = fig1.add_subplot(311)
    ax1.set_xlabel(r'$t$ [h]')
    ax1.set_ylabel(r'$q$ [W/m²]')
    ax1.plot(hours, Q / A_he, 'k-', lw=0.6)
    ax1.plot(hours, Q_m / A_he, 'r-', lw=2)
    ax1.legend(['spezifische Entzugsleistung [W/m²]'],
               prop={'size': font['size'] - 5}, loc='upper right')
    ax1.grid('major')
    
    # Wasser- und Schneebilanzlinie (entspr. Wasserhöhe an der Oberfläche)
    ax2 = fig1.add_subplot(312)
    ax2.set_ylabel('[mm]')
    ax2.plot(hours, m_Rw / A_he, 'b-', lw=0.8)
    ax2.plot(hours, m_Rs / A_he, 'g-', lw=0.8)
    ax2.legend(['Wasserhöhe', 'Schneehöhe'],
                          prop={'size': font['size'] - 5}, loc='upper left')
    ax2.grid('major')

    # Temperaturverläufe
    ax3 = fig1.add_subplot(313)
    ax3.set_ylabel(r'$T$ [°C]')
    # plots
    ax3.plot(hours, Theta_b, 'r-', lw=1.2)
    ax3.plot(hours, Theta_surf, 'c-', lw=0.6)
    ax3.legend(['T_Bohrlochrand', 'T_Oberfläche'],
               prop={'size': font['size'] - 5}, loc='upper right')
    ax3.grid('major')

    # -------------------------------------------------------------------------
    # 8.2) Figure 2
    # -------------------------------------------------------------------------
    
    plt.rc('figure')
    fig2 = plt.figure()

    font = {'weight': 'bold', 'size': 22}
    plt.rc('font', **font)

    # Darstellungen Simulationsmodus
    ax4 = fig2.add_subplot(311)
    ax4.set_xlabel(r'$t$ [h]')
    # ax4.plot(hours, start_sb_counter, 'k--', lw=1.5)
    ax4.plot(hours, sb_active, 'g-', lw=1.3)
    ax4.plot(hours, sim_mod, 'y-', lw=1.3)
    ax4.legend(['sb_active', 'sim_mod'],
               prop={'size': font['size'] - 5}, loc='upper right')
    ax4.grid('major')
    
    # Schneefallrate und -bilanzlinie
    ax5 = fig2.add_subplot(312)
    ax5_2 = ax5.twinx()
    ax5.set_ylabel('Schneefallrate [mm/h]')
    ax5_2.set_ylabel('Schneehöhe [mm]')
    ax5.plot(hours, S_w, 'b-', lw=0.8)
    ax5_2.plot(hours, m_Rs / A_he, 'g-', lw=0.8)
    ax5.legend(['Schneefallrate'],
                          prop={'size': font['size'] - 5}, loc='upper left')
    ax5_2.legend(['Schneehöhe'],
                          prop={'size': font['size'] - 5}, loc='upper right')
    ax5.grid('major')

    # Temperaturverläufe
    ax6 = fig2.add_subplot(313)
    ax6.set_ylabel(r'$T$ [°C]')
    # plots
    ax6.plot(hours, Theta_b, 'r-', lw=1.2)
    ax6.plot(hours, Theta_surf, 'c-', lw=0.6)
    ax6.legend(['T_Bohrlochrand', 'T_Oberfläche'],
               prop={'size': font['size'] - 5}, loc='upper right')
    ax6.grid('major')

    # Beschriftung Achsenwerte
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    ax4.xaxis.set_minor_locator(AutoMinorLocator())
    ax4.yaxis.set_minor_locator(AutoMinorLocator())
    ax5.xaxis.set_minor_locator(AutoMinorLocator())
    ax5.yaxis.set_minor_locator(AutoMinorLocator())
    ax6.xaxis.set_minor_locator(AutoMinorLocator())
    ax6.yaxis.set_minor_locator(AutoMinorLocator())
    # plt.tight_layout()  # Fenstergröße anpassen

    return


# Main function
if __name__ == '__main__':
    main()
