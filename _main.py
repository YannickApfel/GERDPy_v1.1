# -*- coding: utf-8 -*-
""" GERDPy - '_main.py'

    Steuerungsfile des Auslegungstools für das Projekt GERDI

    Legende:
        - Temperaturen:
            - T in Kelvin [K] - für (kalorische) Gleichungen
            - Theta in Grad Celsius [°C] - Input aus dem Wetterdatenfile

    Autor(en): Yannick Apfel
"""
import sys
import matplotlib.pyplot as plt
import time as tim
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
from scipy.constants import pi

# import GERDPy-modules
import boreholes, heatpipes, heating_element, gfunction, load_aggregation
from load_generator import *
from R_th import *
from weather_data import get_weather_data
from geometrycheck import check_geometry
from utilities import Q_moving_average


def main():

    # -------------------------------------------------------------------------
    # 1.) Parametrierung der Simulation (Geometrien, Stoffwerte, etc.)
    # -------------------------------------------------------------------------

    # 1.0) Standort
    h_NHN = 469                                     # Höhe über Normal-Null des Standorts  [m]
    path_wd = \
    './data/Wetterdaten_Hamburg_h.xlsx'    # Dateipfad der Wetterdaten-Datei definieren

    # 1.1) Erdboden
    a = 1.0e-6                                      # Temperaturleitfähigkeit [m²/s]
    lambda_g = 2.0                                  # Wärmeleitfähigkeit [W/mK]
    Theta_g = 10.0                                  # ungestörte Bodentemperatur [°C]

    # 1.2) Erdwärmesondenfeld

    # Geometrie-Import (.txt-Datei)
    boreField = boreholes.field_from_file('./data/custom_field_Eching.txt')

    # Sondenmeter (gesamt)  [m]
    H_field = boreholes.length_field(boreField)

    # Layout-Plot des Erdwärmesondenfelds
    boreholes.visualize_field(boreField)

    # 1.3) Bohrloch

    # Geometrie
    N = 10                                          # Anzahl Heatpipes pro Bohrloch [-]
    r_b = boreField[0].r_b                          # Radius der Erdwärmesondenbohrung [m]
    r_w = 0.07                                      # Radius der Wärmerohr-Mittelpunkte [m]
    r_iso_b = 0.006                                 # Außenradius der Isolationsschicht [m]
    r_pa = 0.006                                    # Außenradius des Wärmerohrs [m]
    r_pi = 0.0053                                   # Innenradius des Wärmerohrs [m]

    # Wärmeleitfähigkeiten [W/mK]
    lambda_b = 2                                    # Hinterfüllung
    lambda_iso = 0.03                               # Isolationsschicht
    lambda_p = 14                                   # Heatpipe

    # Geometrie-Erstellung
    hp = heatpipes.Heatpipes(N, r_b, r_w, r_iso_b, r_pa, r_pi, lambda_b,
                             lambda_iso, lambda_p)
    # Layout-Plot der Wärmerohrkonfiguration
    hp.visualize_hp_config()

    # 1.4) Anbindung zum Heizelement (zusätzliche Größen)

    # Geometrie
    D_iso_An = 0.005                                   # Dicke der Isolationsschicht [m]
    r_iso_An = r_pa + D_iso_An                         # Außenradius der Isolationsschicht [m]

    # Länge der Anbindungen zwischen Bohrlöchern und Heizelement (ab Geländeoberkante) [m]
    ''' l_An * N ergibt die Gesamtlänge an Heatpipe im Bereich der Anbindung
    '''
    l_An = 5

    # 1.5) Heizelement

    # Fläche Heizelement [m2]
    A_he = 35

    # minimaler Oberflächenabstand [m]
    x_min = .015

    # Wärmeleitfähigkeit des Betonelements [W/mk]
    lambda_Bet = 2.1

    # Mittelachsabstand der Kondensatrohre im Heizelement [m]
    s_R = .050

    # Gesamtlänge im Heizelement verbauter Kondensatrohre [m]
    l_R = 1000

    # Betondicke des Heizelements [m]
    D_he = 0.25
    
    # Dicke der Isolationsschicht an der Unterseite [m]
    D_iso_he = 0.03

    # Geometrie-Erstellung
    he = heating_element.HeatingElement(A_he, x_min, lambda_Bet, lambda_p, 
                                        2 * r_pa, 2 * r_pi, s_R, l_R, 
                                        D_he, D_iso_he)

    # 2.) Simulation

    # Simulationsparameter
    ''' dt darf 3600s (eine Stunde) aufgrund des Gültigkeitsbereichs der G-Functions
        nicht unterschreiten
    '''
    dt = 3600.                                   # Zeitschrittweite [s]
    sim_multiple_years = True
    a_sim = 4
    if sim_multiple_years:
        tmax = a_sim * 365 * 24 * dt  # total simulation time (multiple years) [s]
    else:
        tmax = 0.1 * 1 * (8760./12) * 3600.   # Gesamt-Simulationsdauer [s]  # total simulation time [s] (default: 730 h * 3600 s)
    Nt = int(np.ceil(tmax / dt))  # number of time steps [-]
        
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
    # 3.) Ermittlung thermischer Widerstände von Bohrlochrand bis Oberfläche
    # -------------------------------------------------------------------------

    # Gesamtsystem
    R_th = R_th_c(boreField) + R_th_b(lambda_g, boreField, hp) + \
        R_th_hp(boreField, hp) + R_th_he(he)

    # Boden bis Heatpipes
    R_th_ghp = R_th_c(boreField) + R_th_b(lambda_g, boreField, hp) + \
        R_th_hp(boreField, hp)

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
    u_inf, Theta_inf, S_w, B, Phi, RR = get_weather_data(Nt, path_wd)
    ''' u_inf - Windgeschwindigkeit [m/s]
        Theta_inf - Umgebungstemperatur [°C]
        S_w - Schneefallrate (Wasserequivalent) [mm/h]
        B - Bewölkungsgrad [octal units / 8]
        Phi - rel. Luftfeuchte [-]
        RR - Niederschlagsmenge (gesamt) [mm/h]
    '''

    # -------------------------------------------------------------------------
    # 6.) Iterationsschleife (Simulation mit Nt Zeitschritten der Länge dt)
    # -------------------------------------------------------------------------

    time = 0.
    i = -1
    start_sb = False                                # sb - snow balancing

    # Initialisierung Temperaturen [°C]
    ''' Vektoren:
            P[i] - Entzugsleistung f. Zeitschritt i 
            Theta_surf[i] - Oberflächentemp. f. Zeitschritt i
    '''
    Theta_b = np.zeros(Nt)                          # Bohrlochrand
    Theta_surf = np.zeros(Nt)                       # Oberfläche Heizelement
    
    # Initialisierung Vektor für Restwassermenge [kg]
    m_Rw = np.zeros(Nt)
    
    # Initialisierung Vektor für Restschneemenge [kg]
    m_Rs = np.zeros(Nt)

    # Initialisierung Gesamt-Entzugsleistung, Nutzleistung und Verluste [W]
    Q = np.zeros(Nt)
    Q_N = np.zeros(Nt)
    Q_V = np.zeros(Nt)
    
    # temp
    # Q_lat = np.zeros(Nt)
    # Q_sen = np.zeros(Nt)
    # Q_con = np.zeros(Nt)
    # Q_rad = np.zeros(Nt)
    Q_eva = np.zeros(Nt)

    # Hilfsgrößen
    start_sb_counter = np.zeros(Nt)
    sb_active = np.zeros(Nt)
    sim_mod = np.zeros(Nt)
    
    print('-----------------Simulation gestartet-----------------')

    while time < tmax:                              # Iterationsschleife (ein Durchlauf pro Zeitschritt)
        
        # Zeitschritt um 1 inkrementieren
        if start_sb == False:
            time += dt
            i += 1

        LoadAgg.next_time_step(time)

        # Ermittlung der Entzugsleistung im 1. Zeitschritt
        if i == 0:  # Annahme Theta_b = Theta_surf = Theta_g, Heizelementoberfläche trocken und schneefrei
            Q[i], Q_N[i], Q_V[i], calc_T_surf, Theta_surf[i], m_Rw[i], m_Rs[i], sb_active[i], sim_mod[i], Q_eva[i] = \
                load(h_NHN, u_inf[i], Theta_inf[i], S_w[i], he, Theta_g,
                     R_th, R_th_ghp, Theta_g, B[i], Phi[i], RR[i], 0, 0, start_sb, 
                     l_An * N, lambda_p, lambda_iso, r_iso_An, r_pa, r_pi)

        # Ermittlung der Entzugsleistung im Zeitschritt 2, 3, ..., Nt
        if i > 0:
            Q[i], Q_N[i], Q_V[i], calc_T_surf, Theta_surf[i], m_Rw[i], m_Rs[i], sb_active[i], sim_mod[i], Q_eva[i] = \
                load(h_NHN, u_inf[i], Theta_inf[i], S_w[i], he, Theta_b[i-1], 
                     R_th, R_th_ghp, Theta_surf[i-1], B[i], Phi[i], RR[i], m_Rw[i-1], m_Rs[i-1], start_sb, 
                     l_An * N, lambda_p, lambda_iso, r_iso_An, r_pa, r_pi)
                
        # Erhöhung der ermittelten Entzugsleistung um die Verluste an Anbindung (An) und Unterseite des Heizelements (He)
        Q[i] += Q_V[i]

        start_sb = False  # Variable Start-Schneebilanzierung zurücksetzen

        # Aufprägung der ermittelten Entzugsleistung auf Bodenmodell ('load_aggregation.py')
        LoadAgg.set_current_load(Q[i]/H_field)

        # Temperatur am Bohrlochrand [°C]
        deltaTheta_b = LoadAgg.temporal_superposition()
        Theta_b[i] = Theta_g - deltaTheta_b

        # Temperatur an der Oberfläche des Heizelements [°C]
        ''' Theta_surf wird hier nur ermittelt, falls Q. >= 0 (positive Entzugsleistung aus Boden), sonst
            erfolgt deren Ermittlung in 'load_generator.py' mit Hilfe einer vereinfachten Leistungsbilanz an der Oberfläche.
        '''
        if calc_T_surf is False:
            Theta_surf[i] = Theta_b[i] - Q[i] * R_th  # Oberflächentemp.

        # Schneebilanzierung starten
        ''' Zeitschritt wird einmalig wiederholt, falls sich eine Schneeschicht beginnt zu bilden:
            - Oberflächentemperatur Theta_surf[i] < 0 UND
            - Schneefallrate S_w[i] > 0 UND
            - Restschneemenge m_Rs[i] == 0 (noch keine Restschneemenge vorhanden)
            müssen erfüllt sein

            => Schnee bleibt liegen: Schneebilanz an Oberfläche
        '''
        if (Theta_surf[i] < 0 and S_w[i] > 0 and m_Rs[i] == 0):
            start_sb = True
            start_sb_counter[i] = 1

        # Konsolenausgabe des momentanen Zeitschritts
        print(f'Zeitschritt {i+1} von {Nt}')

    # Zeitstempel (Simulationsdauer) [s]
    toc = tim.time()
    print('Total simulation time: {} sec'.format(toc - tic))

    # -------------------------------------------------------------------------
    # 7.) Energiekennzahlen
    # -------------------------------------------------------------------------
    ''' Q_m                 - zeitlich gemittelte Entzugsleistung [W]
        E                   - Gesamtenergiemenge, die dem Boden entzogen wurde [MWh]
        f_N [%] = E_N / E   - Nutzenergiefaktor [%]
            E_N             - zur sensiblen Erwärmung und latenten Schmelze des Schnees verwendete Energie
            E - E_N         - Summe der Verluste durch Konvektion, Strahlung und Verdunstung
    '''

    # 24h-gemittelte Entzugsleistung (gleitender Mittelwert)
    Q_ma = Q_moving_average(Q)

    # Gesamtenergiemenge [MWh]
    E = (np.sum(Q) / len(Q)) * Nt * 1e-6
    print('-----------------Simulation beendet-----------------')
    print(f'Dem Boden wurden {round(E, 4)} MWh entnommen')

    # Nutzenergiefaktor [%]
    # f_N = (np.sum(Q_N) / len(Q_N)) / (np.sum(Q) / len(Q)) * 100
    # print(f'Davon wurden {round(f_N, 2)} % als Nutzenergie zur Schneeschmelze aufgewendet. Der Rest sind Verluste an der Ober- und Unterseite des Heizelements sowie dessen Anbindungsleitungen.')
    # print(50*'-')
    
    # -------------------------------------------------------------------------
    # 8.) Plots
    # -------------------------------------------------------------------------

    # x-Achse aller Plots (Simulationsstunden) [h]
    hours = np.array([(j+1)*dt/3600. for j in range(Nt)])

    # -------------------------------------------------------------------------
    # 8.1) Figure 1 & 2 (Plots für End-User)
    # -------------------------------------------------------------------------

    plt.rc('figure')
    fig1 = plt.figure()
    fig2 = plt.figure()

    font = {'weight': 'bold', 'size': 10}
    plt.rc('font', **font)

    # Lastprofil {Entzugsleistung - Entzugsleistung (gleitender Mittelwert 24h) - Verdunstungsleistung - Verluste (Anbindung + Unterseite Heizelement)}
    ax1 = fig1.add_subplot(411)
    ax1.set_ylabel(r'$q$ [W/m2]')
    ax1.plot(hours, Q / A_he, 'k-', lw=1.2)
    ax1.plot(hours, Q_ma / A_he, 'r--', lw=1.2)
    ax1.plot(hours, Q_V / A_he, 'g-', lw=1.2)
    ax1.plot(hours, Q_eva / A_he, 'y-', lw=1.2)
    ax1.legend(['Entzugsleistung', 'Entzugsleistung-24h-gemittelt',
                'Verluste (Anbindung + Unterseite Heizelement)', 'Verdunstung'],
               prop={'size': font['size'] - 4}, loc='upper left')
    ax1.grid('major')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())

    # Schneefallrate - Schneehöhe - Umgebungstemperatur - Windgeschwindigkeit
    ax2 = fig1.add_subplot(412)
    ax2.set_ylabel('Schneefallrate [mm/h] \n Schneehöhe [H2O-mm]')
    ax2.plot(hours, S_w, 'b-', lw=0.8)
    ax2.plot(hours, m_Rs / (A_he * (997 / 1000)), 'g-', lw=0.8)
    ax2.legend(['Schneefallrate', 'Schneehöhe'],
               prop={'size': font['size'] - 4}, loc='upper left')
    ax2.grid('major')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())

    # Umgebungstemperatur - Windgeschwindigkeit
    ax3 = fig1.add_subplot(413)
    ax3.set_ylabel('$T$ [degC] \n Windgeschwindigkeit [m/s]')
    ax3.plot(hours, Theta_inf, 'k-', lw=0.8)
    ax3.plot(hours, u_inf, 'm--', lw=0.8)
    ax3.legend(['Umgebungstemperatur', 'Windgeschwindigkeit'],
                 prop={'size': font['size'] - 4}, loc='upper right')
    ax3.grid('major')
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax3.yaxis.set_minor_locator(AutoMinorLocator())

    # Temperaturverläufe Bohrlochrand und Oberfläche Heizelement
    ax4 = fig1.add_subplot(414)
    ax4.set_xlabel(r'$t$ [h]')
    ax4.set_ylabel(r'$T$ [degC]')
    ax4.plot(hours, Theta_b, 'r-', lw=1.2)
    ax4.plot(hours, Theta_surf, 'c-', lw=0.6)
    ax4.legend(['T_Bohrlochrand', 'T_Oberflaeche'],
               prop={'size': font['size'] - 4}, loc='upper right')
    ax4.grid('major')
    ax4.xaxis.set_minor_locator(AutoMinorLocator())
    ax4.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Borehole wall temperature annual stacked curves
    # fig2
    if sim_multiple_years:
        ax5 = fig2.add_subplot(211)
        ax5.set_xlabel(r'$t$ [h]')
        ax5.set_ylabel(r'$T$ [degC]')
        colour_map = iter(plt.cm.gist_rainbow(np.linspace(0, 1, a_sim)))
        for j in range(a_sim):
            ax5.plot(np.arange(1, 8761, 1, dtype=int), Theta_b[(0+j*8760):(8760+j*8760)], c=next(colour_map), 
                     lw=0.7, label=f'Borehole wall temperature - Year {j+1}')
    else:
        ax5 = fig2.add_subplot(111)
        ax5.set_xlabel(r'$t$ [h]')
        ax5.set_ylabel(r'$T$ [degC]')
        ax5.plot(hours, Theta_b, 'r-', lw=1.2, label='Borehole wall temperature - Year 1')
    ax5.legend(prop={'size': font['size'] - 2}, loc='lower center')
    ax5.grid('major')
    ax5.xaxis.set_minor_locator(AutoMinorLocator())
    ax5.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Borehole wall temperature at beginning of heating period
    # fig2
    if sim_multiple_years:
        ax6 = fig2.add_subplot(212)
        ax6.set_xlabel(r'$Year$ [a]')
        ax6.set_ylabel(r'$T$ [degC]')
        ax6_x = np.arange(0, a_sim+1, 1, dtype=int)
        # straight connecting lines:
        ax6.plot([ax6_x[0], ax6_x[1]], [Theta_g, Theta_b[5860]], 'b', linewidth=1)
        for j in range(a_sim):
            if j<(a_sim-1):
                ax6.plot([ax6_x[j+1], ax6_x[j+2]], [Theta_b[5860+j*8760], Theta_b[5860+(j+1)*8760]], 'b', linewidth=1, label='_nolegend_')
            else:
                break
        # scatter plot:
        ax6.plot(ax6_x[0], Theta_g, marker='x', markersize=10, markeredgecolor='green', label='Undisturbed ground temperature')
        for j in range(a_sim):
            ax6.plot(ax6_x[j+1], Theta_b[5860+j*8760], color='red', marker='o', markersize=10, markeredgewidth=0.0)
            if j==0:
                ax6.plot(ax6_x[j+1], Theta_b[5860+j*8760], color='red', marker='o', markersize=10, markeredgewidth=0.0, 
                         label='Borehole wall temperature at the beginning of the heating period (01.09.)')
        ax6.legend(prop={'size': font['size'] - 2}, loc='upper right')
        ax6.grid('major')
        ax6.xaxis.set_minor_locator(AutoMinorLocator())
        ax6.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax6.yaxis.set_minor_locator(AutoMinorLocator())
    
    # -------------------------------------------------------------------------
    # 8.2) Figure 2 (zusätzliche Plots)
    # -------------------------------------------------------------------------

    plt.rc('figure')
    fig3 = plt.figure()

    font = {'weight': 'bold', 'size': 10}
    plt.rc('font', **font)

    # Darstellungen Simulationsmodus
    ax7 = fig3.add_subplot(311)
    ax7.plot(hours, start_sb_counter, 'k--', lw=1.5)
    ax7.plot(hours, sb_active, 'g-', lw=1.3)
    ax7.plot(hours, sim_mod, 'y-', lw=1.3)
    ax7.plot(hours, sim_mod, 'y-', lw=1.3)    
    ax7.legend(['sb_active', 'sim_mod'],
                prop={'size': font['size'] - 4}, loc='upper right')
    ax7.grid('major')
    ax7.xaxis.set_minor_locator(AutoMinorLocator())
    ax7.yaxis.set_minor_locator(AutoMinorLocator())
    
    # Wasser- und Schneebilanzlinie (Wasserequivalent)
    ax8 = fig3.add_subplot(312)
    ax8.set_ylabel('[mm]')
    ax8.plot(hours, m_Rw / (A_he * (997 / 1000)), 'b-', lw=0.8)
    ax8.plot(hours, m_Rs / (A_he * (997 / 1000)), 'g-', lw=0.8)
    ax8.legend(['Wasserhoehe', 'Schneehoehe'],
                prop={'size': font['size'] - 4}, loc='upper left')
    ax8.grid('major')
    ax8.xaxis.set_minor_locator(AutoMinorLocator())
    ax8.yaxis.set_minor_locator(AutoMinorLocator())

    # Temperaturverläufe Bohrlochrand und Oberfläche Heizelement
    ax9 = fig3.add_subplot(313)
    ax9.set_xlabel(r'$t$ [h]')
    ax9.set_ylabel(r'$T$ [degC]')
    ax9.plot(hours, Theta_b, 'r-', lw=1.2)
    ax9.plot(hours, Theta_surf, 'c-', lw=0.6)
    ax9.legend(['T_Bohrlochrand', 'T_Oberflaeche'],
                prop={'size': font['size'] - 4}, loc='upper right')
    ax9.grid('major')
    ax9.xaxis.set_minor_locator(AutoMinorLocator())
    ax9.yaxis.set_minor_locator(AutoMinorLocator())

    # plt.tight_layout()
    plt.show()
    
    return


# Main function
if __name__ == '__main__':
    main()
