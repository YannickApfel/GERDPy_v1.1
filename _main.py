# -*- coding: utf-8 -*-
""" GERDPy - Main-File
    Steuerungsfile des Auslegungstools für GERDI

    Autor: Yannick Apfel

    todos:
        - GUI für End-User
        - Lastmodell implementieren
"""
import sys
import matplotlib.pyplot as plt
import time as tim
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from scipy.constants import pi

# import GERDPy-modules
import boreholes, heatpipes, heating_element, gfunction, load_aggregation
from load_generator import load
from geometrycheck import check_geometry
from R_th_tot import R_th_tot


def main():
    # -------------------------------------------------------------------------
    # 1.) Parametrierung der Simulation (Geometrien, Stoffwerte, etc.)
    # -------------------------------------------------------------------------

    # 1.1) Erdboden
    a = 1.0e-6                  # Temperaturleitfähigkeit [m2/s]
    lambda_inf = 2.0            # Wärmeleitfähigkeit [W/mK]
    T_g = 10.0                  # ungestörte Bodentemperatur [degC]

    # 1.2) Erdwärmesondenfeld

    # Geometrie-Import (.txt-Datei)
    boreField = boreholes.field_from_file('./data/custom_field.txt')

    # Sondenmeter (gesamt)
    H_field = boreholes.length_field(boreField)

    # Layout-Plot des Erdwärmesondenfelds
    boreholes.visualize_field(boreField)

    # 1.3) Bohrloch

    # Geometrie
    N = 6                        # Anzahl Heatpipes pro Bohrloch [-]
    r_b = boreField[0].radius()  # Radius der Erdwärmesondenbohrung [m]
    r_w = 0.12                # Radius der Wärmerohr-Mittelpunkte [m]
    r_pa = 0.016                 # Außenradius der Isolationsschicht [m]
    r_iso = 0.016                # Innenradius der Isolationsschicht [m]
    r_pi = 0.015                 # Innenradius des Wärmerohrs [m]

    # Wärmeleitfähigkeiten
    lambda_b = 2                # lambda der Hinterfüllung [W/mK]
    lambda_iso = 0.3            # lambda der Isolationsschicht [W/mK]
    lambda_p = 14               # lambda der Heatpipe [W/mK]

    # Geometrie-Erstellung
    hp = heatpipes.Heatpipes(N, r_b, r_w, r_pa, r_iso, r_pi, lambda_b,
                             lambda_iso, lambda_p)
    # Layout-Plot der Wärmerohrkonfiguration
    hp.visualize_hp_config()

    # 1.4) Heizelement

    # Fläche Heizelement [m2]
    A_he = 165

    # minimaler Oberflächenabstand [mm]
    x_min = 25

    # Wärmeleitfähigkeit des Betonelements [W/mk]
    lambda_Bet = 2.3

    # Geometrie-Erstellung
    he = heating_element.HeatingElement(A_he, x_min, lambda_Bet)

    # 2.) Simulation

    # Simulationsparameter
    dt = 3600.                           # Zeitschrittweite [s]
    tmax = 0.5 * 1 * (8760./12) * 3600.  # Gesamt-Simulationsdauer [s]
    Nt = int(np.ceil(tmax/dt))           # Anzahl Zeitschritte [-]

    # -------------------------------------------------------------------------
    # 2.) Überprüfung der geometrischen Verträglichkeit (Sonden & Heatpipes))
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

    R_th = R_th_tot(lambda_inf, boreField, hp, he)

    # -------------------------------------------------------------------------
    # 4.) Ermittlung der G-Function
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
    LoadAgg.initialize(gFunc/(2*pi*lambda_inf))

    # -------------------------------------------------------------------------
    # 5.) Iterationsschleife (Simulation mit Nt Zeitschritten der Länge dt)
    # -------------------------------------------------------------------------

    time = 0.
    i = -1

    # Initialsierung Temperaturvektoren (ein Eintrag pro Zeitschritt)
    T_b = np.zeros(Nt)
    # T_grout = np.zeros(Nt)
    # T_hp = np.zeros(Nt)
    # T_he = np.zeros(Nt)
    T_surf = np.zeros(Nt)

    Q = np.zeros(Nt)

    while time < tmax:  # Iterationsschleife (ein Durchlauf pro Zeitschritt)

        # Zeitschritt um 1 inkrementieren
        time += dt
        i += 1
        LoadAgg.next_time_step(time)

        # Ermittlung der Entzugsleistung mit 'load_generator.py'
        Q[i] = load(time/3600.)

        # Aufprägung der ermittelten Entzugsleistung mit 'load_aggregation.py'
        LoadAgg.set_current_load(Q[i]/H_field)

        # Temperatur am Bohrlochrand
        deltaT_b = LoadAgg.temporal_superposition()
        T_b[i] = T_g - deltaT_b

        # Temperaturkaskade (T-Verlust bis zur Oberfläche)
        # Temperatur am Bohrlochrand (Hinterfüllmaterial-seitig)
        # T_grout[i] = T_b[i] - Q[i] * R_th_c(boreField)

        # Temperatur an den Wärmerohren im Verdampfer
        # T_hp[i] = T_grout[i] - Q[i] * R_th_b(lambda_inf, boreField, hp)

        # Temperatur an Wärmerohren im Kondensator (Heizelement)
        # T_he[i] = T_hp[i] - Q[i] * R_th_hp(boreField, hp)

        # Temperatur an der Oberfläche des Heizelements
        T_surf[i] = T_b[i] - Q[i] * R_th

    # -------------------------------------------------------------------------
    # 6.) Plots
    # -------------------------------------------------------------------------

    # Zeitstempel (Simulationsdauer)
    toc = tim.time()
    print('Total simulation time: {} sec'.format(toc - tic))

    plt.rc('figure')
    fig = plt.figure()

    font = {'weight': 'bold', 'size': 22}
    plt.rc('font', **font)

    # Lastprofil (thermische Leistung Q. über die Simulationsdauer)
    ax1 = fig.add_subplot(211)
    ax1.set_xlabel(r'$t$ (hours)')
    ax1.set_ylabel(r'$Q$ (W)')
    hours = np.array([(j+1)*dt/3600. for j in range(Nt)])
    ax1.plot(hours, Q, 'b-', lw=1.5)  # plot

    # Temperaturverläufe
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel(r'$t$ (hours)')
    ax2.set_ylabel(r'$T_b$ (degC)')
    # plots
    ax2.plot(hours, T_b, 'r-', lw=2.0)
    # ax2.plot(hours, T_grout, 'y-', lw=0.5)
    # ax2.plot(hours, T_hp, 'm-', lw=0.5)
    # ax2.plot(hours, T_he, 'g-', lw=0.5)
    ax2.plot(hours, T_surf, 'c-', lw=1.0)

    ax2.legend(['T_borehole-wall', 'T_surface'],
               prop={'size': font['size'] - 15}, loc='upper right')

    # Beschriftung Achsenwerte
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    # Fenstergröße anpassen
    # plt.tight_layout()

    return


# Main function
if __name__ == '__main__':
    main()
