# -*- coding: utf-8 -*-
""" Import der notwendigen Wetterdaten aus .csv-File

    Autor: Yannick Apfel
"""


def get_weather_data(Nt):

    import numpy as np
    import pandas as pd

    # Dateipfad der Wetterdaten-Datei definieren
    path = './data/Wetterdaten_München-Riem_h.xlsx'

    # Wetterdaten importieren
    data = pd.read_excel(path)

    # 1.) Windgeschwindigkeit (Vektor mit Nt random floats zwischen 0 und 5) [m/s]
    # u_inf = np.ones(Nt)
    # u_inf = np.random.uniform(0, 5, Nt)
    u_inf = np.array(data.iloc[4:(Nt+5), 6], dtype='float')

    # 2.) Umgebungstemperatur (Vektor mit Nt random floats zwischen -2 und 3) [°C]
    # Theta_inf = np.ones(Nt) * -1
    # Theta_inf = np.random.uniform(-2, 3, Nt)
    # Theta_inf[50:60] = 12
    Theta_inf = np.array(data.iloc[4:(Nt+5), 4], dtype='float')

    # 3.) Schneefallrate (Vektor mit Nt random floats zwischen 0 und 1.5) [mm/h]
    # S_w = np.ones(Nt)
    # S_w = np.random.uniform(0, 1.5, Nt)
    S_w = np.array(data.iloc[4:(Nt+5), 3], dtype='float')
    # Einträge zu null setzen, falls Theta_inf >= 1 °C (Schnee fällt als Regen)
    for i,j in enumerate(Theta_inf):
        if j >= 1:
            S_w[i] = 0

    # 4.) Bewölkungsgrad (Vektor mit Nt ints zwischen 0 und 8) [-]
    # B = np.ones(Nt) * 3
    # B = np.random.randint(9, size=Nt)
    B = np.array(data.iloc[4:(Nt+5), 7], dtype='int') / 8

    # 5.) rel. Luftfeuchte (Vektor mit Nt floats zwischen 0.63 und 0.96) [%]
    # Phi = np.ones(Nt) * 0.7
    # Phi = np.random.uniform(0.63, 0.96, Nt)
    Phi = np.array(data.iloc[4:(Nt+5), 5], dtype='float') / 100

    return u_inf, Theta_inf, S_w, B, Phi
