# -*- coding: utf-8 -*-
""" GERDPy - 'R_th.py'
    
    Modul für die thermischen Widerstände im System
    
    1.) R_th_c - Kontaktwiderstand Erdboden-Bohrlochrand (Hinterfüllung)

        - nach [VDI-Wärmeatlas 2013 - Wärmeübergangskoeffizient Wand-Schüttung]
        
    2.) R_th_b - Bohrlochwiderstand

        - nach [Hellström 1991 - Line source approximation]
        - N Heatpipes gleichmäßig kreisförmig im Bohrloch angeordnet
        
    3.) R_th_hp - thermischer Widerstand des Thermosiphons (Heatpipe)

        - Thermischer Widerstand im Thermosiphon (verdampferseitig) mit
          simpler Charakteristik im Q.-deltaT - Diagramm (linearer Anstieg)
        - Vernachlässigung von Leistungsgrenzen (z. B: Entrainment-Limit)
        - Heatpipe-Charakteristik: Q. = (500 W * N_Anzahl_Heatpipes / 1 K) * deltaT
        
    4.) R_th_he - thermischer Widerstand des Heizelements

        - analytisch ermittelter therm. Widerstand basierend auf der Norm
          VDI 2055-1 [2019] für die thermischen Verluste gedämmter Rohrleitungen
          verlegt in Registeranordnung in ein- oder mehrschichtigen Platten
          (z. B. Fußbodenheizungen) - vereinfacht und angepasst in [Apfel2020]
        - variable Größen:
            - minimaler Oberflächenabstand x_min [m]
            - Wärmeleitfähigkeiten Beton und Kondensatrohre [W/mK]
            - Rohrgeometrie des Kondensator (Durchmesser, Rohrabstand,
                                             Rohrlänge) [m]

    Autor(en): Yannick Apfel
"""
# R_th_c - Kontaktwiderstand Erdboden-Bohrlochrand (Hinterfüllung)
def R_th_c(borefield):

    import math, boreholes
    from scipy.constants import pi, R, Stefan_Boltzmann
    from boreholes import length_field

    # 1.) Stoffwerte und Parameter

    phi = 0.8               # Flächenbedeckungsgrad [-]
    lambda_g = 0.025        # Wärmeleitfähigkeit Gas [W/mK]
    d = 1e-3                # Partikeldurchmesser [m]
    delta = 250 * 1e-6      # Oberflächenrauhigkeit Partikel [m]
    C = 2.8                 # materialabhängige Konstante [-]
    M = 0.02896             # molare Masse Gas [kg/mol]
    T = 283                 # Temperatur der Kontaktzone [K]
    c_pg = 1007             # spez. Wärmekapazität Gas [J/kgK]
    epsilon_S = 0.2         # Emissionskoeffizient Schüttung [-]
    epsilon_W = 0.2         # Emissionskoeffizient Wand [-]
    p = 100000              # Gasdruck [Pa]

    # 2a.) Hilfsparameter

    # Akkomodationskoeffizient
    gamma = (10 ** (0.6 - (1000 / T + 1) / C) + 1) ** -1

    # freie Weglänge der Gasmoleküle

    l_frei = 2 * (2 - gamma) / gamma * math.sqrt(2 * pi * R * T / M) * \
        lambda_g / (p * (2 * c_pg - R / M))

    # Verhältnisse der Emissionsgrade
    C_WS = Stefan_Boltzmann / (1 / epsilon_W + 1 / epsilon_S - 1)

    # 2b.) Wärmeübergangskoeffizient Wand-Schüttung

    # Anteil Wärmeleitung
    alpha_WP = 4 * lambda_g / d * ((1 + 2 * (l_frei + delta) / d) *
                                   math.log(1 + d / (2 * (l_frei + delta)))
                                   - 1)

    # Anteil Strahlung
    alpha_rad = 4 * C_WS * T ** 3

    # Wärmeübergangskoeffizient Wand-Schüttung
    alpha_WS = phi * alpha_WP + alpha_rad

    # 3.) thermischer Kontaktwiderstand des Erwärmesondenfelds

    r_b = borefield[0].r_b
    H_field = length_field(borefield)

    R_th_c = (2 * math.pi * r_b * alpha_WS * H_field) ** -1

    return R_th_c


# R_th_b - Bohrlochwiderstand
def R_th_b(lambda_g, borefield, hp):

    import math
    import numpy as np
    from numpy.linalg import inv
    from scipy.constants import pi
    from boreholes import length_field

    # 1.) Parameter

    # Geometrie Erdwärmesonde
    H_field = length_field(borefield)   # Gesamtlänge Sondenfeld [m]
    r_b = borefield[0].r_b              # Radius Bohrloch [m]

    # Geometrie Heatpipes
    N = hp.N                            # Anzahl Heatpipes im Bohrloch [-]
    r_iso_a = hp.r_iso_a                # Außenradius Ummantelung [m]
    r_pa = hp.r_pa                      # Außenradius Wärmerohr [m]
    r_pi = hp.r_pi                      # Innenradius Wärmerohr [m]

    # Wärmeleitfähigkeiten [W/mK]:
    # lambda_g (importiert)
    lambda_b = hp.lambda_b              # Wärmeleitfähigkeit Verfüllung
    lambda_iso = hp.lambda_iso          # Wärmeleitfähigkeit Iso
    lambda_p = hp.lambda_p              # Wärmeleitfähigkeit Heatpipe

    # 2a.) Koordinaten der Wärmerohre im System Bohrloch

    xy = hp.xy_mat()  # 1. Spalte: x, 2. Spalte: y

    # 2b.) Hilfsgrößen

    # Verhältnis der Wärmeleitfähigkeiten
    sigma = (lambda_b - lambda_g) / (lambda_b + lambda_g)

    # Übergangswiderstand Wärmerohr + Isolationsschicht
    r_pm = math.log(r_iso_a / r_pa) / (2 * pi * lambda_iso) + \
        math.log(r_pa / r_pi) / (2 * pi * lambda_p)
    # r_pm = 0 (falls der thermische Widerstand von Iso und Rohr ignoriert werden soll)

    # koordinatenabhängige Hilfskoeffizienten
    b_m = lambda x_m, y_m: math.sqrt(x_m ** 2 + y_m ** 2) / r_b
    b_mn = lambda x_n, x_m, y_n, y_m: math.sqrt((x_n - x_m) ** 2 + (y_n - y_m) ** 2) / r_b
    b_mn_ = lambda b_m, b_n, b_mn: math.sqrt((1 - b_m ** 2) * (1 - b_n ** 2) + b_mn ** 2)

    # Koeffizientenmatrix des Bohrlochs

    R_mn_0 = np.zeros([N, N])

    # Koeffizientenmatrix befüllen
    for i in range(N):          # iterieren für m
        for j in range(N):      # iterieren für n
            if i == j:
                R_mn_0[i, j] = \
                    (2 * pi * lambda_b) ** -1 * (math.log(r_b / r_pa)
                    - sigma * math.log(1 - b_m(xy[i, 0], xy[i, 1]) ** 2)) + r_pm
            else:
                R_mn_0[i, j] = \
                    - (2 * pi * lambda_b) ** -1 * (math.log(b_mn(xy[j, 0], xy[i, 0], xy[j, 1], xy[i, 1]))
                    - sigma * math.log(b_mn_(b_m(xy[i, 0], xy[i, 1]), b_m(xy[j, 0], xy[j, 1]), b_mn(xy[j, 0], xy[i, 0], xy[j, 1], xy[i, 1]))))

    # 3.) Ermittlung Bohrlochwiderstand

    R_th_b = (sum(sum(inv(R_mn_0))) * H_field) ** -1

    return R_th_b


# R_th_hp - thermischer Widerstand des Thermosiphons (Heatpipe)
def R_th_hp(borefield, hp):

    # 1.) Parameter

    N = hp.N                                # Anzahl Heatpipes pro Bohrloch [-]
    N_b = len(borefield)                    # Anzahl Bohrlöcher

    # Steigung der Charakteristik (delta_y / delta_x)
    delta_y = 500       # Leistung pro Heatpipe [W]
    delta_x = 1         # delta_T [K]

    # 2.) thermischer Widerstand (Entrainment wird vernachlässigt)

    R_th_hp = delta_x / (delta_y * N * N_b)

    return R_th_hp


# R_th_he - thermischer Widerstand des Heizelements
def R_th_he(he):  # thermischer Widerstand [K/W]

    from heating_element_utils import q_l

    # 1.) zusätzliche Parameter des Heizelements
    x_o = he.x_min + 0.5 * he.d_R_a  # Achsabstand Kondensatrohr nach oben [m]
    x_u = he.D - x_o  # Achsabstand Kondensatrohr nach unten [m]

    # 2.) delta-T zwischen Kondensatrohren und Oberfläche zu 1 K definieren
    Theta_R = 1
    Theta_inf_o = 0

    # 3.) thermischer Widerstand des Heizelements in [K/W] --> R_th = 1 K / (q_l * l_R)
    ''' andere Seite wird als thermisch isoliert betrachtet: state_u_insul=True
    '''
    R_th_he = (he.l_R * q_l(x_o, x_u, he.d_R_a, he.d_R_i, he.lambda_B, he.lambda_R, he.s_R, Theta_R, Theta_inf_o, 
                            state_u_insul=True)) ** -1

    return R_th_he
