# -*- coding: utf-8 -*-
""" Bohrlochwiderstand
    [Hellström 1991 - Line source approximation]

    N Heatpipes gleichmäßig kreisförmig im Bohrloch angeordnet

    Autor: Yannick Apfel
"""


def R_th_b(lambda_inf, borefield, hp):

    import math
    import numpy as np
    from numpy.linalg import inv
    from scipy.constants import pi
    from boreholes import length_field

    # %% 1.) Parameter

    # Geometrie Erdwärmesonde
    H_field = length_field(borefield)   # Gesamtlänge Sondenfeld [m]
    r_b = borefield[0].radius()         # Radius Bohrloch [m]

    # Geometrie Heatpipes
    N = hp.N                            # Anzahl Heatpipes im Bohrloch [-]
    r_pa = hp.r_pa                      # Außenradius Wärmerohr (mit Iso) [m]
    r_iso = hp.r_iso                    # Innenradius Ummantelung [m]
    r_pi = hp.r_pi                      # Innenradius Wärmerohr [m]

    # Wärmeleitfähigkeiten [W/mK]:
    # lambda_inf (importiert)
    lambda_b = hp.lambda_b              # Wärmeleitfähigkeit Verfüllung
    lambda_iso = hp.lambda_iso          # Wärmeleitfähigkeit Iso
    lambda_p = hp.lambda_p              # Wärmeleitfähigkeit Heatpipe

    # %% 2a.) Koordinaten der Wärmerohre im System Bohrloch

    xy = hp.xy_mat()  # 1. Spalte: x, 2. Spalte: y

    # %% 2b.) Hilfsgrößen

    # Verhältnis der Wärmeleitfähigkeiten
    sigma = (lambda_b - lambda_inf) / (lambda_b + lambda_inf)

    # Übergangswiderstand Wärmerohr + Isolationsschicht
    r_pm = math.log(r_pa / r_iso) / (2 * pi * lambda_iso) + \
        math.log(r_iso / r_pi) / (2 * pi * lambda_p)
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

    # %% 3.) Ermittlung Bohrlochwiderstand

    R_th_b = (sum(sum(inv(R_mn_0))) * H_field) ** -1

    return R_th_b
