# -*- coding: utf-8 -*-
""" Autor: Yannick Apfel

    Modul zur Prüfung geometrischer Verträglichkeit:
        - Abstand der Erwärmesondenbohrungen zueinander muss > 2 * r_b
        - identische Außenradien der Erdwärmesondenbohrungen
        - Abstände der Heatpipes in der Bohrung
"""
import numpy as np


def check_geometry(borefield, heatpipes):
    # errors_exist mit 0 initialisieren
    errors_exist = 0
    
    # Test auf Bohrloch-Duplikate
    duplicate_pairs = []
    for i in range(len(borefield)):
        borehole_1 = borefield[i]
        for j in range(i, len(borefield)):  # Schleife
            borehole_2 = borefield[j]
            if i == j:  # kein Vergleich mit sich selbst
                continue
            else:
                dist = borehole_1.distance(borehole_2)
            if abs(dist - borehole_1.r_b) < borehole_1.r_b:
                duplicate_pairs.append((i+1, j+1))
    if duplicate_pairs:
        print(f'Geometric conflict between the following borehole pairs:')
        print(*duplicate_pairs, sep = ", ")
        print('Please adjust!')
        errors_exist = 1

    # Test: Bohrlochradien müssen identisch sein
    for i in range(len(borefield) - 1):
        if borefield[i].r_b != borefield[i+1].r_b:
            print('Borehole radii must be identical!')
            print('Please adjust!')
            errors_exist = 1
            break

    # Heatpipe-Geometrie prüfen
    xy = heatpipes.xy_mat()
    heatpipe_distance = np.sqrt((xy[1, 0] - xy[2, 0])**2 + (xy[1, 1] - xy[2, 1])**2)
    if ((heatpipes.r_pa + heatpipes.r_w) >= heatpipes.r_b):
        print("Heatpipe circle too big!")
        errors_exist = 1
    elif heatpipe_distance <= (2 * heatpipes.r_pa):
        print("Too many heat pipes!")
        errors_exist = 1

    return errors_exist
