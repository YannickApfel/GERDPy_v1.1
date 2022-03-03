# -*- coding: utf-8 -*-
""" GERDPy - 'R_th_tot.py'
    
    Modul für den thermischen Gesamtwiderstand des Systems
    
    Reihenschaltung aus Teilwiderständen:
        - R_th_tot: Gesamtwiderstand Boden-Umgebung: 
            = R_th_c + R_th_b + R_th_hp + R_th_he
        - R_th_g_hp: Teil-Widerstand Boden-Heatpipe:
            = R_th_c + R_th_b + R_th_hp = R_th_tot - R_th_he

    Autor(en): Yannick Apfel
"""
from R_th_c import R_th_c
from R_th_b import R_th_b
from R_th_hp import R_th_hp
from R_th_he import R_th_he


def R_th_tot(lambda_g, borefield, hp, he):  # [K/W]

    R_th_tot = R_th_c(borefield) + R_th_b(lambda_g, borefield, hp) + \
        R_th_hp(borefield, hp) + R_th_he(he)
    # + weitere thermische Widerstände

    return R_th_tot


def R_th_g_hp(lambda_g, borefield, hp):  # [K/W]
    
    R_th_g_hp = R_th_c(borefield) + R_th_b(lambda_g, borefield, hp) + \
        R_th_hp(borefield, hp)
        
    return R_th_g_hp
