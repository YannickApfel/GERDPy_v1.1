# -*- coding: utf-8 -*-
""" Gesamtwiderstand des Erdwärmesondenfelds
    Summe aus 4 Einzelwiderständen

    Autor: Yannick Apfel
"""
from R_th_c import R_th_c
from R_th_b import R_th_b
from R_th_hp import R_th_hp
from R_th_he import R_th_he


def R_th_tot(lambda_inf, borefield, hp, he):

    R_th_total = R_th_c(borefield) + R_th_b(lambda_inf, borefield, hp) + \
        R_th_hp(borefield, hp) + R_th_he(he)
    # + weitere thermische Widerstände

    return R_th_total
