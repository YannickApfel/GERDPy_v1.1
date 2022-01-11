# -*- coding: utf-8 -*-
""" Autor: Yannick Apfel
"""


class HeatingElement(object):
    """
    Enthält die Klasse für das Heizelement.

    Attribute
    ----------
    A_he:           float [m²]
                    Fläche des Heizelements.
    x_min:          float [m]
                    minimaler Oberflächenabstand der Verrohrung (Überdeckung)
    lambda_B:       float [W/mK]
                    Wärmeleitfähigkeit des Betonelements
    lambda_R:       float [W/mK]
                    Wärmeleitfähigkeit der Kondensatrohre
    d_R_a:          float [m]
                    Außendurchmesser eines Kondensatrohrs
    d_R_i:          float [m]
                    Innendurchmesser eines Kondensatrohrs
    s_R:            float [m]
                    Mittelachsabstand der Kondensatrohre des Rohrregisters
    l_R:            float [m]
                    Gesamtlänge im Heizelement verbauter Kondensatrohre
    D:              float [m]
                    Betondicke des Heizelements
    D_iso:          float [m]
                    Dicke der Isolationsschicht an der Unterseite des Heizelements
    """

    def __init__(self, A_he, x_min, lambda_B, lambda_R, d_R_a, d_R_i, s_R,
                 l_R, D, D_iso):
        self.A_he = float(A_he)
        self.x_min = float(x_min)
        self.lambda_B = float(lambda_B)
        self.lambda_R = float(lambda_R)
        self.d_R_a = float(d_R_a)
        self.d_R_i = float(d_R_i)
        self.s_R = float(s_R)
        self.l_R = float(l_R)
        self.D = float(D)
        self.D_iso = float(D_iso)

    def __repr__(self):
        s = ('HeatingElement(A_he={self.A_he}, x_min={self.x_min}, \
             lambda_B={self.lambda_B}, lambda_R={self.lambda_R},'
             ' d_R_a={self.d_R_a}, d_R_i={self.d_R_i}, \
                     s_R={self.s_R}, l_R={self.l_R}, D={self.D}, \
                         D_iso={self.D_iso}').format(self=self)
        return s
