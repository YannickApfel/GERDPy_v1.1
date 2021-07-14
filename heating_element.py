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
    x_min:          float [mm]
                    minimaler Oberflächenabstand der Verrohrung (Überdeckung)
    lambda_Bet:     float [W/mK]
                    Wärmeleitfähigkeit des Betonelements

    """

    def __init__(self, A_he, x_min, lambda_Bet):
        self.A_he = float(A_he)
        self.x_min = float(x_min)
        self.lambda_Bet = float(lambda_Bet)

    def __repr__(self):
        s = ('HeatingElement(A_he={self.A_he}, x_min={self.x_min}, \
             lambda_Bet={self.lambda_Bet}').format(self=self)
        return s

    def A_he(self):
        return self.A_he

    def x_min(self):
        return self.x_min

    def lambda_Bet(self):
        return self.lambda_Bet
