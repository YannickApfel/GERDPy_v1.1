# -*- coding: utf-8 -*-
""" Autor: Yannick Apfel
"""
import numpy as np
from scipy.constants import pi


class Heatpipes(object):
    """
    Enthält alle Informationen über das Heatpipe-Layout im Bohrloch und dessen
    Wärmeleitfähigkeiten.
    Dies muss für alle Bohrlöcher identisch sein.

    Attribute
    ----------
    N:          int [-]
                Anzahl der Heatpipes im Bohrloch (kreisförmig angeordnet).
    r_b:        float [m]
                Radius der Erdwärmesondenbohrung
    r_w:        float [m]
                Radius der Wärmerohr-Mittelpunkte.
    r_pa:       float [m]
                Außenradius der Isolationsschicht.
    r_iso:      float [m]
                Innenradius der Isolationsschicht.
    r_pi:       float [m]
                Innenradius des Wärmerohrs.
    lambda_b:   float [W/mK]
                Wärmeleitfähigkeit des Hinterfüllmaterials
    lambda_iso: float [W/mK]
                Wärmeleitfähigkeit der Isolationsschicht.
    lambda_p:   float [W/mK]
                Wärmeleitfähigkeit der Heatpipe

    """

    def __init__(self, N, r_b, r_w, r_pa, r_iso, r_pi, lambda_b, lambda_iso, lambda_p):
        self.N = int(N)
        self.r_b = float(r_b)
        self.r_w = float(r_w)
        self.r_pa = float(r_pa)
        self.r_iso = float(r_iso)
        self.r_pi = float(r_pi)
        self.lambda_b = float(lambda_b)
        self.lambda_iso = float(lambda_iso)
        self.lambda_p = float(lambda_p)

    def __repr__(self):
        s = ('Heatpipes(N={self.N}, r_w={self.r_w}, r_pa={self.r_pa}, r_iso={self.r_iso},'
             ' r_pi={self.r_pi}, lambda_iso={self.lambda_iso},'
             ' lambda_p={self.lambda_p})').format(self=self)
        return s

    def N(self):
        return self.N

    def r_b(self):
        return self.r_b

    def r_w(self):
        return self.r_w

    def r_pa(self):
        return self.r_pa

    def r_iso(self):
        return self.r_iso

    def r_pi(self):
        return self.r_pi

    def lambda_b(self):
        return self.lambda_b

    def lambda_iso(self):
        return self.lambda_iso

    def lambda_p(self):
        return self.lambda_p

    def xy_mat(self):
        """
        Gibt eine Nx2 Matrix mit den Koordinaten der Heatpipe-Mittelpunkte
        im Bohrlochkoordinatensystem (Ursprung = Bohrlochmittelpunkt) zurück.
        Die N Heatpipes werden gleichmäßig auf einem Kreis mit
        Radius r_w verteilt.

        """

        xy_mat = np.zeros([self.N, 2])  # 1. Spalte: x, 2. Spalte: y

        for i in range(self.N):
            xy_mat[i, 0] = self.r_w * np.cos(2 * pi * i / self.N)
            xy_mat[i, 1] = self.r_w * np.sin(2 * pi * i / self.N)

        return xy_mat

    def visualize_hp_config(self):
        """
        Plot des Querschnitts eines Bohrlochs (visualisiert Heatpipe-Layout).

        Rückgabe
        -------
        fig : figure
            Figure object (matplotlib).

        """
        import matplotlib.pyplot as plt
        from matplotlib.ticker import AutoMinorLocator

        # Initialisierung
        LW = .5    # Line width
        FS = 12.    # Font size

        plt.rc('figure')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Umriss der Bohrlochwand
        borewall = plt.Circle((0., 0.), radius=self.r_b,
                              fill=False, linestyle='--', linewidth=LW)
        ax.add_patch(borewall)

        # Heat Pipes
        for i in range(self.N):
            # Koordinaten
            xy = self.xy_mat()
            (x, y) = (xy[i, 0], xy[i, 1])

            # Plot der Ummantellung und der Heatpipes
            hp_iso = plt.Circle((x, y), radius=self.r_pa,
                                fill=False, linestyle='-', linewidth=LW)
            hp_itself_outer = plt.Circle((x, y), radius=self.r_iso,
                                fill=False, linestyle='-', linewidth=LW)
            hp_itself_inner = plt.Circle((x, y), radius=self.r_pi,
                                fill=False, linestyle='-', linewidth=LW)
            ax.text(x, y, i + 1,
                    ha="center", va="center", size=FS)

            ax.add_patch(hp_iso)
            ax.add_patch(hp_itself_outer)
            ax.add_patch(hp_itself_inner)

        # Achsen
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title('Heatpipe Layout')
        plt.axis('equal')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        plt.tight_layout()

        return fig
