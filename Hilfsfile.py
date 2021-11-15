# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 14:22:53 2021

@author: jai_n

Plotting, Prototyping
"""
import math, sys
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from load_generator import *

#%% Vergleich der WÜK-Korrelationen von Bentz und VDI

# theta_l = 13.3e-6   kin. Vis. von Luft [m²/s] -- > bereits importiert
num = 50000

Re_vec = np.linspace(1e2, 1e6, num=num)
alpha_VDI_vec = np.zeros(num)
alpha_Bentz_vec = np.zeros(num)

for i,j in enumerate(Re_vec):
    alpha_Bentz_vec[i] = alpha_kon_Bentz(j * theta_l)
    # alpha_VDI_vec[i] = alpha_kon_VDI(10, j * theta_l, -5, 0)
    

fig = plt.figure()
ax = fig.add_subplot(111)
# ax.plot(Re_vec, alpha_VDI_vec)
ax.plot(Re_vec, alpha_Bentz_vec)
ax.set_xlabel('Re [-]')
ax.set_ylabel('alpha [W/m²K]')

ax.grid('on')

ax.legend(['alpha_VDI = fct(A = 10 m², u_inf, T_inf = -5 °C, T_o = 0 °C)', 'alpha_Bentz = fct(u_inf)'], 
          prop={'size': 10}, loc='upper left')

ax.set_xscale('log')

#%% Überlegungen: Verbesserung des Lösealgorithmus: F_Q

from sympy.plotting import plot
from load_generator import *
from load_generator_utils import *

h_NHN = 520
u_inf = 3
Theta_inf = -4
S_w = 1
A_he = 8
Theta_b = 10
R_th = 0.0051
Theta_surf = 0
S_w = 2
B = 3/8
Phi = 0.8
RR = 0.4
m_Rw_0 = 1
m_Rs_0 = 0
ssb = 0

# q = np.linspace(-10000, 10000, 1000)
# F_q = load(h_NHN, v, Theta_inf, S_w, A_he, Theta_b_0, R_th, Theta_surf, B, Phi, RR, m_Rw_0, m_Rs_0, ssb)
# plt.plot(q, F_q)
# plt.show()

 # 2.1) Pre-Processing
R_f = 1  # free-area ratio

# 2.2) Schmelz- und Verlustleistungsterme (parametriert) - Parameter Q
Q = sp.symbols('Q')  # thermische Leistung Q als Parameter definieren

# Q_Konvektion
Q_con = alpha_kon_Bentz(u_inf) * (Theta_b - Q * R_th - Theta_inf) * A_he

# Q_Strahlung
# Theta_surf_0 statt Theta_surf_, da sonst 4 Nullstellen
Q_rad = sigma * epsilon_surf('Beton') * ((Theta_b - Q * R_th + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

# Q_Verdunstung
''' Voraussetzungen: 
        - Theta_surf >= 0 °C
        - Oberfläche ist nass (Abfrage der Restwassermenge m_Rw)
'''
Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he
if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
    Q_eva = 0  

# Q_sensibel
# Q_sen = 0
Q_sen = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_b - Q * R_th - Theta_Schm)) * (3.6e6)**-1 * A_he

# Q_latent
# Q_lat = 0
Q_lat = rho_w * S_w * h_Ph_sl * (3.6e6)**-1 * A_he

# 2.3) stationäre Leistungbilanz (Erdboden + Oberfläche + Umgebung) & Auflösung nach Q
F_Q = sp.Eq(Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva) - Q, 0)
Q_sol = np.array(sp.solve(F_Q, Q))
print(type(Q_sol))

# print(len(Q_sol))
print(Q_sol)

F_Q_ = Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva) - Q
fig = plot(F_Q_, (Q, -3000, 3000))
fig.show()

#%% Überlegungen: Verbesserung des Lösealgorithmus: F_T


#%% Verdunstung

from load_generator import *

fig = plt.figure()
ax = fig

h_NHN = 520
u_inf = np.array([5, 3, 1])
Phi = np.array([0.51, 0.75])

T_inf_vec = np.linspace(-12, 0, num=50)
erg1 = np.zeros(len(T_inf_vec))
erg2 = np.zeros(len(T_inf_vec))
erg3 = np.zeros(len(T_inf_vec))
erg4 = np.zeros(len(T_inf_vec))
erg5 = np.zeros(len(T_inf_vec))
erg6 = np.zeros(len(T_inf_vec))
for i, j in enumerate(T_inf_vec):
    erg1[i] = load(h_NHN, u_inf[0], j, 2, 1, 0, 0, 0, 0.125, Phi[0], 2, 1, 0, 0)[0]
    erg2[i] = load(h_NHN, u_inf[0], j, 2, 1, 0, 0, 0, 0.125, Phi[1], 2, 1, 0, 0)[0]
    erg3[i] = load(h_NHN, u_inf[1], j, 2, 1, 0, 0, 0, 0.125, Phi[0], 2, 1, 0, 0)[0]
    erg4[i] = load(h_NHN, u_inf[1], j, 2, 1, 0, 0, 0, 0.125, Phi[1], 2, 1, 0, 0)[0]
    erg5[i] = load(h_NHN, u_inf[2], j, 2, 1, 0, 0, 0, 0.125, Phi[0], 2, 1, 0, 0)[0]
    erg6[i] = load(h_NHN, u_inf[2], j, 2, 1, 0, 0, 0, 0.125, Phi[1], 2, 1, 0, 0)[0]
    
plt.plot(T_inf_vec, erg1)
plt.plot(T_inf_vec, erg2)
plt.plot(T_inf_vec, erg3)
plt.plot(T_inf_vec, erg4)
plt.plot(T_inf_vec, erg5)
plt.plot(T_inf_vec, erg6)

ax.legend(['u=5, Phi=0.51', 'u=5, Phi=0.75', 'u=3, Phi=0.51', 'u=3, Phi=0.75', 'u=1, Phi=0.51', 'u=1, Phi=0.75',
         'u=1, Phi=0.75'])
plt.grid()
