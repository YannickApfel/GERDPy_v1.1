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
import time as tim
from load_generator import *
from load_generator_utils import *

tic = tim.time()

h_NHN = 520
u_inf = 5
Theta_inf = -4
S_w = 1
A_he = 3
Theta_b_0 = 10
R_th = 0.001
Theta_surf_0 = 5
B = 3/8
Phi = 0.8
RR = 1
m_Rw_0 = 1
m_Rs_0 = 0
ssb = 0

con = True
rad = True
eva = True
sen = True
lat = True


# 2.1) Pre-Processing
R_f = 1  # free-area ratio

# Iterativer Lösungsalgorithmus:

step_refine = 0  # Hilfsvariable zur Verfeinerung der Schrittweite
step = 10  # doppelter Startwert als Iterationsschrittweite für Q
res = 0.01  # zulässiges Residuum für F_Q (Restfehler)

Theta_surf = 0  # Startwert für Q

error = abs(F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he))

while error > res:
    if F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) > 0:
        step_refine += 1
        while F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) > 0:
            Theta_surf -= (step / (2 * step_refine))  # Halbierung der Schrittweite für eine weitere Überschreitung/Unter- des Zielwerts
    elif F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) < 0:
        step_refine += 1
        while F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he) < 0:
            Theta_surf += (step / (2 * step_refine))  # Halbierung der Schrittweite für eine weitere Überschreitung/Unter- des Zielwerts

    error = abs(F_T(R_f, lat, S_w, Theta_surf, sen, Theta_inf, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi, B, A_he))

Theta_surf_sol = Theta_surf

print(error)
print(type(Theta_surf_sol))
print(Theta_surf_sol)

toc = tim.time()
print('Total simulation time: {} sec'.format(toc - tic))

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
