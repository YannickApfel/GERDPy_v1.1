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

#%% Plot F(q_str) = sigma * epsilon_surf('Beton') * ((T_b_0 - q * R_th + 273.15) ** 4 - T_MS(S_w, T_inf, B, Phi) ** 4)

from load_generator import T_MS

T_b_0 = 8
R_th = 0.1
S_w = 2
T_inf = -6
B = 3/8
Phi = 0.8

q_str = np.linspace(0, 100, 1000)
F_q_str = 5.67e-8 * 0.94 * ((T_b_0 - q_str * R_th + 273.15) ** 4 - T_MS(S_w, T_inf, B, Phi) ** 4)
plt.plot(q_str, F_q_str)
plt.show()

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
    erg1[i] = load(h_NHN, u_inf[0], j, 2, 1, 0, 0, 0, 0.125, Phi[0], 2, 1)[0]
    erg2[i] = load(h_NHN, u_inf[0], j, 2, 1, 0, 0, 0, 0.125, Phi[1], 2, 1)[0]
    erg3[i] = load(h_NHN, u_inf[1], j, 2, 1, 0, 0, 0, 0.125, Phi[0], 2, 1)[0]
    erg4[i] = load(h_NHN, u_inf[1], j, 2, 1, 0, 0, 0, 0.125, Phi[1], 2, 1)[0]
    erg5[i] = load(h_NHN, u_inf[2], j, 2, 1, 0, 0, 0, 0.125, Phi[0], 2, 1)[0]
    erg6[i] = load(h_NHN, u_inf[2], j, 2, 1, 0, 0, 0, 0.125, Phi[1], 2, 1)[0]
    
plt.plot(T_inf_vec, erg1)
plt.plot(T_inf_vec, erg2)
plt.plot(T_inf_vec, erg3)
plt.plot(T_inf_vec, erg4)
plt.plot(T_inf_vec, erg5)
plt.plot(T_inf_vec, erg6)

ax.legend(['u=5, Phi=0.51', 'u=5, Phi=0.75', 'u=3, Phi=0.51', 'u=3, Phi=0.75', 'u=1, Phi=0.51', 'u=1, Phi=0.75',
         'u=1, Phi=0.75'])
plt.grid()
