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
u_inf = 30
Theta_inf = -4
S_w = 1
A_he = 8
Theta_b_0 = 10
R_th = 0.051
Theta_surf_0 = 5
S_w = 0.1
B = 3/8
Phi = 0.8
RR = 0.4
m_Rw_0 = 1
m_Rs_0 = 0
ssb = 0

con = True
rad = True
eva = True
sen = True
lat = True


def F_Q(lat, S_w, A_he, Q, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi):

    F_Q = Q_lat(lat, S_w, A_he) \
        + Q_sen(Q, sen, S_w, Theta_inf, Theta_b_0, R_th, A_he) \
            + R_f \
                * (Q_con(Q, con, u_inf, Theta_b_0, R_th, Theta_inf, A_he) \
                    + Q_rad(Q, rad, Theta_b_0, R_th, S_w, Theta_inf, B, Phi, A_he) \
                        + Q_eva_pos(Q, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Theta_b_0, R_th, Phi, A_he)) \
                            - Q

    return F_Q

def Q_con(Q, con, u_inf, Theta_b_0, R_th, Theta_inf, A_he):  # Q_Konvektion
    Q_con = 0
    if con:
        Q_con = alpha_kon_Bentz(u_inf) * (Theta_b_0 - Q * R_th - Theta_inf) * A_he
        
    return Q_con


def Q_rad(Q, rad, Theta_b_0, R_th, S_w, Theta_inf, B, Phi, A_he):  # Q_Strahlung
    Q_rad = 0
    if rad:
        Q_rad = sigma * epsilon_surf('Beton') * ((Theta_b_0 - Q * R_th + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he
        
    return Q_rad


def Q_eva_pos(Q, eva, Theta_surf_0, m_Rw_0, Theta_inf, u_inf, h_NHN, Theta_b_0, R_th, Phi, A_he):  # Q_Verdunstung (nur positiv)
    Q_eva = 0
    if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
        Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_b_0 - Q * R_th, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he

    if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
        Q_eva = 0
    
    return Q_eva


def Q_sen(Q, sen, S_w, Theta_inf, Theta_b_0, R_th, A_he):  # Q_sensibel
    Q_sen = 0
    if sen:
        Q_sen = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_b_0 - Q * R_th - Theta_Schm)) * (3.6e6)**-1 * A_he
        
    return Q_sen


def Q_lat(lat, S_w, A_he):  # Q_latent
    Q_lat = 0
    if lat:
        Q_lat = rho_w * S_w * h_Ph_sl * (3.6e6)**-1 * A_he
        
    return Q_lat


# 2.1) Pre-Processing
R_f = 1  # free-area ratio

# # 2.2) Schmelz- und Verlustleistungsterme (parametriert) - Parameter Q
# Q = sp.symbols('Q')  # thermische Leistung Q als Parameter definieren

# # Q_Konvektion
# Q_con = 0
# if con:
#     Q_con = alpha_kon_Bentz(u_inf) * (Theta_b_0 - Q * R_th - Theta_inf) * A_he

# # Q_Strahlung
# Q_rad = 0
# if rad:  # Theta_surf_0 statt Theta_surf_, da sonst 4 Nullstellen
#     Q_rad = sigma * epsilon_surf('Beton') * ((Theta_b_0 - Q * R_th + 273.15) ** 4 - T_MS(S_w, Theta_inf, B, Phi) ** 4) * A_he

    
# # Q_Verdunstung
# Q_eva = 0
# if (eva and Theta_surf_0 >= 0 and m_Rw_0 > 0):
#     Q_eva = rho_l * beta_c(Theta_inf, u_inf, h_NHN) * (X_D_sat_surf(Theta_surf_0, h_NHN) - X_D_inf(Theta_inf, Phi, h_NHN)) * h_Ph_lg * A_he
# if Q_eva < 0:  # Verdunstungswärmestrom ist definitorisch positiv! <--> Kondensation wird vernachlässigt
#     Q_eva = 0  

# # Q_sensibel
# Q_sen = 0
# if sen:
#     Q_sen = rho_w * S_w * (c_p_s * (Theta_Schm - Theta_inf) + c_p_w * (Theta_b_0 - Q * R_th - Theta_Schm)) * (3.6e6)**-1 * A_he

# # Q_latent
# Q_lat = 0
# if lat:
#     Q_lat = rho_w * S_w * h_Ph_sl * (3.6e6)**-1 * A_he 

# # 2.3) stationäre Leistungbilanz (Erdboden + Oberfläche + Umgebung) & Auflösung nach Q
# F_Q = sp.Eq(Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva) - Q, 0)
# Q_sol = float(np.array(sp.solve(F_Q, Q))[0])  # Nullstelle Nr. 1


# Iterativer Lösungsalgorithmus:

step_refine = 0  # Hilfsvariable zur Verfeinerung der Schrittweite
step = 30  # doppelter Startwert für Iterationsschrittweite für Q_var
res = 0.5  # zulässiges Residuum für F_Q (Restfehler)

Q_var = 0  # Startwert für Leistung

error = abs(F_Q(lat, S_w, A_he, Q_var, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi))

while error > res:
    if F_Q(lat, S_w, A_he, Q_var, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi) > 0:
        step_refine += 1
        while F_Q(lat, S_w, A_he, Q_var, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi) > 0:
            Q_var += (step / (2 * step_refine))
    elif F_Q(lat, S_w, A_he, Q_var, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi) < 0:
        step_refine += 1
        while F_Q(lat, S_w, A_he, Q_var, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi) < 0:
            Q_var -= (step / (2 * step_refine))

    error = abs(F_Q(lat, S_w, A_he, Q_var, sen, Theta_inf, Theta_b_0, R_th, con, u_inf, rad, eva, Theta_surf_0, m_Rw_0, h_NHN, Phi))
    
Q_sol = Q_var

print(type(Q_sol))
print(Q_sol)

# Q = sp.symbols('Q')  # thermische Leistung Q als Parameter definieren
# F_Q = Q_lat + Q_sen + R_f * (Q_con + Q_rad + Q_eva) - Q
# fig = plot(F_Q, (Q, -10000, 10000))
# fig.show()


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
