"""
test explicit and implicit radial FD heat flow eq.

"""

__author__ = 'elco'

import pdb
import numpy as np
import matplotlib.pyplot as pl
import lib.pyBHTlib as pyBHTlib

# total length
length = 5.0

# cell radius
dr = 0.01

# borehole radius
br = 0.5

# borehole temeprature
T_borehole = 72.0

#
T_formation = 104.0

density = 2650.0

#
circ_time = 5.0 * 60.0 * 60.0
recovery_time = 20.0 * 60.0 * 60.0
dt = 5.0

# calculate number of timesteps
nt_circulation = int(circ_time / dt)
nt_recovery = int(recovery_time / dt)

# set up arrays
r = np.arange(0, length + dr, dr)
r_mid = (r[1:] + r[:-1]) / 2.0
T = np.ones_like(r) * T_formation
T_1D = T.copy()
rho = np.ones_like(r) * density
cp = np.ones_like(r) * 900.0
K = np.ones(r.shape[0]-1) * 2.5
q = np.zeros(r.shape[0] + 1)
nr = len(r)

# set up 2D arrays
T_2D = np.ones((nr, nr)) * T_formation
rho_2D = np.ones_like(T_2D) * density
cp_2D = np.ones_like(T_2D) * 900.0
Kh = np.ones((T_2D.shape[0]-1, T_2D.shape[0])) * 2.5
Kv = np.ones((T_2D.shape[0], T_2D.shape[0]-1)) * 2.5
qh = np.zeros((T_2D.shape[0] + 1, T_2D.shape[0]))
qv = np.zeros((T_2D.shape[0], T_2D.shape[0] + 1))

rx = np.resize(r, [nr, nr])
ry = np.resize(r, [nr, nr]).T
r_2D = np.sqrt(rx**2 + ry**2)

# set fixed T at r=0
borehole_index = r < (br + dr)
borehole_index_2D = r_2D < (br + dr)

# set initial temp in borehole
T[borehole_index] = T_borehole
T_1D[borehole_index] = T_borehole
T_2D[borehole_index_2D] = T_borehole

Tt = T.copy()

print 'circulation:'

for i in xrange(nt_circulation):

    T = pyBHTlib.radial_explicit_heat_flow(T, q, r, r_mid, dr, K, rho, cp, dt)
    #Tt = pyBHTlib.radial_explicit_heat_flow_test(Tt, q, r, r_mid, dr, K, rho, cp, dt)
    T_1D = pyBHTlib.explicit_heat_flow_1d(T_1D, q, dr, K, rho, cp, dt)
    T_2D = pyBHTlib.explicit_heat_flow_2d(T_2D, qh, qv, dr, dr,
                                          Kh, Kv, rho_2D, cp_2D, dt)
    # impose fixed T in borehole
    T[borehole_index] = T_borehole
    T_1D[borehole_index] = T_borehole
    T_2D[borehole_index_2D] = T_borehole

    if i / 100 == i / 100.0:
        print i, T[borehole_index].mean(), np.max(T - T_1D)


xmax = 1.0
ind = r+dr < xmax

fig = pl.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(r[ind], T_2D[0, :][ind], color='blue', lw=1.0, label='2D')
ax.plot(r[ind], T_1D[ind], color='red', lw=1.0, label='1D')
ax.plot(r[ind], T[ind], color='black', lw=1.0, label='radial')

ax.legend(loc='lower right')

fig.savefig('fig/test_radial_eq.png')


print 'recovery:'

for i in xrange(nt_recovery):

    T = pyBHTlib.radial_explicit_heat_flow(T, q, r, r_mid, dr, K, rho, cp, dt)
    #T_1D = pyBHTlib.explicit_heat_flow_1D(T_1D, q, dr, K, rho, cp, dt)

    if i / 1000 == i / 1000.0:
        print i, T[borehole_index].mean(), np.max(T - T_1D)


print 'done'
