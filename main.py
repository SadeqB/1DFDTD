""" One Dimensional FDTD """

import numpy as np
import matplotlib.pyplot as plt

# Constants
mu0 = 3*np.pi*1e-7
eps0 = 8.85e-12
c = 1/np.sqrt(mu0*eps0)

# Problem Defining Variables
wvlngt = 1e-6
k = 2*np.pi/wvlngt
omega = k*c
f = omega/2*np.pi
T = 1/f

L = 10*wvlngt

dz = wvlngt/20
dt = dz/c/2

Tf = 30*T

eps_r = 1
mu_r = 1
conductance = 1e3
# Field and material Allocations

Nz = int(L/dz + 1)
Nt = int(Tf/dt + 1)

eps = eps0*eps_r*np.ones((Nz,))
mu = mu0*mu_r*np.ones((Nz,))
sigma = conductance*np.ones((Nz,))

c1 = -dt/dz/mu
c2 = -dt/dz/eps/(1+sigma*dt/eps/2)
c3 = (1 - sigma*dt/eps/2)/(1 + sigma*dt/eps/2)

Ex = np.zeros((Nz, Nt))
Hy = np.zeros((Nz, Nt))

Jx = np.zeros((Nz, Nt))

# Initialization


# Time Evolution

for nt in range(1,Nt):
    Ex[int(Nz/2), nt-1] = np.sin(omega*nt*dt)
    for nz in range(Nz-1):

        Hy[nz, nt] = c1[nz]*(Ex[nz+1,nt-1] - Ex[nz, nt-1]) + Hy[nz, nt-1]
        Ex[nz, nt] = c2[nz]*(Hy[nz,nt] - Hy[nz-1, nt] + Jx[nz,nt]) + \
                     c3[nz]*Ex[nz, nt-1]

plt.plot(Ex[:,Nt-1])