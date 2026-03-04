import numpy as np
from matplotlib import pyplot as plt
from scipy import constants
from Plasma_func import *

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h

u_2019 = 0.75
epb = np.array([250,150,80]) # energy per beam in J
n = 12 # number of beams
E_total = epb*n + 120 # total energy in J, beam 47 is at 120 J
u_2026 = u_2019*((E_total/(2*2500))**(1/5))

t = np.linspace(0,10,50)
# 1st order correction to the 2019 design
fig, ax = plt.subplots()
for i in range(len(epb)):
    ax.plot(t, u_2026[i]*t, label=f'{u_2026[i]:.4f} mm/ns ({epb[i]}J, 2026)')

ax.plot(t, u_2019*t, label=f'{u_2019:.4f} mm/ns (2019)')
ax.axvline(x=8, color='red', linestyle='--', label=' 8ns')
ax.axhline(y=5.5, color='red', linestyle='--', label=' TSS pointing at 8 mm offset')

ax.set_xlabel('time (ns)')
ax.set_ylabel('Distance from Foil (mm)')

plt.minorticks_on()
ax.grid(which='minor', linestyle='--', linewidth=0.5)
ax.grid(which='major', linestyle='-', linewidth=0.8)
plt.legend()
plt.show()
