import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi


def gj_n(p):
    p = p*(6894.76)
    V = 150*(10**(-6))  # cylinder volume in m^3
    T = 273 
    rho = p/(kb*T) # number density in m^-3
    return rho/(10**6)

def gj_dens(p : float,M : float,D : float,x : float):
    """
    Calculates the gas number density from an LLE gas jet at distance x away from the nozzle exit.
    The calculations are taken from A. Hansen, Rev. Sci. Instrum. 89, 10C103 (2018), https://doi.org/10.1063/1.5036645

    Parameters
    ----------
    p : float
        Pressure in cylinder (psi)
    M : float
        Mach number of LLE nozzle
    D : float
        Exit diameter of LLE nozzle (mm)
    x : float
        Distance from nozzle (mm)
    
    Returns
    -------
    rho_x : float
        gas number density at distance x away from the nozzle exit (cm^-3)
    """
    y = 5/3      # adiabatic constant for monoatomic gas
    rho_c = gj_n(p)   # Number dens inside cylinder
    rho_ex = rho_c/(1+((y-1)/2)*(M**2))**(1/(y-1)) # Number density at nozzle exit
    R = D + 2*(x/M)   # fwhm of gas at distance x away from nozzle 
    rho_x = rho_ex*((D/R)**2)
    return rho_x

# y = 2
# M = 2
# rho_c = 1
# rho_ex = rho_c/(1+((y-1)/2)*(M**2))**(1/(y-1))
# print(rho_ex)

## PLot for density of different gas jet nozzles 20mm and 5mm comparison
# fig, ax = plt.subplots()
# x = np.linspace(0,30,100)
# ax.plot(x,1000*mp*gj_dens(1500,8,20,x),color='MediumSeaGreen',linewidth=2,label='M8, 20 mm')
# ax.plot(x,1000*mp*gj_dens(1500,5,10,x),color='orange',linewidth=2,label='M5, 10 mm')
# ax.plot(x,1000*mp*gj_dens(1500,5,5,x),color='orange', linestyle = ':',linewidth=2,label='M5, 5 mm')
# ax.plot(x,1000*mp*gj_dens(1500,5,2,x),color='orange', linestyle = '--',linewidth=2,label='M5, 2 mm')
# ax.plot(x,1000*mp*gj_dens(1500,6.5,15,x),color='blue',linewidth=2,label='M6.5, 15 mm')

# ax.set_xlim([0,max(x)])
# ax.set_ylabel(r"$\rho (g/cm^3)$")
# ax.set_xlabel("Distance from nozzle exit (mm)")

# plt.legend()
# plt.show()

## Plot for comparing density from 2019 and 2026
fig, ax = plt.subplots()
x = np.linspace(0,15,100)
ax.plot(x,1000*2*mp*gj_dens(696,5,10,x),color='MediumSeaGreen',linewidth=2,label='2019')
ax.plot(x,1000*2*mp*gj_dens(1500,5,10,x),color='orange',linewidth=2,label='2026')

ax.set_xlim([0,max(x)])
ax.set_ylabel(r"$\rho (g/cm^3)$")
ax.set_xlabel("Distance from nozzle exit (mm)")

plt.legend()
plt.show()