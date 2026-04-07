import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from Plasma_func import v_th, lam_db, omg_p

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h

def cs(m, Z, Te):
    m = m*mp
    Te = e*Te
    speed = ((5/3)*Z*Te/m)**(0.5)
    return speed

def cs_ij(a, n0,f1, f2, m1, m2, Z1, Z2, Te, Ti=0):
    """
    Calculates ion acoustic speed in a mixed thermalized plasma, using one of 3 ways -
    (i) My naive fluid approach
    (ii) E.A. Williams fluid formula from his paper (Physics of Fluids, 1995, https://doi.org/10.1063/1.871101)

    Parameters
    ----------
    a : int
        0 for my naive formula, 1 for E.A. Williams fluid formula, 2 for E.A. Williams kinetic formula
    f1, f2 : float
         Fractions of the ion species.
    m1, m2 : float
        Masses of the ion species in atomic mass units.
    Z1, Z2 : float
        Charge states of species in e units
    Te : float
        Electron temprature in eV
   
    Returns
    -------
    cs_ij : float
        Ion acoustic wave speed(sound speed) in a binary plasma in the cold limit Ti << Te
    """
    # n1 = n1*10**6
    # n2 = n2*10**6
    # n0 = Z1*n1 + Z2*n2
    # f1 = Z1*n1/n0
    # f2 = Z2*n2/n0
    ldb = lam_db(Te, n0)
    k = np.linspace(0.1,10,100)/ldb
    q = 1 + (k**2)*(ldb**2)
    
    if a == 0:
        # My naive approach
        cs1 = cs(m1,Z1,Te)**2  # This function takes Te in eV itself, so no need to convert
        cs2 = cs(m2,Z2,Te)**2
        cs_ij = (f1*Z1*cs1 + f2*Z2*cs2)**(0.5)

    elif a == 1:
        # Implement E.A. Williams fluid formula
        Te = e*Te # This function takes Te in eV, but we need it in Joules for this formula
        Z_avg = Z1*f1 + Z2*f2
        Z2_A = np.average([Z1**2/m1, Z2**2/m2], weights=[f1,f2])
        y1 = 3
        y2 = 3
        A = (y1/m1 + y2/m2)*(Ti/mp)**2 + (Z2_A*Te)/(Z_avg*mp*q)
        B1 = ((y1/m1 - y2/m2)**2)*(Ti/mp)**2
        B2 = 2*(y1/m1 - y2/m2)*(f1*(Z1**2)/m1 - f2*(Z2**2)/m2)*(Ti*Te/(mp**2*Z_avg*q))
        B3 = (Z2_A*Te/(Z_avg*mp*q))**2
        B = B1 + B2 + B3
        cs_ij_fast = (0.5*(A + B**(0.5)))**(0.5)
        cs_ij_slow = (0.5*(A - B**(0.5)))**(0.5)
        cs_ij = [cs_ij_fast, cs_ij_slow]
        pass

    return cs_ij

if __name__ == "__main__":
    a = 1
    f1 = 0.99
    f2 = 0.01
    m1 = 1
    m2 = 131
    Z1 = 1
    Z2 = 40
    Te = 1000
    n0 = 10**20
    ldb = lam_db(Te, ne=n0)
    k = np.linspace(0.1,10,100)/ldb

    cs_ij_H_Xe_fast, cs_ij_H_Xe_slow = cs_ij(a, n0=n0, f1=f1, f2=f2, m1=m1, m2=m2, Z1=Z1, Z2=Z2, Te=Te)
    
    fig, ax = plt.subplots()
    ax.plot(k*ldb, cs_ij_H_Xe_fast, label='H-Xe fast mode')
    ax.plot(k*ldb, cs_ij_H_Xe_slow, label='H-Xe slow mode')
    ax.set_xlabel(r'$k \lambda_D$')
    ax.set_ylabel('Ion Acoustic Speed (m/s)')
    ax.legend()
    plt.show()