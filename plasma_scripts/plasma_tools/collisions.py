import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import astropy.units as u

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h

def mfp_ii(mi, mj, Zi, Zj, Ti, nj):
    """
    Calculates ion-ion mean free path for m1 colliding with m2 in a thermalized plasma

    Parameters
    ----------
    mi, mj : float
        Masses of the ion species in atomic mass units.
    Ti : float
        Temprature of the projectile species in KeV
    nj : float
        Density of the target species in 10^(20) cm^-3    
    Zi, Zj : float
        Charge states of species in e units

    Returns
    -------
    t_ij : float
        Ion-Ion mean free path in thermalized plasma in mm

    Examples
    --------
    >>> mfp_ii(1,1,1,1,1,1)
    0.22999999999999998
    """

    lnA = 10  # Coulomb Logarithm
    c = 1.15 # mm, constant factor from Hans's ICF workshop
    mr = (mi*mj)/(mi+mj) # reduced mass
    num = (mi*(Ti)**(2))
    den = ((mr)*nj*((Zi*Zj)**2)*lnA)
    mfp_ij = (num/den)*c
    return mfp_ij

def col_t_ii(mi = 1, mj = 1, Ti = 1, nj = 1, Zi = 1, Zj = 1):
    """
    Calculates ion-ion collision time for m1 colliding with m2 in a thermalized plasma

    Parameters
    ----------
    mi, mj : float
        Masses of the ion species in atomic mass units.
    Ti : float
        Temprature of the projectile species in KeV
    nj : float
        Density of the target species in 10^(20) cm^-3    
    Zi, Zj : float
        Charge states of species in e units

    Returns
    -------
    t_ij : float
        Ion-Ion collision time in s

    Examples
    --------
    >>> col_t_ij(1,1,1,1,1,1)
    6.593658098434409e-10
    """

    lnA = 10   # Coulomb Logarithm
    c1 = (3*(pi)**(0.5)/(2**(0.5)))*((4*pi*eps)**2)*(1/e**4)*(1/(4*pi)) # constants from the MIT OCW formula
    cn = 10**(-26)  # 10^20 cm^-3 to m^-3 in denominator
    ct = (1000*e)**(1.5) # 1 KeV to J
    cm = (mp)**(0.5) # mass factor 
    c = c1*cn*ct*cm # total constant conversion factor
    mr = (mi*mj)/(mi+mj) # reduced mass
    num = (mi*(Ti)**(1.5))
    den = ((mr)**(0.5)*nj*((Zi*Zj)**2)*lnA)
    t_ij = (num/den)*c
    return t_ij
