import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c


def v_th(T,m):
    """
    Calculates the thermal velocity of particle species

    Parameters
    ----------
    T : float
        Temprature of the species in eV
    m : float
        mass of the species in kg
    
    Returns
    -------
    v_th : float
        thermal velocity of the particle species in m/s
    """
    T = T*e
    vt = (T/m)**(0.5)
    return vt

def lam_db(Te, ne):
    """
    Calculates the Debye length for a plasma

    Parameters
    ----------
    Te : float
        Temprature of the electrons in eV
    ne : float
        Electron density of the plasma in m^-3  

    Returns
    -------
    lam : float
        Debye length in m

    Examples
    --------
    >>> print(lam_db(10**18,1000))
    0.00023508188704793946
    """
    kT = Te*e
    omg = omg_p(ne)
    vt = (kT/me)**(0.5)
    lam = vt/omg
    return lam

def omg_p(ne):
    """
    Calculates the electron plasma frequency

    Parameters
    ----------
    ne : float
        Electron density of the plasma in m^-3  

    Returns
    -------
    omg : float
        plasma frequency in s^-1

    Examples
    --------
    >>> print(omg_p(10**18))
    56414602254.29499
    """
    omg = ((ne*(e**2))/(me*eps))**(0.5)
    return omg



def col_t_ii(mi, mj, Ti, nj, Zi, Zj):
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

def cs(m, Z, Te):
    m = m*mp
    Te = e*Te
    speed = ((5/3)*Z*Te/m)**(0.5)
    return speed

def nc(lam):
    """
    Calculates crirtical density of a laser with wavelength lam
    Parameters
    ----------
    lam: float
        Laser wavelength in nm
   
    Returns
    -------
    nc : float
        Plasma critical density in cm^-3

    Examples
    --------
    >>> nc(351)
    
    """
    lam = lam*(10**(-9))
    w = 2*pi*c/lam
    nc = (me*w**(2)*eps)/(e**(2))
    return nc/10**(6)

def cs_ij(n1, n2, m1, m2, Z1, Z2, Te):
    """
    Calculates ion acoustic speed in a mixed thermalized plasma
    Parameters
    ----------
    n1, n2 : float
         Density of the ion species in 10^20 cm^-3.
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

    Examples
    --------
    >>> cs_ij(1,1,1,1,1,1)
    0.22999999999999998
    """
    n1 = n1*10**6
    n2 = n2*10**6
    n0 = Z1*n1 + Z2*n2
    f1 = Z1*n1/n0
    f2 = Z2*n2/n0
    cs1 = cs(m1,Z1,Te)**2
    cs2 = cs(m2,Z2,Te)**2
    cs_ij = (f1*cs1 + f2*cs2)**(0.5)
    return cs_ij

def f_rho(f1, Z1, Z2, ne):
    """
    Calculates actual ion number density based on the electron density and the fractions provided from the Thomson analysis
    This is for a binary plasma
    Parameters
    ----------
    f1 : float
         Fraction of the ion species 
    Z1, Z2 : float
        Charge states of species in e units
    ne : float
        Electron density in whatever units you want
   
    Returns
    -------
    [n1,n2] : float
        Ion number densities of first and second species
    """
    f2 = 1-f1
    ni = ne/(Z1*f1 + Z2*f2)
    n1 = f1*ni
    n2 = f2*ni
    return [n1,n2]

def RH(M,rho1,u1,T1):
    """
    Calculates the Rankine-Hugoniot jump factors in a strong shock
    Parameters
    ----------
    M : float
        Mach number of the shock
    rho1, u1, T1 : float
        The density, velocity and temperature of the unshocked region in any units
   
    Returns
    -------
    [rho2,u2,T2] : float
        The density, velocity and temperature of the shocked region
    """
    y = 5/3
    v_fact = ((y-1)/(y+1)) + (2/((M**2)*(y+1)))
    u2 = v_fact*u1

    rho_fact = 1/v_fact
    rho2 = rho_fact*rho1

    P_fact = ((y+1)*rho_fact - (y-1))/((y+1)-((y-1)*rho_fact))
    T2 = (P_fact/rho_fact)*T1
    #y = (g+1)*x**2/((g-1)*(x**2)+2)
    return [rho2,u2,T2]

def blast_wave_r_v(E,rho1,t):
    """
    Calculates the ideal Sedov-Taylor expansion of a blast wave without any radiation losses
    Parameters
    ----------
    E : float
        Initial energy in the blast wave
    rho1 : float
        The density of the ambient medium which the blast wave is pushing against
    t : float
        Time after initial explosion
    Returns
    -------
    [r_t,v_t] : float
        The radius and velocity of the blast wave
    """
    eps = 1  # constant which holds significance, but I dunno how to calculate
    r_t = eps*((E*(t**2))/rho1)**(1/5)
    v_t = (2/5)*eps*(E/(rho1*(t**3)))**(1/5)
    return [r_t,v_t]

def power_calc(n,E,t,A):
    """
    Calculates the total power from laser beams on a certain surface area
    Parameters
    ----------
    n : float
        No. of beams
    E : float
        Energy in each beam (J)
    t : float
        Pulse duration (ps)
    A : float
        Area over which the beams are focused (cm^2)
    Returns
    -------
    P : float
        The power in TW/cm2
    """
    t = t*1e-12 # ps to s
    A = A*1e-4 # cm2 to m2
    P = (n*E)/ (t*A) # W/m2
    P = P*1e-16 # TW/cm2
    return P

#plt.plot(x,RH_T(x))
#plt.show()

#print("Speed of sound in H+/He2+ plasma = ",cs_ij(0.067,0.033,1,4,1,2,50))
#print("Speed of sound in H+ plasma = ",cs(1,1,50))

# c = np.linspace(0.01,0.99,50)
# f = 1-c
# df = 0.2
# fig1, ax1 = plt.subplots()
# ax1.plot(c,cs_ij(df*2*c,df*2*f,1,14,1,2,55),label =r"$N^{2+}$")
# ax1.plot(c,cs_ij(df*2*c,df*2*f,1,14,1,3,55),label =r"$N^{3+}$")
# ax1.plot(c,cs_ij(df*2*c,df*2*f,1,14,1,4,55),label =r"$N^{4+}$")
# ax1.set_xlabel(r"$H_{2}$(%)")
# ax1.set_ylabel(r"$cs_{H^{+}/N^{Z+}}$")
# ax1.set_title(r"Sound speed for $H_{2}/N_{2}$ gas mixture")
# ax1.legend()
# plt.show()

#print(((mp/me)**(0.5))*mfp_hh)


## ION DENSITIES FROM THE ELECTRON DENSITY
# fi, ax = plt.subplots()
# ax.plot(hconc,np.squeeze(f_rho(hconc,1,2,5)[0]),label="H")
# ax.plot(hconc,np.squeeze(f_rho(hconc,1,2,5)[1]),label="He")
# ax.set_xlabel("H fraction")
# ax.set_ylabel(r"$n_{i}$")
# ax.set_title(r"H and He number densities for $n_{e} = 5e19 cm^{-3}$")
# plt.legend()
# plt.show()

# T = np.linspace(0.01,0.5,50)
# fig, ax = plt.subplots()
# ax.plot(T,mfp_ii(1,1,1,1,T,0.3))
# ax.set_xlabel("T(KeV)")
# ax.set_ylabel(r"$MFP(\mu m)$")
# ax.set_title("H+/H+ MFP in 86801")
# plt.show()

## RANKINE HUGONIOT CHECK
# M = np.linspace(1,30,100)
# rho1 = 1
# u1 = 1
# P1 = 1
# [rho2,u2,P2] = RH(M,rho1,u1,P1)
# fig, ax = plt.subplots(3,1)
# ax[0].plot(M,rho2)
# ax[1].plot(M,u2)
# ax[2].plot(M,P2)
# plt.show()

# Sedov blast wave expansion
# t = np.linspace(1e-11,1e-8,100)
# E = 2500 # energy of the laser was 2.5 KJ
# E2 = 2000
# E3 = 1000
# rho1 = 2e19*mp*(10**6)
# [r,v] = blast_wave_r_v(E,rho1,t)
# fig, ax = plt.subplots(2,1)
# ax[0].plot(t*1e9,r/1000,label=r"$r ~ (\frac{Et^{2}}{\rho_{1}})^{1/5}$")

# ax[0].set_title('Radius')
# ax[0].set_xlabel(r"$t (ns)$")
# ax[0].set_xlabel(r"$r (mm)$")
# plt.legend()

# ax[1].plot(v,t*1e9,label=r"$v \sim (\frac{E}{\rho_{1}t^{3}})^{1/5}$")
# ax[1].set_title('velocity')
# ax[1].set_xlabel(r"$t (ns)$")
# ax[1].set_xlabel(r"$v (mm)$")

# plt.legend()
# plt.show()

A = pi*((0.003)**2)  # area in cm2
print(power_calc(1,3,0.03,A))

a = cs_ij(0.413431,0.103357,1,4,1,2,50)
b = cs(1,1,50)
#print(a)
#print(b)
