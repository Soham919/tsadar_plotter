import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants
from Omega60_func import gj_dens, gj_n
import sys
# ----- Mac -----#
sys.path.append("/Users/soham/Documents/Kinshock/Scripts")
# ----- Windows -----#
#sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")
from h5_helpers.h5_helper import h5_to_dict

#--------- Constants ---------- #
mp = constants.proton_mass # mass of H(kg)
eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c

def GJS_analyze(fp,field):
    # --------- Load H5 file ---------- #
    with h5py.File(fp, 'r') as f:
        full_data = h5_to_dict(f)
        #f.visititems(print_h5)

    ## ------ Data to visualize -------- ##
    q = (full_data[field]) #/1000)*(80/1500)  # quantity to plot
    r = full_data['r']*1000 # mm
    z = full_data['z']*1000 # mm
    return [r, z, q]

# ---- Run for seeing plot ---- #
if __name__ == "__main__":
   # --------- Mac ---------- #
    #baseDir = Path("/Users/soham/Documents/Kinshock")
   # --------- Windows ---------- #
    baseDir = Path(r"\\profiles\Users$\sban\Documents\Kinshock")
    fp = baseDir/"Kinshock-26A"/"GasJet"
    file = "D-GJ-C-232_H2He_1500psi.h5"
    fp = fp/file

    r, z, q = GJS_analyze(fp, "density")
    # ------ Plot -------- #

    ## Density contours at specific z/r locations
    target = [5, 7.25, 8, 10, 15]  # mm in meters
    idx = np.zeros(len(target), dtype=int)
    for i in range(len(target)):
        idx[i] = np.abs(z - target[i]).argmin()

    f1 = 8/10  # fraction of H
    f2 = 2/10  # fraction other gas
    A1 = 1 # atomic mass of H
    A2 = 4 # atomic mass of other gas
    Z1 = 1 # ionization state of H
    Z2 = 2 # ionization state of the other gas
    n = q
    n = q/(f1*(A1*mp) + f2*(A2*mp))  # convert from mass density to average number density
    n = n/(10**6)  # convert to ni/cm^3 
    ne = ((Z1*f1*n) + (Z2*f2*n))  # no. of electrons per cm^3 
    #q2 = (n*(0.8)*A1*mp + n*(0.2)*A2*mp)/(10**3)  # convert back to mass density for
    #q = (3*q)/((5*mp)*(10**6))  # my attempt at going from kg/m3 -> /cm3 for H/He
    #q = (1*q)/((2.3*mp)*(10**6))  # my attempt at going from kg/m3 -> ni/cm3 for H/N
    #q = (1.4*q)/((2.3*mp)*(10**6))  # my attempt at going from kg/m3 -> ne/cm3 for H/N

    n_hansen = gj_dens(600,3,5,5)  # number density at z = 5 mm from Hansen's formula

    fig, ax = plt.subplots(1,2, figsize=(12,5))
    im = ax[0].pcolormesh(r, z, q, cmap='viridis',shading='auto')
    ax[0].axhline(y=5, color='red', linestyle='--', label='z = 5 mm')
    ax[0].axhline(y=7.25, color='red', linestyle='--', label='z = 7.25 mm')
    ax[0].axhline(y=8, color='red', linestyle='--', label='z = 8 mm')
    ax[0].axhline(y=10, color='red', linestyle='--', label='z = 10 mm')

    ax[0].legend()
    ax[0].set_xlabel('r (mm)')
    ax[0].set_ylabel('z (mm)')
    plt.colorbar(im, label='rho (kg/m^-3)', ax=ax[0])

    for i in idx:
        ax[1].plot(r, q[i,:], linestyle='-',label=f'z = {r[i]} mm')
    
    #ax[1].axhline(y=n_hansen, color='red', linestyle='--', label='z = 5 mm Hansen')  # expectation from Hansen
    #ax[1].axvline(x=5, color='black', linestyle='--', label='foil location')
    #ax[1].axvline(x=7.5, color='grey', linestyle='--')
    #ax[1].axvline(x=2.5, color='grey', linestyle='--')
    ax[1].set_xlabel('r (mm)')
    ax[1].set_ylabel('rho (kg/m^-3)')
    ax[1].legend()

    fig.suptitle("rho for H2(0.8)/He2(0.2) Gas Jet at 1500 psi")
    plt.show()