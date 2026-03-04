import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants
from h5_import import h5_to_dict

def GJS_analyze(file,field):
    #--------- File path ---------- #
    baseDir = Path().resolve().parent
    fp = baseDir/"Kinshock"/"Kinshock-26A"/"GasJet"
    #file = "D-GJ-C-232_H2He_1500psi.h5"
    fp = fp/file

    #--------- Constants ---------- #
    mp = constants.proton_mass # mass of H(kg)
    eps = constants.epsilon_0 # epislon naught SI
    kb = constants.Boltzmann # Boltzmann constant SI
    e = constants.e # charge of electron(Coulomb)
    me = constants.electron_mass # mass electron(kg)
    mp = constants.proton_mass # mass of H(kg)
    pi = constants.pi # pi
    c = constants.c


    # ---------Load H5 file---------- #
    with h5py.File(fp, 'r') as f:
        full_data = h5_to_dict(f)
        #f.visititems(print_h5)

    # Density contours at specific z/r locations
    target = [5e-3, 7.25e-3, 8e-3, 10e-3, 15e-3]  # mm in meters
    arr = full_data['z']
    idx = np.zeros(len(target), dtype=int)
    for i in range(len(target)):
        idx[i] = np.abs(arr - target[i]).argmin()

    # ------ Data to visualize -------- #
    q = (full_data[field]) #/1000)*(80/1500)  # quantity to plot
    r = full_data['r']*1000 # mm
    z = full_data['z']*1000 # mm
    
    f1 = 9/10  # fraction of H
    f2 = 1/10  # fraction other gas
    A1 = 1 # atomic mass of H
    A2 = 14 # atomic mass of other gas
    Z = 3 # ionization state of the other gas
    n=q
    #n = q/(f1*(A1*mp) + f2*(A2*mp))
    #n = ((f1*n) + (Z*f2*n))/(10**6)  # convert ni/m^3 to ne/cm^3 
    #t = full_data['Density']  # temperature in K, from the ideal gas law
    #q = q/(kb*t*(10**6))
    #q = (3*q)/((5*mp)*(10**6))  # my attempt at going from kg/m3 -> /cm3 for H/He
    #q = (1*q)/((2.3*mp)*(10**6))  # my attempt at going from kg/m3 -> ni/cm3 for H/N
    #q = (1.4*q)/((2.3*mp)*(10**6))  # my attempt at going from kg/m3 -> ne/cm3 for H/N



    # ------ Plot -------- #
    # fig, ax = plt.subplots(1,2, figsize=(12,5))
    # im = ax[0].pcolormesh(r, z, n, cmap='viridis',shading='auto')
    # ax[0].axhline(y=5, color='red', linestyle='--', label='z = 5 mm')
    # ax[0].axhline(y=7.25, color='red', linestyle='--', label='z = 7.25 mm')
    # ax[0].axhline(y=8, color='red', linestyle='--', label='z = 8 mm')
    # ax[0].axhline(y=10, color='red', linestyle='--', label='z = 10 mm')

    # ax[0].legend()
    # ax[0].set_xlabel('r (mm)')
    # ax[0].set_ylabel('z (mm)')
    # plt.colorbar(im, label='n (g/cm^-3)', ax=ax[0])

    # for i in idx:
    #     ax[1].plot(r, n[i,:], linestyle='-',label=f'z = {full_data["r"][i]*1000} mm')

    # #ax[1].avhline(x=5, color='black', linestyle='--', label='foil location')
    # #ax[1].avhline(x=7.5, color='grey', linestyle='--')
    # #ax[1].avhline(x=2.5, color='grey', linestyle='--')
    # ax[1].set_xlabel('r (mm)')
    # ax[1].set_ylabel('n (g/cm^-3)')
    # ax[1].legend()

    # fig.suptitle("n for H2/N2 Gas Jet at 80 psi")
    # plt.show()
    return [r, z, q]


#GJS_analyze("D-GJ-C-232_H2N2_1500psi.h5", "density")