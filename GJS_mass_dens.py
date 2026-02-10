import numpy
import matplotlib.pyplot as plt
from Omega60_func import gj_dens
from scipy import constants
mp = constants.proton_mass  # kg
n = gj_dens(696, 5, 10, 5)
f_H2 = 1/9
f_He = 1-f_H2
n_H = 2*f_H2*n
n_He = f_He*n
rho = (n_H*mp + n_He*4*mp)*1000  # kg/m3
print(rho)  
