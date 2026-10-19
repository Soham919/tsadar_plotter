import numpy as np
import opacplot2 as opp
import os
from pathlib import Path
# --------- Mac ----------- #
#flash_folder = Path("/Users/soham/Documents/Flash/2DCylindrical/Si3N4_cyl_5TW_1300psiGJS")

# ---- Windows ---- #
flash_folder = Path(r"C:\Simulation_data\FLASH\2D_Cylindrical\FV4.8\cvisc_0.3_carbuncle_debug")  # directory containing the FLASH plot files


amu = 1.66053906660e-24

A = 4.00
mpi = A * amu

fp = flash_folder / "he-imx-005.cn4"

op = opp.OpacIonmix(fp, mpi, man = True)

print(op.__dict__.keys())

print("Temperature grid:")
print(op.temps)

print("Density grid:")
print(op.dens)