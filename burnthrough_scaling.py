import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")
from fitter_1D import fit_1d_curve


#burnthrough scaling v.s mass
# m = np.array([0.5, 1, 1.5, 2, 2.5, 3])*2.33
# bt = np.array([225.1, 350, 500.1, 625, 775, 925])

# def bt_scaling(mass,a, b):
#     return a*mass**(b)

# result = fit_1d_curve(m, bt, bt_scaling, title="Burnthrough scaling v.s mass", xlabel="Mass (g/cm^2)", ylabel="Burnthrough time (ps)", data_label="LILAC", fit_label=f"Fit = C*mass^{{0.85}}")

#burnthrough scaling v.s intensity
# I = np.arange(300, 1300, 100)
# print(I)
# bt = np.array([475, 425, 400, 375, 375, 350, 350, 350, 325, 325])

# def bt_scaling(intensity,a, b):
#     return a*intensity**(b)

# result = fit_1d_curve(I, bt, bt_scaling, title="Burnthrough scaling v.s intensity", xlabel="Intensity (TW/cm^2)", ylabel="Burnthrough time (ps)", data_label="LILAC", fit_label=f"Fit = C*I^{{-0.256}}")

# burnthrough scaling v.s. intensity for FLASH
I = np.array([2, 3, 4.16])

bt = np.array([3.3010359643e-01, 3.0170952909e-01, 2.8291859837e-01])

def bt_scaling(intensity,a, b):
    return a*intensity**(b)

result = fit_1d_curve(I, bt, bt_scaling, title="Burnthrough scaling v.s intensity", xlabel="Intensity (TW/cm^2)", ylabel="Burnthrough time (ps)", data_label="FLASH", fit_label=f"Fit = C*I^{{-0.44}}")

# burnthrough scaling v.s. thickness for FLASH
# m = np.array([0.5, 1, 1.5, 2, 2.5, 3])*2.33

# bt = np.array([ 1.4000154305e-01, 2.0001536052e-01,3.2000053413e-01, 4.0000698714e-01, 5.00e-01, 6.6001074187e-01])

# def bt_scaling(mass,a, b):
#     return a*mass**(b)

# result = fit_1d_curve(m, bt, bt_scaling, title="Burnthrough scaling v.s mass", xlabel="Mass (g/cm^2)", ylabel="Burnthrough time (ps)", data_label="FLASH", fit_label=f"Fit = C*mass^{{1.05}}")