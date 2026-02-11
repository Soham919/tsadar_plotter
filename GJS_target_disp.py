import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants
from h5_import import h5_to_dict
from scipy.interpolate import RegularGridInterpolator

def interp_Prz(P_rz, r_grid, z_grid,
                          method="linear",
                          bounds_error=False,
                          fill_value=np.nan):
    """
    P_rz: 2D array sampled on (z_grid, r_grid).
          Expected shape = (len(z_grid), len(r_grid)).
          If you have (len(r_grid), len(z_grid)), it auto-transposes.
    r_grid: 1D increasing array of r values
    z_grid: 1D increasing array of z values

    Returns: function P_of_xyz(x,y,z) that interpolates P at (r=sqrt(x^2+y^2), z)
    """
    P_rz = np.asarray(P_rz)
    r_grid = np.asarray(r_grid)
    z_grid = np.asarray(z_grid)

    # fix common orientation mismatch
    if P_rz.shape == (len(r_grid), len(z_grid)):
        P_rz = P_rz.T
    if P_rz.shape != (len(z_grid), len(r_grid)):
        raise ValueError(
            f"P_rz must be shape (len(z_grid), len(r_grid)) = {(len(z_grid), len(r_grid))}, got {P_rz.shape}"
        )

    interp = RegularGridInterpolator(
        (z_grid, r_grid), P_rz,
        method=method,
        bounds_error=bounds_error,
        fill_value=fill_value
    )

    def P_of_zr(z, r):
        z = np.asarray(z)
        r = np.asarray(r)

        z, r = np.broadcast_arrays(z, r)
        pts = np.stack([z, r], axis=-1)   # (..., 2) as (z, r)
        return interp(pts)

    return P_of_zr

def force_on_foil(P_of_zr, x0, y1, y2, z1, z2, nz=100, ny=200):
    z = np.linspace(z1, z2, nz)
    y = np.linspace(y1, y2, ny)
    Z, Y = np.meshgrid(z, y, indexing="xy")

    R = np.sqrt(x0**2 + Y**2)
    

    # RegularGridInterpolator wants points (..., 2) in (z, r) order
    P = P_of_zr(Z,R)  # shape (nz, ny)

    # integrate dy then dz
    Fy = np.trapezoid(P, y, axis=0)   # integrate over y
    F  = np.trapezoid(Fy, z, axis=0)  # integrate over z
    return F

def GJS_target_disp(F,t,m):
    a = F/m
    d = 0.5*a*(t**2)
    return d

###----------------------------------------------------------------------------------------------------###

#---------File path---------- #
baseDir = Path().resolve().parent
fp = baseDir/"Kinshock-26A"/"GasJet"
file = "D-GJ-C-232_H2He_1500psi.h5"
fp = fp/file

# ---------Load H5 file---------- #
with h5py.File(fp, 'r') as f:
    full_data = h5_to_dict(f)
    #f.visititems(print_h5)

#---------- Interpolate Pressure -----------#
r = full_data['r']*1000 # mm
z = full_data['z']*1000 # mm
p = full_data['total-pressure'] # Pa

P_zr = interp_Prz(p,r,z)

#---------- Calculate Force and Displacement -----------#
F_foil = force_on_foil(P_zr, 15, -2.5, 2.5, 2.5, 7.5)/(10**6) # converting from mm^2 to m^2  

#P = 2e5  # Pressure in Pa for upper threshold
t = (160)*(10**(-6)) # time in seconds
A1 = 4e-6   # Area of the Si3N4 window in m^2
l1 = 10**(-6) # thickness of the Si3N4 window in m
A2 = (25-4)*(10**(-6))   # Area of the Si frame in m^2
l2 = 200*(10**(-6)) # thickness of the Si frame in m
rho1 = (3.2)*(10**(3)) # mass density of Si3N4 in kg/m^3
rho2 = (2.33)*(10**(3)) # mass density of Si in kg/m^3

m1 = rho1*A1*l1 # mass of window in kg
m2 = A2*l2*rho2 # mass of frame in kg

m = m1 + m2 # ~ total mass of target in kg
print(f"\nMass of target = {m*1000} g")

print(f"Average pressure on target = {F_foil/(A1+A2)} Pa")

print(f"\n\nTarget is displaced by {GJS_target_disp(F_foil,t,m)*1000} mm")