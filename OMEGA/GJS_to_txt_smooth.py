import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from scipy.ndimage import gaussian_filter

# ----- Mac -----
# sys.path.append("/Users/soham/Documents/Kinshock/Scripts")

# ----- Windows -----
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")

from OMEGA.GJS_analyze import GJS_analyze
from OMEGA.GJS_target_disp import interp_Prz
from h5_helpers.h5_helper import h5_to_dict


# ============================================================
# USER OPTIONS
# ============================================================

# "axial"  -> centerline profile copied across transverse x
# "full2d" -> original full 2D gas-jet profile
profile_mode = "full2d"

# If True and profile_mode == "axial", replace centerline
# with a manually chosen Gaussian. No fitting is performed.
fit = True

# Manual Gaussian parameters
gaussian_center_mm = 0.0
gaussian_sigma_mm = 5.0

# If True, use original centerline peak.
# Otherwise use gaussian_peak below.
use_original_peak = False
gaussian_peak = 0.20  # kg/m^3

# Optional Gaussian smoothing in GRID CELLS
# 0 = OFF
smooth_sigma = 0

# Pressure scaling
pressure_scale = 1300/1500


# ============================================================
# FILE PATH
# ============================================================

# ----- Mac -----
# baseDir = Path("/Users/soham/Documents/Kinshock")

# ----- Windows -----
baseDir = Path(r"\\profiles\Users$\sban\Documents\Kinshock")
fp = baseDir / "Kinshock-26A" / "GasJet" / "D-GJ-C-232_H2He_1500psi.h5"


# ============================================================
# READ GAS JET DATA
# ============================================================

r, z, q = GJS_analyze(fp, "density")
q = q * pressure_scale
P_rz, P_xyz = interp_Prz(q, r, z)

print("r range:", r.min(), r.max())
print("z range:", z.min(), z.max())


# ============================================================
# CREATE 2D PLANE FOR FLASH
# ============================================================

x = np.linspace(-10, 10, 1000)  # mm
y = np.linspace(-10, 10, 1000)  # mm
Z = 5  # distance from nozzle upward, mm

X, Y = np.meshgrid(x, y)
rho_interp = P_xyz(X, Y, Z)

print("rho_interp shape:", rho_interp.shape)
print("Original rho min:", np.nanmin(rho_interp))
print("Original rho max:", np.nanmax(rho_interp))


# ============================================================
# BUILD REQUESTED PROFILE
# ============================================================

profile_mode = profile_mode.lower()
ix_axis = np.abs(x).argmin()
rho_axis_original = rho_interp[:, ix_axis]

if profile_mode == "axial":

    print("Using AXIAL-ONLY gas profile")
    print("Using axis location x =", x[ix_axis], "mm")

    rho_axis_final = rho_axis_original.copy()
    output_tag = "axial"

    if fit:
        print()
        print("Using MANUAL GAUSSIAN axial profile")
        print("No curve fitting is being performed.")

        if use_original_peak:
            rho_peak = np.nanmax(rho_axis_original)
        else:
            rho_peak = gaussian_peak

        # Correct Gaussian: exponent must contain power 2, not 3
        rho_axis_gaussian = rho_peak * np.exp(-0.5 * ((y - gaussian_center_mm) / gaussian_sigma_mm) ** 4)

        rho_axis_final = rho_axis_gaussian.copy()
        output_tag = "axial_gaussian"

        print(f"Gaussian center = {gaussian_center_mm:.3f} mm")
        print(f"Gaussian sigma  = {gaussian_sigma_mm:.3f} mm")
        print(f"Gaussian FWHM   = {2.355 * gaussian_sigma_mm:.3f} mm")
        print(f"Gaussian peak   = {rho_peak:.6e} kg/m^3")

    rho_final = np.repeat(rho_axis_final[:, np.newaxis], len(x), axis=1)

elif profile_mode == "full2d":

    print("Using FULL 2D gas profile")
    rho_final = rho_interp.copy()
    output_tag = "full2d"

    if fit:

        print()
        print("Using MANUAL 2D GAUSSIAN gas profile")
        print("No curve fitting is being performed.")

        if use_original_peak:
            rho_peak = np.nanmax(rho_interp)
        else:
            rho_peak = gaussian_peak

        rho_final = rho_peak * np.exp(-0.5 * (((X - gaussian_center_mm)**4 + (Y - gaussian_center_mm)**4) / gaussian_sigma_mm**4))

        output_tag = "full2d_gaussian"

        print(f"Gaussian center = ({gaussian_center_mm:.3f}, {gaussian_center_mm:.3f}) mm")
        print(f"Gaussian sigma  = {gaussian_sigma_mm:.3f} mm")
        print(f"Gaussian FWHM   = {2.355 * gaussian_sigma_mm:.3f} mm")
        print(f"Gaussian peak   = {rho_peak:.6e} kg/m^3")

else:
    raise ValueError("profile_mode must be either 'axial' or 'full2d'")


# ============================================================
# OPTIONAL GAUSSIAN SMOOTHING
# ============================================================

if smooth_sigma > 0:
    print()
    print(f"Applying Gaussian smoothing with sigma = {smooth_sigma} grid cells")
    rho_final = gaussian_filter(rho_final, sigma=smooth_sigma, mode="nearest")
    output_tag += f"_smooth_sigma{smooth_sigma}"

else:
    print()
    print("Gaussian smoothing OFF")
    output_tag += "_nosmooth"


# ============================================================
# PRINT USEFUL INFORMATION
# ============================================================

dx_mm = x[1] - x[0]
dy_mm = y[1] - y[0]

print()
print("Grid spacing:")
print(f"dx = {dx_mm:.6f} mm")
print(f"dy = {dy_mm:.6f} mm")

if smooth_sigma > 0:
    print(f"Approximate smoothing length in x = {smooth_sigma * dx_mm:.6f} mm")
    print(f"Approximate smoothing length in y = {smooth_sigma * dy_mm:.6f} mm")

print()
print("Final rho min:", np.nanmin(rho_final))
print("Final rho max:", np.nanmax(rho_final))


# ============================================================
# WRITE FLASH TABLE
# ============================================================

outDir = baseDir / "Kinshock-26A" / "GasJet" / "FLASH_tables"
outDir.mkdir(parents=True, exist_ok=True)

outFile = outDir / f"M510mm_1300psi_10mmFoilOffset_full2d_smooth.txt"

# Convert mm -> cm.
X_cm = X / 10.0

# Shift y from [-10,10] mm -> [0,20] mm, then convert to cm.
Y_cm = (Y + 10.0) / 10.0

with open(outFile, "w") as f:
    f.write("# Gas jet interpolated density table\n")
    f.write(f"# Profile mode: {profile_mode}\n")
    f.write(f"# Manual Gaussian replacement: {fit}\n")

    if fit and profile_mode == "axial":
        f.write(f"# Gaussian center mm: {gaussian_center_mm}\n")
        f.write(f"# Gaussian sigma mm: {gaussian_sigma_mm}\n")
        f.write(f"# Gaussian peak kg/m3: {rho_peak}\n")

    f.write(f"# Gaussian smoothing sigma: {smooth_sigma} grid cells\n")
    f.write("# Columns: x_cm y_cm rho_g_cm3\n")
    f.write(f"# nx = {len(x)}, ny = {len(y)}\n")
    f.write("# x_min_cm x_max_cm y_min_cm y_max_cm\n")
    f.write(f"# {X_cm.min():.8e} {X_cm.max():.8e} {Y_cm.min():.8e} {Y_cm.max():.8e}\n")

    for j in range(len(y)):
        for i in range(len(x)):
            rho_gcc = rho_final[j, i] / 1000.0
            f.write(f"{X_cm[j, i]:.8e} {Y_cm[j, i]:.8e} {rho_gcc:.8e}\n")

print()
print("Wrote table to:")
print(outFile)


# ============================================================
# VERIFY OUTPUT
# ============================================================

fig, ax = plt.subplots(1, 3, figsize=(15, 4))

# Original 2D profile
im0 = ax[0].pcolormesh(X, Y, rho_interp / 1000.0, shading="auto")
ax[0].set_xlabel("x (mm)")
ax[0].set_ylabel("y (mm)")
ax[0].set_title("Original GJS profile")
fig.colorbar(im0, ax=ax[0], label=r"$\rho$ (g/cm$^3$)")

# Final profile
im1 = ax[1].pcolormesh(X, Y, rho_final / 1000.0, shading="auto")
ax[1].set_xlabel("x (mm)")
ax[1].set_ylabel("y (mm)")
ax[1].set_title(f"FLASH profile\n{output_tag}")
fig.colorbar(im1, ax=ax[1], label=r"$\rho$ (g/cm$^3$)")

# Axial lineout comparison
ax[2].plot(y, rho_axis_original / 1000.0, label="Original centerline", linewidth=2)
ax[2].plot(y, rho_final[:, ix_axis] / 1000.0, label="Final profile", linestyle="--", linewidth=2)

if profile_mode == "axial" and fit:
    ax[2].plot(y, rho_axis_gaussian / 1000.0, linestyle=":", linewidth=2, label="Gaussian guess")

ax[2].set_xlabel("y (mm)")
ax[2].set_ylabel(r"$\rho$ (g/cm$^3$)")
ax[2].set_title("Centerline density")
ax[2].legend()
ax[2].grid()

plt.tight_layout()
plt.show()