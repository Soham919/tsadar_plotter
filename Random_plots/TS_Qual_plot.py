import numpy as np
import matplotlib.pyplot as plt
import sys

# ----- Mac ----- #
# sys.path.append("/Users/soham/Documents/Kinshock/Scripts")

# ----- Windows ----- #
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")

from scipy.constants import c, m_e as me, m_p as mp
from plasma_tools import *
from lineout_tool import InteractiveLineout

# ============================================================
# FIGURE SETTINGS
# ============================================================

FIGSIZE = (6, 5)

FONT_SIZE = 16
LABEL_SIZE = 18
TITLE_SIZE = 18
TICK_SIZE = 16
LEGEND_SIZE = 14

AXIS_WIDTH = 1.2
TICK_WIDTH = 1.2
TICK_LENGTH = 5
LINE_WIDTH = 2.0

TOTAL_COLOR = "navy"
PROBE_COLOR = "mediumseagreen"

plt.rcParams.update({
    "font.family": "serif",
    "font.size": FONT_SIZE,
    "mathtext.fontset": "cm",

    "axes.labelsize": LABEL_SIZE,
    "axes.titlesize": TITLE_SIZE,
    "axes.linewidth": AXIS_WIDTH,

    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": TICK_WIDTH,
    "ytick.major.width": TICK_WIDTH,
    "xtick.major.size": TICK_LENGTH,
    "ytick.major.size": TICK_LENGTH,

    "legend.fontsize": LEGEND_SIZE,
    "lines.linewidth": LINE_WIDTH,
})

# ============================================================
# Probe geometry
# ============================================================

probe_wavelength_nm = 526.5
scattering_angle_deg = 60.0

probe_wavelength_m = probe_wavelength_nm * 1e-9
w0 = 2.0 * np.pi * c / probe_wavelength_m

scattering_angle_rad = np.deg2rad(scattering_angle_deg)
k = (4.0 * np.pi / probe_wavelength_m) * np.sin(scattering_angle_rad / 2.0)

# ============================================================
# Frequency / phase-velocity axis
# ============================================================

phase_velocity = np.linspace(-7.0e7, 7.0e7, 10000)

# ============================================================
# Plasma conditions
# ============================================================

ne_cm3 = 1e19

Te_eV = 100
TH_eV = 100

ud_e = 0.0
ud_H = 0.0

ZH = 0.1

# ============================================================
# Calculate Thomson spectrum
# ============================================================

electron = PlasmaSpecies(
    name="electron", Z=-1, m=me,
    T_eV=Te_eV, n_cm3=ne_cm3, ud=ud_e * 1e4
)

hydrogen = PlasmaSpecies(
    name="hydrogen", Z=ZH, m=mp,
    T_eV=TH_eV, n_cm3=ne_cm3 / ZH, ud=ud_H * 1e4
)

ions = [hydrogen]

Result = thomson_structure_factor(
    k=k,
    phase_velocity=phase_velocity,
    electron=electron,
    ions=ions
)

spectra = np.asarray(Result["S_total"])

# ============================================================
# Convert phase velocity -> wavelength
# ============================================================

lam = (2.0 * np.pi * c / (w0 + k * phase_velocity)) * 1e9  # nm

# Wavelength decreases with increasing frequency, so sort it
sort_idx = np.argsort(lam)
lam = lam[sort_idx]
spectra = spectra[sort_idx]

# ============================================================
# Normalize spectrum to unit area in wavelength
# ============================================================

spectra_norm = spectra / np.trapezoid(spectra, lam)

# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(figsize=FIGSIZE)

ax.plot(lam, spectra_norm, color=TOTAL_COLOR)

ax.set_xlabel(r"$\lambda$ (nm)")
ax.set_ylabel("Normalized spectral density")
ax.set_xlim(475, 576)

# Probe wavelength
ax.axvline(
    probe_wavelength_nm,
    color=PROBE_COLOR,
    linestyle="--",
    linewidth=LINE_WIDTH,
    alpha=0.9
)

# ============================================================
# Custom wavelength ticks
# ============================================================

xticks = np.arange(475, 576, 25)           # multiples of 25 nm
xticks = xticks[~np.isclose(xticks, 525)]  # remove 525 nm
xticks = np.sort(np.append(xticks, probe_wavelength_nm))

ax.set_xticks(xticks)

# Make 526.5 display explicitly
ax.set_xticklabels([
    "526.5" if np.isclose(x, probe_wavelength_nm) else f"{x:.0f}"
    for x in xticks
])

ax.tick_params(which="both", direction="in")

fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.96)

# ============================================================
# Save
# ============================================================

# fig.savefig("thomson_features.pdf")

plt.show()