import numpy as np
import matplotlib.pyplot as plt
import sys

# ----- Mac ----- #
# sys.path.append("/Users/soham/Documents/Kinshock/Scripts")

# ----- Windows ----- #
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")

from scipy.constants import c, m_e as me, m_p as mp
from plasma_tools import *

# ============================================================
# FIGURE SETTINGS
# ============================================================

FIGSIZE = (12, 9)

FONT_SIZE = 16
LABEL_SIZE = 18
TITLE_SIZE = 18
TICK_SIZE = 15
LEGEND_SIZE = 14

AXIS_WIDTH = 1.2
TICK_WIDTH = 1.2
TICK_LENGTH = 5
LINE_WIDTH = 2.0

TOTAL_COLOR = "navy"
H_COLOR = "orangered"
HE_COLOR = "mediumseagreen"
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
# Common phase-velocity axis
# ============================================================

phase_velocity = np.linspace(-7.0e7, 7.0e7, 10000)

# Convert phase velocity -> wavelength
lam = (2.0 * np.pi * c / (w0 + k * phase_velocity)) * 1e9

# ============================================================
# PANEL (a): Full Thomson spectrum
# ============================================================

ne_dummy = 1e19
Te_dummy = 100
TH_dummy = 100
ZH_dummy = 0.1

electron_dummy = PlasmaSpecies(
    name="electron", Z=-1, m=me,
    T_eV=Te_dummy, n_cm3=ne_dummy, ud=0.0
)

hydrogen_dummy = PlasmaSpecies(
    name="hydrogen", Z=ZH_dummy, m=mp,
    T_eV=TH_dummy, n_cm3=ne_dummy / ZH_dummy, ud=0.0
)

result_dummy = thomson_structure_factor(
    k=k,
    phase_velocity=phase_velocity,
    electron=electron_dummy,
    ions=[hydrogen_dummy]
)

spectra_dummy = np.asarray(result_dummy["S_total"])

# Sort by wavelength
sort_full = np.argsort(lam)
lam_full = lam[sort_full]
spectra_dummy = spectra_dummy[sort_full]

# Normalize full spectrum to unit area
spectra_dummy_norm = spectra_dummy / np.trapezoid(spectra_dummy, lam_full)

# ============================================================
# PANELS (b-d): H/He composition scan
# ============================================================

ne_cm3 = 1e20

Te_eV = 500
TH_eV = 300
THe_eV = 300

ZH = 1.0
ZHe = 2.0

ud_e = 0.0
ud_H = 0.0
ud_He = 0.0

fH = np.linspace(0.01, 0.99, 500)
results = []

for f in fH:

    Zbar = ZH * f + ZHe * (1.0 - f)

    electron = PlasmaSpecies(
        name="electron", Z=-1, m=me,
        T_eV=Te_eV, n_cm3=ne_cm3, ud=ud_e * 1e4
    )

    hydrogen = PlasmaSpecies(
        name="hydrogen", Z=ZH, m=mp,
        T_eV=TH_eV, n_cm3=f * (ne_cm3 / Zbar), ud=ud_H * 1e4
    )

    helium = PlasmaSpecies(
        name="helium", Z=ZHe, m=4.0 * mp,
        T_eV=THe_eV, n_cm3=(1.0 - f) * (ne_cm3 / Zbar), ud=ud_He * 1e4
    )

    results.append(
        thomson_structure_factor(
            k=k,
            phase_velocity=phase_velocity,
            electron=electron,
            ions=[hydrogen, helium]
        )
    )

spectra = np.array([result["S_total"] for result in results])
S_H = np.array([result["S_ions"]["hydrogen"] for result in results])
S_He = np.array([result["S_ions"]["helium"] for result in results])

# ============================================================
# Ion-feature wavelength window
# ============================================================

ion_mask = np.abs(phase_velocity) <= 1.0e6

lam_ion = lam[ion_mask]
spectra_ion = spectra[:, ion_mask]
H_ion = S_H[:, ion_mask]
He_ion = S_He[:, ion_mask]

# Wavelength runs opposite frequency, so sort it
sort_ion = np.argsort(lam_ion)

lam_ion = lam_ion[sort_ion]
spectra_ion = spectra_ion[:, sort_ion]
H_ion = H_ion[:, sort_ion]
He_ion = He_ion[:, sort_ion]

# ============================================================
# Normalize each composition to unit area
#
# H and He are divided by the SAME total-spectrum area so
# their relative contributions are preserved.
# ============================================================

areas = np.trapezoid(spectra_ion, lam_ion, axis=1)

spectra_norm = spectra_ion / areas[:, None]
H_norm = H_ion / areas[:, None]
He_norm = He_ion / areas[:, None]

# Selected compositions
idx20 = np.argmin(np.abs(fH - 0.2))
idx80 = np.argmin(np.abs(fH - 0.8))

# ============================================================
# Make 2 x 2 figure
# ============================================================

fig, ax = plt.subplots(2, 2, figsize=FIGSIZE)

ax_a = ax[0, 0]
ax_b = ax[1, 0]
ax_c = ax[0, 1]
ax_d = ax[1, 1]

# ============================================================
# (a) Full Thomson spectrum
# ============================================================

ax_a.plot(lam_full, spectra_dummy_norm, color=TOTAL_COLOR)

ax_a.axvline(
    probe_wavelength_nm,
    color=PROBE_COLOR,
    linestyle="--",
    linewidth=LINE_WIDTH,
    alpha=0.9
)

ax_a.set_xlim(475, 576)
ax_a.set_xlabel(r"$\lambda$ (nm)")
ax_a.set_ylabel(r"$S(k, \omega)$ (norm.)")
ax_a.set_title("(a) Electron and ion features")

xticks = np.arange(475, 576, 25)
xticks = xticks[~np.isclose(xticks, 525)]
xticks = np.sort(np.append(xticks, probe_wavelength_nm))

ax_a.set_xticks(xticks)
ax_a.set_xticklabels([
    "526.5" if np.isclose(x, probe_wavelength_nm) else f"{x:.0f}"
    for x in xticks
])

# ============================================================
# (b) 20% Hydrogen
# ============================================================

ax_b.plot(lam_ion, spectra_norm[idx20], color=TOTAL_COLOR, label="Total")
ax_b.plot(lam_ion, H_norm[idx20], color=H_COLOR, linestyle="--", label="H")
ax_b.plot(lam_ion, He_norm[idx20], color=HE_COLOR, linestyle="--", label="He")

ax_b.axvline(probe_wavelength_nm, color=PROBE_COLOR, linestyle=":", linewidth=1.5, alpha=0.8)

ax_b.set_xlabel(r"$\lambda$ (nm)")
ax_b.set_xticks([525, 526, 527, 528])
ax_b.set_ylabel(r"$S(k, \omega)$ (norm.)")
ax_b.set_title(r"(c) $f_H=0.2$")
ax_b.legend(frameon=False)

# ============================================================
# (c) Composition scan
# ============================================================

pcm = ax_c.pcolormesh(
    lam_ion,
    fH,
    spectra_norm,
    shading="auto",
    cmap="turbo"
)

ax_c.axhline(0.2, color="white", linestyle="--", linewidth=1.8)
ax_c.axhline(0.8, color="white", linestyle="--", linewidth=1.8)
ax_c.axvline(probe_wavelength_nm, color="white", linestyle=":", linewidth=1.5, alpha=0.9)

ax_c.set_xlabel(r"$\lambda$ (nm)")
ax_c.set_ylabel(r"$f_H$")
ax_c.set_title("(b) Composition dependence")

cbar = fig.colorbar(pcm, ax=ax_c, pad=0.02)
cbar.set_label(r"$S(k, \omega)$ (norm.)")
cbar.ax.tick_params(labelsize=TICK_SIZE)

# ============================================================
# (d) 80% Hydrogen
# ============================================================

ax_d.plot(lam_ion, spectra_norm[idx80], color=TOTAL_COLOR, label="Total")
ax_d.plot(lam_ion, H_norm[idx80], color=H_COLOR, linestyle="--", label="H")
ax_d.plot(lam_ion, He_norm[idx80], color=HE_COLOR, linestyle="--", label="He")

ax_d.axvline(probe_wavelength_nm, color=PROBE_COLOR, linestyle=":", linewidth=1.5, alpha=0.8)

ax_d.set_xlabel(r"$\lambda$ (nm)")
ax_d.set_xticks([525, 526, 527, 528])
ax_d.set_ylabel(r"$S(k, \omega)$ (norm.)")
ax_d.set_title(r"(d) $f_H=0.8$")
ax_d.legend(frameon=False)

# ============================================================
# Common formatting
# ============================================================

for a in ax.flat:
    a.tick_params(which="both", direction="in")

# Give the three ion-feature panels identical x limits
ION_XLIM = (lam_ion.min(), lam_ion.max())

ax_b.set_xlim(ION_XLIM)
ax_c.set_xlim(ION_XLIM)
ax_d.set_xlim(ION_XLIM)

fig.subplots_adjust(
    left=0.10,
    right=0.96,
    bottom=0.10,
    top=0.95,
    wspace=0.30,
    hspace=0.32
)

# ============================================================
# Save
# ============================================================

fig.savefig("thomson_scattering_overview.pdf")
plt.show()