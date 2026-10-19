import numpy as np
import matplotlib.pyplot as plt
import sys

# ----- Mac ----- #
# sys.path.append("/Users/soham/Documents/Kinshock/Scripts")

# ----- Windows ----- #
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")

from plasma_chi import *
from lineout_tool import InteractiveLineout
from matplotlib.colors import LogNorm

# ============================================================
# Probe geometry
# ============================================================
probe_wavelength_nm = 526.5
scattering_angle_deg = 60.0

probe_wavelength_m = probe_wavelength_nm * 1e-9
scattering_angle_rad = np.deg2rad(scattering_angle_deg)

k = (4.0 * np.pi / probe_wavelength_m) * np.sin(scattering_angle_rad / 2.0)

# ============================================================
# Frequency / phase-velocity axis
# ============================================================
phase_velocity = np.linspace(-7.0e7, 7.0e7, 10000)

# ============================================================
# Plasma conditions
# ============================================================
ne_cm3 = 1e20
ne = ne_cm3 * 1e6

Te_eV = 500
TH_eV = 300
THe_eV = 300

ud_e = 0.0
ud_H = 0.0
ud_He = 0.0

ZH = 1.0
ZHe = 2.0

fH = np.linspace(0.01, 0.99, 500)
fHe = 1 - fH

results = []

# ============================================================
# Calculate Thomson spectra
# ============================================================
for f in fH:

    Zbar = ZH * f + ZHe * (1 - f)

    electron = PlasmaSpecies(name='electron', Z=-1, m=me, T_eV=Te_eV, n_cm3=ne_cm3, ud=ud_e*1e4)
    hydrogen = PlasmaSpecies(name='hydrogen', Z=ZH, m=mp, T_eV=TH_eV, n_cm3=f*(ne_cm3/Zbar), ud=ud_H*1e4)
    helium = PlasmaSpecies(name='helium', Z=ZHe, m=4.0*mp, T_eV=THe_eV, n_cm3=(1-f)*(ne_cm3/Zbar), ud=ud_He*1e4)

    ions = [hydrogen, helium]

    Result = thomson_structure_factor(k=k, phase_velocity=phase_velocity, electron=electron, ions=ions)
    results.append(Result)

spectra = np.array([result["S_total"] for result in results])
eps = np.array([result["epsilon"] for result in results])
S_i = np.array([result["S_ions"] for result in results])

# ============================================================
# Plotting settings
# ============================================================
ion_mask = np.abs(phase_velocity) <= 1.0e6

# Convert m/s to units of 10^6 cm/s
v_plot = phase_velocity[ion_mask] / 1e4

H_COLOR = "orangered"
HE_COLOR = "mediumseagreen"
TOTAL_COLOR = "navy"

plt.rcParams.update({
    #"font.family": "sans-serif",
    #"font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 18,
    "axes.labelsize": 16,
    "axes.titlesize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "axes.linewidth": 1.2,
})

fig, ax = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

# ============================================================
# Middle panel: Spectra versus hydrogen fraction
# ============================================================
pcm = ax[1].pcolormesh(v_plot, fH, spectra[:, ion_mask]*1e13, shading='auto', cmap='turbo')

ax[1].set_xlabel(r'Phase velocity ($10^6$ cm/s)')
ax[1].set_ylabel(r'$f_H$')
ax[1].set_title('Thomson Spectra')
ax[1].axhline(0.1, linestyle = '--', linewidth = 2, color = 'w')
ax[1].axhline(0.9, linestyle = '--', linewidth = 2, color = 'w')

fig.colorbar(pcm, ax=ax[1], label='Spectra')

# ============================================================
# Left panel: 20% hydrogen
# ============================================================
idx = np.argmin(abs(fH - 0.2))

h = S_i[idx]['hydrogen']
he = S_i[idx]['helium']

ax[0].plot(v_plot, spectra[idx, ion_mask]*1e13, color=TOTAL_COLOR, linewidth=2, label="Total")
ax[0].plot(v_plot, h[ion_mask]*1e13, color=H_COLOR, linewidth=2, linestyle = '--',label="H")
ax[0].plot(v_plot, he[ion_mask]*1e13, color=HE_COLOR, linewidth=2, linestyle = '--',label="He")

ax[0].set_xlabel(r'Phase velocity ($10^6$ cm/s)')
ax[0].set_ylabel('Spectra (a.u.)')
ax[0].set_title('20% Hydrogen')
ax[0].set_xlim([-100, 100])
ax[0].legend(frameon=False, fontsize=14)

# ============================================================
# Right panel: 80% hydrogen
# ============================================================
idx = np.argmin(abs(fH - 0.8))

h = S_i[idx]['hydrogen']
he = S_i[idx]['helium']

ax[2].plot(v_plot, spectra[idx, ion_mask]*1e13, color=TOTAL_COLOR, linewidth=2, label="Total")
ax[2].plot(v_plot, h[ion_mask]*1e13, color=H_COLOR, linewidth=2, linestyle = '--',label="H")
ax[2].plot(v_plot, he[ion_mask]*1e13, color=HE_COLOR, linewidth=2, linestyle = '--',label="He")

ax[2].set_xlabel(r'Phase velocity ($10^6$ cm/s)')
ax[2].set_ylabel('Spectra (a.u.)')
ax[2].set_title('80% Hydrogen')
ax[2].set_xlim([-100, 100])
ax[2].legend(frameon=False, fontsize=14)

plt.show()