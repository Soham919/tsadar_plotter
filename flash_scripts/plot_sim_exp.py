import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants
eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h

# ============================================================
# SIMULATIONS
# ============================================================

simulations = [
    #{"path": Path(r"C:\Simulation_data\FLASH\2D\FV4.7.1_carbuncle_debugs\19A_GJS_2.5TW\plots\shock_results.npy"), "label": "Ablation driven", "color":"blue"},
    {"path": Path(r"C:\Simulation_data\FLASH\2D\FV4.7.1_carbuncle_debugs\26A_GJS_5TW\plots\shock_results.npy"), "label": "FLASH", "color":"tab:blue"},
]


# ============================================================
# EXPERIMENTAL POINTS
#
# Each point can have its own time and values/errors.
# ============================================================

experimental = [
    # {
    #     "label": "Shot 1",
    #     "time": 4.0,
    #     "position": 5.1, "position_err": 0.1,
    #     "velocity": 29,  "velocity_err": 2,
    #     "tele": 500,     "tele_err": 100,
    #     "tion": 700,     "tion_err": 50,
    #     "color": "blue"
    # },
    # {
    #     "label": "Shot 2",
    #     "time": 6.0,
    #     "position": 6.6, "position_err": 0.1,
    #     "velocity": 30,  "velocity_err": 2,
    #     "tele": 400,     "tele_err": 50,
    #     "tion": 600,     "tion_err": 50,
    #     "color": "blue"
    # },
    {
        "label": "Shot 3",
        "time": 9.0,
        "position": 10.0, "position_err": 0.1,
        "velocity": 35,  "velocity_err": 2,
        "tele": 300,     "tele_err": 40,
        "tion": 300,     "tion_err": 60,
        "color": "k"
    },
]


# ============================================================
# CREATE FIGURE
# ============================================================

plt.rcParams.update({
    #"font.family": "sans-serif",
    #"font.sans-serif": ["Arial", "DejaVu Sans"],

    "font.size": 16,
    "axes.labelsize": 18,
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

fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)

ax_pos = axes[0, 0]
ax_vel = axes[0, 1]
ax_te  = axes[1, 0]
ax_ti  = axes[1, 1]


# ============================================================
# PLOT SIMULATIONS
# ============================================================

for sim in simulations:

    shock_res = np.load(Path(sim["path"]), allow_pickle=True)

    time = np.array([r["Time"] for r in shock_res])
    position = np.array([r["shock_coord"] for r in shock_res])/1000
    velocity = np.array([r["postshock_vely"] for r in shock_res]) / 1e6   # cm/s -> km/s
    tele = np.array([r["postshock_tele"] for r in shock_res]) * kb/e
    tion = np.array([r["postshock_tion"] for r in shock_res]) * kb/e

    ax_pos.plot(time, position, lw=2, label=sim["label"], color = sim["color"])
    ax_vel.plot(time, velocity, lw=2, label=sim["label"],color = sim["color"])
    ax_te.plot(time, tele, lw=2, label=sim["label"], color = sim["color"])
    ax_ti.plot(time, tion, lw=2, label=sim["label"], color = sim["color"])


# ============================================================
# PLOT EXPERIMENTAL POINTS
# ============================================================

for exp in experimental:

    ax_pos.errorbar(exp["time"], exp["position"], yerr=exp["position_err"],
                    fmt="o", color=exp["color"], ecolor="black", capsize=4, zorder=10)

    ax_vel.errorbar(exp["time"], exp["velocity"], yerr=exp["velocity_err"],
                    fmt="o", color=exp["color"], ecolor="black", capsize=4, zorder=10)

    ax_te.errorbar(exp["time"], exp["tele"], yerr=exp["tele_err"],
                   fmt="o", color=exp["color"], ecolor="black", capsize=4, zorder=10)

    ax_ti.errorbar(exp["time"], exp["tion"], yerr=exp["tion_err"],
                   fmt="o", color=exp["color"], ecolor="black", capsize=4, zorder=10)


# ============================================================
# LABELS
# ============================================================

ax_pos.set_ylabel(r"$x_s$ (mm)")
ax_vel.set_ylabel(r"$v_{f}$ ($\times 10^{6}$ cm/s)")
ax_te.set_ylabel(r"$T_e$ (eV)")
ax_ti.set_ylabel(r"$T_i$ (eV)")
ax_pos.set_title(r"$Shock$ $position$")
ax_vel.set_title(r"$Flow$ $velocity$")
ax_te.set_title(r"$T_e$ ")
ax_ti.set_title(r"$T_i$ ")
#ax_ti.set_ylim([0,1000])
ax_pos.set_xlim([2,9.5])
ax_vel.set_xlim([2,9.5])
ax_vel.set_ylim([0,100])
ax_te.set_xlim([2,9.5])
ax_te.set_ylim([0,3000])
ax_ti.set_xlim([2,9.5])
ax_ti.set_ylim([0,1000])

ax_te.set_xlabel("Time (ns)")
ax_ti.set_xlabel("Time (ns)")

for ax in axes.flat:
    ax.grid(alpha=0.25)
    ax.tick_params(direction="in")

ax_pos.legend(frameon=False)

fig.tight_layout()
plt.show()