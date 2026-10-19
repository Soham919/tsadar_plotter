
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants

# ============================================================
# SETTINGS
# ============================================================

BASE_DIR = Path(r"C:\Simulation_data\FLASH\2D\FV4.7.1_carbuncle_debugs\26A_GJS_5TW\plots")
SPHERICAL_PATH = Path(r"C:\Simulation_data\FLASH\1D_spherical\iFP_bubble_T_7.5e10\plots\shock_results.npy")
COMPARE_PATH = BASE_DIR / "shock_results.npy"
SAVE_NAME = BASE_DIR / "shock_sedov_exp_comparison.pdf"

FIT_START, LASER_OFF = 1.0, 1.0
COMPARE_XLIM = (2.0, 9.5)
POSITION_OFFSET_2D, POSITION_OFFSET_1D = -2.0, 0.0

BLUE, ORANGE, RED, GREEN = "tab:blue", "tab:orange", "tab:red", "tab:green"

plt.rcParams.update({
    "font.family": "serif", "font.size": 16, "mathtext.fontset": "cm",
    "axes.labelsize": 18, "axes.titlesize": 18, "axes.linewidth": 1.2,
    "xtick.labelsize": 16, "ytick.labelsize": 16,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.major.width": 1.2, "ytick.major.width": 1.2,
    "xtick.major.size": 5, "ytick.major.size": 5,
    "legend.fontsize": 12, "lines.linewidth": 2
})

# ============================================================
# EXPERIMENTAL POINTS
# ============================================================

experimental = [
    # {"time": 9.0, "position": 10.0, "position_err": 0.1,
    #  "velocity": 0.35, "velocity_err": 0.02,
    #  "tele": 300, "tele_err": 40, "tion": 300, "tion_err": 60},

    # {"time": 4.0, "position": 3.1, "position_err": 0.1,
    #  "velocity": 0.29, "velocity_err": 0.02,
    #  "tele": 450, "tele_err": 50, "tion": 380, "tion_err": 20},

    # {"time": 6.0, "position": 4.5, "position_err": 0.1,
    #  "velocity": 0.29, "velocity_err": 0.01,
    #  "tele": 375, "tele_err": 25, "tion": 600, "tion_err": 100}
]

# ============================================================
# FUNCTIONS
# ============================================================

def load_shock_results(path):
    s = np.load(path, allow_pickle=True)
    if s.ndim == 0: s = [s.item()]

    t = np.array([r["Time"] for r in s], dtype=float)
    x = np.array([r["shock_coord"] for r in s], dtype=float) / 1000
    v = np.array([r["postshock_vely"] for r in s], dtype=float) / 1e8
    te = np.array([r["postshock_tele"] for r in s], dtype=float) * constants.Boltzmann / constants.e / 1000
    ti = np.array([r["postshock_tion"] for r in s], dtype=float) * constants.Boltzmann / constants.e / 1000

    order = np.argsort(t)
    return t[order], x[order], v[order], te[order], ti[order]


def get_shock_speed(t, x):
    mask = np.isfinite(t) & np.isfinite(x)
    t, x = t[mask], x[mask]
    if len(t) < 2: return t, np.full_like(t, np.nan)
    return t, np.gradient(x, t)


def format_ax(ax):
    ax.grid(True, alpha=0.22, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(which="both", direction="in")


# ============================================================
# LOAD DATA
# ============================================================

t, x, v, te, ti = load_shock_results(COMPARE_PATH)
ts, xs, vs, tes, tis = load_shock_results(SPHERICAL_PATH)

tu, u = get_shock_speed(t, x)
tus, us = get_shock_speed(ts, xs)

mask = (tu > FIT_START) & (u > 0) & np.isfinite(u)
tf, uf = tu[mask], u[mask]

if len(tf) >= 2:
    alpha, logA = np.polyfit(np.log(tf), np.log(uf), 1)
    ufit = np.exp(logA) * tf**alpha
    usedov = uf[0] * (tf / tf[0])**(-0.6)

# ============================================================
# FIGURE LAYOUT
# ============================================================

fig = plt.figure(figsize=(10, 10))
outer = fig.add_gridspec(2, 1, height_ratios=[1, 2], hspace=0.16)
top = outer[0].subgridspec(1, 6)
lower = outer[1].subgridspec(2, 2, hspace=0.12, wspace=0.20)

ax_speed = fig.add_subplot(top[0, 1:5])
ax_pos = fig.add_subplot(lower[0, 0])
ax_vel = fig.add_subplot(lower[0, 1])
ax_te = fig.add_subplot(lower[1, 0])
ax_ti = fig.add_subplot(lower[1, 1])

# ============================================================
# (a) SHOCK SPEED
# ============================================================

ax_speed.scatter(tu, u, s=40, color=BLUE, label="2D FLASH", zorder=3)
ax_speed.plot(tus, us, "o--", color=ORANGE, markersize=4, label="1D spherical FLASH")

if len(tf) >= 2:
    ax_speed.plot(tf, ufit, "--", color=RED, label=fr"Best fit: $u_s \propto t^{{{alpha:.2f}}}$")
    ax_speed.plot(tf, usedov, ":", color=GREEN, label=r"3D Sedov: $u_s \propto t^{-3/5}$")

ax_speed.axvline(LASER_OFF, color="0.4", linestyle="--", linewidth=1.5, label="Laser off")
ax_speed.set(xlabel="Time (ns)", ylabel=r"$u_s$ (mm/ns)")
ax_speed.legend(frameon=False, loc="upper right")
format_ax(ax_speed)

# ============================================================
# (b-e) SHOCK POSITION, FLOW VELOCITY, Te, Ti
# ============================================================

axes = [ax_pos, ax_vel, ax_te, ax_ti]
data_2d = [x + POSITION_OFFSET_2D, v, te, ti]
data_1d = [xs + POSITION_OFFSET_1D, vs, tes, tis]
labels = [r"$x_s$ (mm)", r"$v_f$ (mm/ns)", r"$T_e$ (keV)", r"$T_i$ (keV)"]

for ax, y2d, y1d, ylabel in zip(axes, data_2d, data_1d, labels):
    ax.plot(t, y2d, color=BLUE, label="2D FLASH")
    ax.plot(ts, y1d, "--", color=ORANGE, label="1D spherical FLASH")
    ax.set_ylabel(ylabel)
    ax.set_xlim(*COMPARE_XLIM)
    format_ax(ax)

# ============================================================
# EXPERIMENTAL OVERLAY
# ============================================================

for exp in experimental:
    common = dict(fmt="o", color="k", ecolor="k", capsize=4, zorder=10)
    ax_pos.errorbar(exp["time"], exp["position"], yerr=exp["position_err"], **common)
    ax_vel.errorbar(exp["time"], exp["velocity"], yerr=exp["velocity_err"]/10, **common)
    ax_te.errorbar(exp["time"], exp["tele"]/1000, yerr=exp["tele_err"]/1000, **common)
    ax_ti.errorbar(exp["time"], exp["tion"]/1000, yerr=exp["tion_err"]/1000, **common)

# ============================================================
# FINAL FORMATTING
# ============================================================

ax_pos.tick_params(labelbottom=False)
ax_vel.tick_params(labelbottom=False)
ax_te.set_xlabel("Time (ns)")
ax_ti.set_xlabel("Time (ns)")

for ax in (ax_vel, ax_te, ax_ti): ax.set_ylim(0, 1)

ax_pos.legend(frameon=False, loc="best")

for ax, label in zip([ax_speed] + axes, ["(a)", "(b)", "(c)", "(d)", "(e)"]):
    ax.text(0.02 if ax is ax_speed else 0.1, 0.94, label, transform=ax.transAxes,
            ha="left", va="top", fontsize=16, fontweight="bold")

fig.subplots_adjust(left=0.11, right=0.97, bottom=0.08, top=0.98)
fig.savefig(SAVE_NAME, bbox_inches="tight")
plt.show()
