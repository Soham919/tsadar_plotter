import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants

# ============================================================
# USER SETTINGS
# ============================================================
BASE_DIR = Path(r"C:\Simulation_data\FLASH\2D_Cylindrical\FV4.7.1_carbuncle_debugs\19A_hll1_2.5TW\plots")
BASE_DIR2 = Path(r"C:\Simulation_data\FLASH\1D_spherical\iFP_bubble_T_7.5e10\plots\shock_results.npy")

# Sedov/shock-speed simulation (kept exactly from your first script)
SEDOV_PATH = BASE_DIR / "shock_results.npy"

# Simulation used for the four experiment-comparison panels
COMPARE_PATH = BASE_DIR / "shock_results.npy"

SAVE_NAME = BASE_DIR / "shock_sedov_exp_comparison.pdf"
FIGSIZE = (10, 10)

FONT_SIZE = 16
LABEL_SIZE = 18
TITLE_SIZE = 18
TICK_SIZE = 16
LEGEND_SIZE = 14
AXIS_WIDTH = 1.2
TICK_WIDTH = 1.2
TICK_LENGTH = 5
LINE_WIDTH = 2.0
MARKER_SIZE = 40

FLASH_COLOR = "tab:blue"
FIT_COLOR = "tab:red"
SEDOV_COLOR = "tab:green"
EXP_COLOR = "k"

# These were different in your original speed script, so I kept both explicit.
FIT_START = 1.0       # ns: start of power-law fit
LASER_OFF = 1.0       # ns: location of vertical "Laser off" line

COMPARE_XLIM = (2.0, 9.5)

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
# EXPERIMENTAL POINTS
# ============================================================

experimental = [
# {
#     "label": "Experiment", #2026 - 117830
#     "time": 9.0,
#     "position": 10.0, "position_err": 0.1,
#     "velocity": 0.35,   "velocity_err": 0.02,
#     "tele": 300,      "tele_err": 40,
#     "tion": 300,      "tion_err": 60,
# },

# {
#     "label": "Experiment", #2019 - 92534
#     "time": 4.0,
#     "position": 3.1, "position_err": 0.1,
#     "velocity": 0.29,   "velocity_err": 0.02,
#     "tele": 450,      "tele_err": 50,
#     "tion": 380,      "tion_err": 20,
# },

# {
#     "label": "Experiment", #2019 - 92537
#     "time": 6.0,
#     "position": 4.5, "position_err": 0.1,
#     "velocity": 0.29,   "velocity_err": 0.01,
#     "tele": 375,      "tele_err": 25,
#     "tion": 600,      "tion_err": 100,
# }


]

# ============================================================
# HELPERS
# ============================================================

def load_shock_results(path):
    s = np.load(path, allow_pickle=True)
    time = np.array([r["Time"] for r in s])
    position = np.array([r["shock_coord"] for r in s]) / 1000.0
    flow_velocity = np.array([r["postshock_vely"] for r in s]) / 1e8
    tele = np.array([r["postshock_tele"] for r in s]) * constants.Boltzmann / constants.e
    tion = np.array([r["postshock_tion"] for r in s]) * constants.Boltzmann / constants.e
    return time, position, flow_velocity, tele, tion

def format_ax(ax):
    ax.grid(True, alpha=0.22, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(which="both", direction="in")

# ============================================================
# LOAD DATA
# ============================================================

# --- Data for Sedov/shock-speed panel ---
t_sedov, x_sedov, _, _, _ = load_shock_results(SEDOV_PATH)
u_s = np.gradient(x_sedov, t_sedov)  # mm/ns

p = load_shock_results(BASE_DIR2)

mask = (t_sedov > FIT_START) & (u_s > 0) & np.isfinite(u_s)
t_fit = t_sedov[mask]
u_fit_data = u_s[mask]

coeff = np.polyfit(np.log(t_fit), np.log(u_fit_data), 1)
alpha = coeff[0]
A = np.exp(coeff[1])

u_fit = A * t_fit**alpha
u_sedov = u_fit_data[0] * (t_fit / t_fit[0])**(-0.6)

# --- Data for FLASH-vs-experiment panels ---
t, position, flow_velocity, tele, tion = load_shock_results(COMPARE_PATH)

# ============================================================
# FIGURE LAYOUT
# ============================================================

fig = plt.figure(figsize=FIGSIZE)

# Main grid: top Sedov panel + entire lower 2x2 block
outer = fig.add_gridspec(
    2, 1,
    height_ratios=[1.0, 2.0],
    hspace=0.16          # <-- gap between TOP plot and lower block
)

# Top panel -- make it narrower than full figure
top = outer[0].subgridspec(
    1, 6,
    wspace=0
)

ax_speed = fig.add_subplot(top[0, 1:5])

# Lower 2x2 block
lower = outer[1].subgridspec(
    2, 2,
    hspace=0.12,         # <-- gap BETWEEN lower rows
    wspace=0.20          # <-- gap between left/right columns
)

ax_pos = fig.add_subplot(lower[0, 0])
ax_vel = fig.add_subplot(lower[0, 1])
ax_te  = fig.add_subplot(lower[1, 0])
ax_ti  = fig.add_subplot(lower[1, 1])

# ============================================================
# (a) SHOCK SPEED / SEDOV
# ============================================================

ax_speed.scatter(t_sedov, u_s, s=MARKER_SIZE, color=FLASH_COLOR, label="FLASH", zorder=3)
ax_speed.plot(t_fit, u_fit, "--", color=FIT_COLOR, label=fr"Best fit: $u_s \propto t^{{{alpha:.2f}}}$")
ax_speed.plot(t_fit, u_sedov, ":", color=SEDOV_COLOR, label=r"3D Sedov: $u_s \propto t^{-3/5}$")
ax_speed.axvline(LASER_OFF, color="0.4", linestyle="--", linewidth=1.5, label="Laser off")

ax_speed.set_ylabel(r"$u_s$ (mm/ns)")
ax_speed.set_xlabel("Time (ns)")
ax_speed.legend(
    frameon=False,
    ncol=1,
    loc="upper right",
    fontsize=12
)
format_ax(ax_speed)

# ============================================================
# (b-e) FLASH vs EXPERIMENT
# ============================================================

ax_pos.plot(t, position-2, color=FLASH_COLOR, label="FLASH")
ax_vel.plot(t, flow_velocity, color=FLASH_COLOR)
ax_te.plot(t, tele/1000, color=FLASH_COLOR)
ax_ti.plot(t, tion/1000, color=FLASH_COLOR)

for exp in experimental:
    common = dict(fmt="o", color=EXP_COLOR, ecolor=EXP_COLOR, capsize=4, zorder=10)
    ax_pos.errorbar(exp["time"], exp["position"], yerr=exp["position_err"], **common)
    ax_vel.errorbar(exp["time"], exp["velocity"], yerr=exp["velocity_err"]/10, **common)
    ax_te.errorbar(exp["time"], exp["tele"]/1000, yerr=exp["tele_err"]/1000, **common)
    ax_ti.errorbar(exp["time"], exp["tion"]/1000, yerr=exp["tion_err"]/1000, **common)

ax_pos.set_ylabel(r"$x_s$ (mm)", labelpad=0)
ax_vel.set_ylabel(r"$v_f$ (mm/ns)", labelpad=-5)
ax_te.set_ylabel(r"$T_e$ (keV)", labelpad=0)
ax_ti.set_ylabel(r"$T_i$ (keV)", labelpad=0)

ax_te.set_xlabel("Time (ns)")
ax_ti.set_xlabel("Time (ns)")

for ax in (ax_pos, ax_vel, ax_te, ax_ti):
    ax.set_xlim(*COMPARE_XLIM)
    format_ax(ax)

# Hide duplicated x tick labels on middle row
ax_pos.tick_params(labelbottom=False)
ax_vel.tick_params(labelbottom=False)

ax_vel.set_ylim(0, 1)
ax_te.set_ylim(0, 1)
ax_ti.set_ylim(0, 1)

# A single legend is enough for the lower four panels
ax_pos.legend(frameon=False, loc="best")

# Optional panel labels -- easier to reference in the LaTeX caption
for ax, label in zip(
    [ax_pos, ax_vel, ax_te, ax_ti],
    ["(b)", "(c)", "(d)", "(e)"]
):
    ax.text(0.1, 0.94, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=FONT_SIZE, fontweight="bold")

ax_speed.text(0.02, 0.94, "(a)", transform=ax_speed.transAxes, ha="left", va="top",
               fontsize=FONT_SIZE, fontweight="bold")
# ============================================================
# SAVE
# ============================================================

fig.subplots_adjust(left=0.11, right=0.97, bottom=0.08, top=0.98)
fig.savefig(SAVE_NAME)
plt.show()
