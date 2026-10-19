import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")
from Plasma_func.Plasma_func import RH


# -------------------------------------------------
# MATLAB-like plotting defaults
# -------------------------------------------------
plt.rcParams.update({
    "font.family": "Arial",       # close to MATLAB Helvetica
    "font.size": 14,
    "axes.labelsize": 18,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "lines.linewidth": 1.5,
})


# -------------------------------------------------
# Rankine-Hugoniot conditions
# -------------------------------------------------
d = 1
v = 1
T = 1

d2, v2, T2 = RH(10, d, v, T)


# -------------------------------------------------
# Discontinuous profiles
# -------------------------------------------------
x = np.array([-2.0, -1.0, 1.0])

rho = np.array([d2, d, d])
temp = np.array([T2, T, T])


# -------------------------------------------------
# Plot
# -------------------------------------------------
fig, ax_rho = plt.subplots(figsize=(5, 5))

# MATLAB/reference-style colors
blue = (0.0, 0.0, 1.0)
red  = (1.0, 0.0, 0.0)


# -------------------------
# Density — LEFT axis
# -------------------------
ax_rho.step(
    x,
    rho,
    where="post",
    color=blue,
    linewidth=1.5,
    label = r"$n_e$"
)

ax_rho.set_ylabel(r"$n_e$", color=blue)
ax_rho.set_ylim(1, 5)

ax_rho.tick_params(
    axis="y",
    colors=blue,
    direction="out",
)

ax_rho.spines["left"].set_color(blue)


# -------------------------
# Temperature — RIGHT axis
# -------------------------
ax_temp = ax_rho.twinx()

ax_temp.step(
    x,
    temp,
    where="post",
    color=red,
    linewidth=1.5,
    label = r"$T_e$"
)

ax_temp.set_ylabel(r"$T_e$", color=red)
ax_temp.set_ylim(1, 35)

ax_temp.tick_params(
    axis="y",
    colors=red,
    direction="out",
)

ax_temp.spines["right"].set_color(red)


# -------------------------
# X axis
# -------------------------
ax_rho.set_xlabel(r"$x$")

ax_rho.tick_params(
    axis="x",
    direction="out",
)


# Remove numerical tick labels like your reference figure
ax_rho.set_xticklabels([])
ax_rho.set_yticklabels([])
ax_temp.set_yticklabels([])


# Keep MATLAB-style box
ax_rho.spines["top"].set_visible(True)
ax_rho.spines["bottom"].set_visible(True)

ax_rho.spines["top"].set_color("black")
ax_rho.spines["bottom"].set_color("black")

ax_temp.spines["top"].set_visible(False)
ax_temp.spines["bottom"].set_visible(False)

lines1, labels1 = ax_rho.get_legend_handles_labels()
lines2, labels2 = ax_temp.get_legend_handles_labels()

legend = ax_rho.legend(
    lines1 + lines2,
    labels1 + labels2,
    loc="upper right",
    frameon=True,
    fancybox=False,
    framealpha=1.0,
    edgecolor="black",
    facecolor="white",
    fontsize=14,
    handlelength=1.6,
    handletextpad=0.5,
    borderpad=0.35,
)

legend.get_frame().set_linewidth(0.8)
plt.tight_layout()
plt.show()