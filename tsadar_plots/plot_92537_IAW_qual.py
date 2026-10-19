import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from TS_aux import shot_data, get_calibration


def plot_92537_IAW(data):

    # ============================================================
    # FIGURE SETTINGS
    # ============================================================

    FIGSIZE = (5, 10)

    FONT_SIZE = 16
    LABEL_SIZE = 17
    TICK_SIZE = 15
    LEGEND_SIZE = 14

    AXIS_WIDTH = 1.2
    TICK_WIDTH = 1.2
    TICK_LENGTH = 5
    MARKER_SIZE = 22

    H_COLOR = "orangered"
    HE_COLOR = "mediumseagreen"

    XMIN, XMAX = 3.55, 4.75

    plt.rcParams.update({
        "font.family": "serif",
        "font.size": FONT_SIZE,
        "mathtext.fontset": "cm",
        "axes.labelsize": LABEL_SIZE,
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
    })

    # ============================================================
    # LOAD IAW DATA
    # ============================================================

    shot = "92537"
    ts_type = "imaging"

    fp = Path("//profiles/Users$/sban/Documents/Kinshock/Kinshock-19A/Data") / shot
    f_iaw = f"IAW_CCD-s{shot}.h5"

    with h5py.File(fp / f_iaw, "r") as f:
        iaw = f["Streak_array"][0, :, :]
        iaw = np.rot90(iaw, 1)
        iaw = np.fliplr(iaw)

    x_epw, x_iaw, y_epw, y_iaw = get_calibration(shot, ts_type)

    # Convert calibrated spatial axis to experimental x position
    x_iaw_plot = shot_data.loc[int(shot), "pointing"] - x_iaw / 1e3

    # ============================================================
    # LOAD FIT PARAMETERS
    # ============================================================

    im1 = data["92537-5"]["92537_5_nersc.csv"]

    r = shot_data.loc[92537, "pointing"] - im1["Radius (\\mum)"] / 1e3

    window = 2

    # Fractions
    f_H  = im1["fract_ion-1"].rolling(window=window, center=True).mean()
    f_He = im1["fract_ion-2"].rolling(window=window, center=True).mean()

    std_f_H  = im1["fract_ion-1"].rolling(window=window, center=True).std()
    std_f_He = im1["fract_ion-2"].rolling(window=window, center=True).std()

    # Temperatures
    T_H  = im1["Ti_ion-1"].rolling(window=window, center=True).mean()
    T_He = im1["Ti_ion-2"].rolling(window=window, center=True).mean()

    std_T_H  = im1["Ti_ion-1"].rolling(window=window, center=True).std()
    std_T_He = im1["Ti_ion-2"].rolling(window=window, center=True).std()

    # Velocities
    V_H  = im1["Va_ion-1"].rolling(window=window, center=True).mean()
    V_He = im1["Va_ion-2"].rolling(window=window, center=True).mean()

    std_V_H  = im1["Va_ion-1"].rolling(window=window, center=True).std()
    std_V_He = im1["Va_ion-2"].rolling(window=window, center=True).std()

    # ============================================================
    # FIGURE
    # ============================================================

    fig, ax = plt.subplots(
        4, 1,
        figsize=FIGSIZE,
        sharex=True,
        gridspec_kw={
            "height_ratios": [1.35, 1, 1, 1],
            "hspace": 0.08
        }
    )

    # ============================================================
    # (1) IAW SPECTRUM
    # ============================================================

    ax[0].pcolormesh(
        x_iaw_plot,
        y_iaw,
        iaw,
        cmap="turbo",
        shading="auto",
        rasterized=True
    )

    ax[0].axhline(
        526.5,
        color="white",
        linestyle="--",
        linewidth=1.8
    )

    ax[0].set_ylabel(r"$\lambda$ (nm)")
    ax[0].set_ylim(524, 528)

    # No x labels/tick labels on this panel
    ax[0].tick_params(labelbottom=False)

    # ============================================================
    # (2) VELOCITIES
    # ============================================================

    ax[1].scatter(r, V_H, color=H_COLOR, s=MARKER_SIZE, label="H", zorder=3)
    ax[1].scatter(r, V_He, color=HE_COLOR, s=MARKER_SIZE, label="He", zorder=3)

    ax[1].fill_between(
        r, V_H-std_V_H, V_H+std_V_H,
        color=H_COLOR, alpha=0.15, linewidth=0
    )

    ax[1].fill_between(
        r, V_He-std_V_He, V_He+std_V_He,
        color=HE_COLOR, alpha=0.15, linewidth=0
    )

    ax[1].set_ylabel(r"$V_i$ ($10^6$ cm/s)")
    ax[1].legend(frameon=False, loc="best")
    ax[1].tick_params(labelbottom=False)

    # ============================================================
    # (3) FRACTIONS
    # ============================================================

    ax[2].scatter(r, f_H, color=H_COLOR, s=MARKER_SIZE, zorder=3)
    ax[2].scatter(r, f_He, color=HE_COLOR, s=MARKER_SIZE, zorder=3)

    ax[2].fill_between(
        r, f_H-std_f_H, f_H+std_f_H,
        color=H_COLOR, alpha=0.15, linewidth=0
    )

    ax[2].fill_between(
        r, f_He-std_f_He, f_He+std_f_He,
        color=HE_COLOR, alpha=0.15, linewidth=0
    )

    ax[2].set_ylabel("Ion fraction")
    ax[2].set_ylim(0, 1)
    ax[2].tick_params(labelbottom=False)

    # ============================================================
    # (4) TEMPERATURES
    # ============================================================

    ax[3].scatter(r, T_H, color=H_COLOR, s=MARKER_SIZE, zorder=3)
    ax[3].scatter(r, T_He, color=HE_COLOR, s=MARKER_SIZE, zorder=3)

    ax[3].fill_between(
        r, T_H-std_T_H, T_H+std_T_H,
        color=H_COLOR, alpha=0.15, linewidth=0
    )

    ax[3].fill_between(
        r, T_He-std_T_He, T_He+std_T_He,
        color=HE_COLOR, alpha=0.15, linewidth=0
    )

    ax[3].set_ylabel(r"$T_i$ (keV)")
    ax[3].set_xlabel(r"$x$ (mm)")
    ax[3].set_ylim(0.2, 0.85)

    # ============================================================
    # COMMON FORMATTING
    # ============================================================

    for a in ax:
        a.set_xlim(XMIN, XMAX)
        a.tick_params(which="both", direction="in")

    # Grid only on inferred-parameter panels
    for a in ax[1:]:
        a.grid(True, alpha=0.20, linewidth=0.8)
        a.set_axisbelow(True)

    fig.subplots_adjust(
        left=0.22,
        right=0.97,
        bottom=0.08,
        top=0.99,
        hspace=0.08
    )

    # ============================================================
    # SAVE
    # ============================================================

    fig.savefig("IAW_spectrum_and_fits.pdf")

    plt.show()

    return fig, ax