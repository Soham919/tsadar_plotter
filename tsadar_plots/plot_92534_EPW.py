### PLOT LEARNED PARAMETERS FOR ELECTRONS ###
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
from pathlib import Path
import re
import pandas as pd
from TS_aux import load_file
from TS_aux import shot_data

def plot_92534_EPW(data):

    # ============================================================
    # FIGURE SETTINGS
    # ============================================================

    FIGSIZE = (5, 5)

    FONT_SIZE = 16
    LABEL_SIZE = 18
    TICK_SIZE = 16

    AXIS_WIDTH = 1.2
    TICK_WIDTH = 1.2
    TICK_LENGTH = 5
    MARKER_SIZE = 22

    TE_COLOR = "tab:red"
    NE_COLOR = "tab:blue"

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
    })

    # ============================================================
    # Data
    # ============================================================

    lp = data["learned_parameters.csv"]

    r = shot_data.loc[92534, "pointing"] - lp["Radius (\\mum)"] / 1e3

    # Rolling mean / standard deviation
    window = 5

    Te = lp["Te_electron"].rolling(window=window, center=True).mean()
    ne = lp["ne_electron"].rolling(window=window, center=True).mean()

    std_Te = lp["Te_electron"].rolling(window=window, center=True).std()
    std_ne = lp["ne_electron"].rolling(window=window, center=True).std()

    # ============================================================
    # Figure
    # ============================================================

    fig, ax1 = plt.subplots(figsize=FIGSIZE)

    # ------------------------------------------------------------
    # Electron temperature
    # ------------------------------------------------------------

    ax1.scatter(
        r, Te,
        color=TE_COLOR,
        s=MARKER_SIZE,
        label=r"$T_e$",
        zorder=3
    )

    ax1.fill_between(
        r,
        Te - std_Te,
        Te + std_Te,
        color=TE_COLOR,
        alpha=0.15,
        linewidth=0
    )

    ax1.set_xlabel(r"$x$ (mm)")
    ax1.set_ylabel(r"$T_e$ (keV)", color=TE_COLOR)

    ax1.tick_params(axis="y", colors=TE_COLOR)
    ax1.spines["left"].set_color(TE_COLOR)

    # ------------------------------------------------------------
    # Electron density
    # ------------------------------------------------------------

    ax2 = ax1.twinx()

    ax2.spines["left"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax1.spines["left"].set_color(TE_COLOR)
    ax1.spines["left"].set_linewidth(AXIS_WIDTH)

    ax2.spines["right"].set_color(NE_COLOR)
    ax2.spines["right"].set_linewidth(AXIS_WIDTH)
    ax2.scatter(
        r, ne,
        color=NE_COLOR,
        s=MARKER_SIZE,
        label=r"$n_e$",
        zorder=3
    )

    ax2.fill_between(
        r,
        ne - std_ne,
        ne + std_ne,
        color=NE_COLOR,
        alpha=0.15,
        linewidth=0
    )

    ax2.set_ylabel(
        r"$n_e$ ($10^{20}$ cm$^{-3}$)",
        color=NE_COLOR
    )

    ax2.tick_params(axis="y", colors=NE_COLOR)
    ax2.spines["right"].set_color(NE_COLOR)

    # ============================================================
    # Formatting
    # ============================================================

    # Match EPW spatial range
    ax1.set_xlim(2.25, 4.0)

    ax1.set_xticks([2.50, 3.00, 3.50, 4.00])

    # Light grid behind data
    ax1.grid(True, alpha=0.20, linewidth=0.8)
    ax1.set_axisbelow(True)

    # Keep both vertical spines visible and colored
    ax1.spines["left"].set_linewidth(AXIS_WIDTH)
    ax2.spines["right"].set_linewidth(AXIS_WIDTH)

    # No title -- use LaTeX caption

    # ============================================================
    # Margins
    # ============================================================

    # Slightly more right margin than EPW because of second y-axis
    fig.subplots_adjust(
    left=0.18,
    right=0.90,
    bottom=0.16,
    top=0.97
    )
    # ============================================================
    # Save
    # ============================================================

    fig.savefig("electron_fit.pdf")

    plt.show()

    return fig, (ax1, ax2)
