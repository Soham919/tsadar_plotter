import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
from TS_aux import shot_data
from TS_aux import get_calibration


def plot_92534_EPW(data):

    # ============================================================
    # FIGURE SETTINGS
    # ============================================================

    FIGSIZE = (5, 7)

    FONT_SIZE = 16
    LABEL_SIZE = 18
    TICK_SIZE = 16

    AXIS_WIDTH = 1.2
    TICK_WIDTH = 1.2
    TICK_LENGTH = 5
    MARKER_SIZE = 22

    TE_COLOR = "tab:red"
    NE_COLOR = "tab:blue"

    XMIN, XMAX = 2.25, 4.00

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
    # LOAD EPW SPECTRUM
    # ============================================================

    shot = "92534"
    ts_type = "imaging"

    fp = Path(
        "//profiles/Users$/sban/Documents/Kinshock/Kinshock-19A/Data"
    ) / shot

    f_epw = f"EPW_CCD-s{shot}.h5"

    with h5py.File(fp / f_epw, "r") as f:
        epw = f["Streak_array"][0, :, :]

        # Same orientation corrections as before
        epw = np.rot90(epw, 1)
        epw = np.fliplr(epw)
        epw = np.flipud(epw)

    # Calibration
    x_epw, x_iaw, y_epw, y_iaw = get_calibration(
        shot,
        ts_type
    )

    # Convert spatial calibration to experimental position
    x_epw_plot = (
        shot_data.loc[int(shot), "pointing"]
        - x_epw / 1e3
    )

    # ============================================================
    # LOAD FITTED ELECTRON PARAMETERS
    # ============================================================

    lp = data["learned_parameters.csv"]

    r = (
        shot_data.loc[92534, "pointing"]
        - lp["Radius (\\mum)"] / 1e3
    )

    # Rolling mean / standard deviation
    window = 5

    Te = lp["Te_electron"].rolling(
        window=window,
        center=True
    ).mean()

    ne = lp["ne_electron"].rolling(
        window=window,
        center=True
    ).mean()

    std_Te = lp["Te_electron"].rolling(
        window=window,
        center=True
    ).std()

    std_ne = lp["ne_electron"].rolling(
        window=window,
        center=True
    ).std()

    ## -- Save the erros and stuff to csv -- ##
    output_dir = Path("EPW_errors")
    output_dir.mkdir(parents=True, exist_ok=True)

    results = pd.DataFrame({
        "x_mm" : r,
        "Te_fit" : lp["Te_electron"],
        "Te_mean" : Te,
        "Te_stddev" : std_Te,
        "ne_fit": lp["ne_electron"],
        "ne_mean" : ne,
        "np_stddev" : std_ne,
    })

    results.to_csv(output_dir / f"EPW_92534-10_window{window}.csv", index=False)
    # ============================================================
    # FIGURE
    # ============================================================

    fig, (ax_spec, ax_fit) = plt.subplots(
        2, 1,
        figsize=FIGSIZE,
        sharex=True,
        gridspec_kw={
            "height_ratios": [1.25, 1],
            "hspace": 0.08
        }
    )

    # ============================================================
    # (1) EPW SPECTRUM
    # ============================================================

    ax_spec.pcolormesh(
        x_epw_plot,
        y_epw,
        epw,
        cmap="turbo",
        shading="auto",
        rasterized=True
    )

    # Probe wavelength
    ax_spec.axhline(
        526.5,
        color="white",
        linestyle="--",
        linewidth=1.8
    )

    ax_spec.set_ylabel(r"$\lambda$ (nm)")

    # No x-axis labels here
    ax_spec.tick_params(labelbottom=False)

    # ============================================================
    # (2) ELECTRON TEMPERATURE
    # ============================================================

    ax_fit.scatter(
        r,
        Te,
        color=TE_COLOR,
        s=MARKER_SIZE,
        label=r"$T_e$",
        zorder=3
    )

    ax_fit.fill_between(
        r,
        Te - std_Te,
        Te + std_Te,
        color=TE_COLOR,
        alpha=0.15,
        linewidth=0
    )

    ax_fit.set_xlabel(r"$x$ (mm)")
    ax_fit.set_ylabel(
        r"$T_e$ (keV)",
        color=TE_COLOR
    )

    ax_fit.tick_params(
        axis="y",
        colors=TE_COLOR
    )

    ax_fit.spines["left"].set_color(TE_COLOR)
    ax_fit.spines["left"].set_linewidth(AXIS_WIDTH)
    ax_fit.spines["right"].set_visible(False)

    # ============================================================
    # ELECTRON DENSITY -- SECOND Y AXIS
    # ============================================================

    ax_ne = ax_fit.twinx()

    ax_ne.scatter(
        r,
        ne,
        color=NE_COLOR,
        s=MARKER_SIZE,
        label=r"$n_e$",
        zorder=3
    )

    ax_ne.fill_between(
        r,
        ne - std_ne,
        ne + std_ne,
        color=NE_COLOR,
        alpha=0.15,
        linewidth=0
    )

    ax_ne.set_ylabel(
        r"$n_e$ ($10^{20}$ cm$^{-3}$)",
        color=NE_COLOR
    )

    ax_ne.tick_params(
        axis="y",
        colors=NE_COLOR
    )

    # Each axis owns only its appropriate spine
    ax_ne.spines["left"].set_visible(False)
    ax_ne.spines["right"].set_color(NE_COLOR)
    ax_ne.spines["right"].set_linewidth(AXIS_WIDTH)

    # ============================================================
    # COMMON FORMATTING
    # ============================================================

    ax_spec.set_xlim(XMIN, XMAX)
    ax_fit.set_xlim(XMIN, XMAX)

    ax_fit.set_xticks([
        2.50,
        3.00,
        3.50,
        4.00
    ])

    ax_fit.grid(
        True,
        alpha=0.20,
        linewidth=0.8
    )

    ax_fit.set_axisbelow(True)

    # ============================================================
    # MARGINS
    # ============================================================

    # Need extra room on right for ne axis
    fig.subplots_adjust(
        left=0.19,
        right=0.85,
        bottom=0.10,
        top=0.99,
        hspace=0.08
    )

    # ============================================================
    # SAVE
    # ============================================================

    fig.savefig("EPW_spectrum_and_fits.pdf")

    plt.show()

    return fig, (ax_spec, ax_fit, ax_ne)