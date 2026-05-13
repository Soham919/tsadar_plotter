import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path


# ============================================================
# User options
# ============================================================

cache_file1 = Path("ni_2D_1umAl.npy")
cache_file2 = Path("ne_2D_1umAl.npy")

field_name = r"FLASH 2D ion density"
field_name2 = r"FLASH 2D"

# Domain in FLASH cgs units [cm]
# Change these to match your FLASH parfile.
xmin = -500.0e-04
xmax = 500.0e-04   # cm

ymin = 0.0
ymax = 300.0e-04   # cm

# Convert cm to microns for plotting
useMicrons = True

# Downsampling factors
# Larger number = fewer points = faster/lighter plotting
sx = 1      # downsample x direction
sy = 20     # downsample y direction

# Log plot?
useLog = True

# Optional manual color limits
vmin_manual = None
vmax_manual = None


# ============================================================
# Helper functions
# ============================================================

def make_edges(xmin, xmax, n):
    """
    Make cell-edge coordinates from domain limits and number of cells.
    """
    return np.linspace(xmin, xmax, n + 1)

def downsample_for_pcolormesh(arr2d, x_edges, y_edges, sx=1, sy=1):
    """
    Downsample 2D field and matching edge arrays for pcolormesh.

    arr2d shape: (nx, ny)
    pcolormesh with arr_plot.T expects:
        x_edges length = arr_plot.shape[0] + 1
        y_edges length = arr_plot.shape[1] + 1
    """

    nx, ny = arr2d.shape

    # Cell indices we keep
    ix = np.arange(0, nx, sx)
    iy = np.arange(0, ny, sy)

    arr_plot = arr2d[np.ix_(ix, iy)]

    # Corresponding edges: use kept cell starts plus the final outer edge
    x_plot_edges = x_edges[np.r_[ix, nx]]
    y_plot_edges = y_edges[np.r_[iy, ny]]

    return arr_plot, x_plot_edges, y_plot_edges

def choose_norm(arr_plot, useLog=True, vmin_manual=None, vmax_manual=None):
    """
    Choose LogNorm or Normalize safely.
    """
    finite = np.isfinite(arr_plot)

    if not np.any(finite):
        raise ValueError("Array contains no finite values.")

    if useLog:
        positive = arr_plot[np.isfinite(arr_plot) & (arr_plot > 0)]

        if positive.size == 0:
            print("Warning: no positive values found. Using linear scale.")
            vmin = np.nanmin(arr_plot[finite]) if vmin_manual is None else vmin_manual
            vmax = np.nanmax(arr_plot[finite]) if vmax_manual is None else vmax_manual
            return Normalize(vmin=vmin, vmax=vmax)

        vmin = np.nanmin(positive) if vmin_manual is None else vmin_manual
        vmax = np.nanmax(arr_plot[finite]) if vmax_manual is None else vmax_manual

        return LogNorm(vmin=vmin, vmax=vmax)

    else:
        vmin = np.nanmin(arr_plot[finite]) if vmin_manual is None else vmin_manual
        vmax = np.nanmax(arr_plot[finite]) if vmax_manual is None else vmax_manual

        return Normalize(vmin=vmin, vmax=vmax)


# ============================================================
# Main
# ============================================================

print(f"Loading {cache_file1}...")
arr2d = np.load(cache_file1, mmap_mode="r")
arr2d_2 = np.load(cache_file2, mmap_mode="r")

print("array shape:", arr2d.shape)
print("array dtype:", arr2d.dtype)

nx, ny = arr2d.shape

x_edges = make_edges(xmin, xmax, nx)
y_edges = make_edges(ymin, ymax, ny)

if useMicrons:
    x_edges_plot = x_edges * 1e4
    y_edges_plot = y_edges * 1e4
    xlabel = r"$x$ [$\mu$m]"
    ylabel = r"$y$ [$\mu$m]"
else:
    x_edges_plot = x_edges
    y_edges_plot = y_edges
    xlabel = r"$x$ [cm]"
    ylabel = r"$y$ [cm]"

arr_plot, x_plot_edges, y_plot_edges = downsample_for_pcolormesh(
    arr2d,
    x_edges_plot,
    y_edges_plot,
    sx=sx,
    sy=sy,
)

arr_plot_2, x_plot_edges, y_plot_edges = downsample_for_pcolormesh(
    arr2d_2,
    x_edges_plot,
    y_edges_plot,
    sx=sx,
    sy=sy,
)
# Convert the downsampled view to a normal array.
# This keeps only the smaller plotted array in memory.
arr_plot = np.asarray(arr_plot)
arr_plot_2 = np.asarray(arr_plot_2)

print("downsampled shape:", arr_plot.shape)

norm = choose_norm(
    arr_plot,
    useLog=useLog,
    vmin_manual=vmin_manual,
    vmax_manual=vmax_manual,
)

fig, ax = plt.subplots(figsize=(7, 5))
# y = np.linspace(ymin, ymax, arr_plot.shape[1])
# ax.plot(y*10**4,arr_plot[512,:],lw = 2, color = "b", label=r"$n_{i}$")
# ax.plot(y*10**4,arr_plot_2[512,:],lw = 2, color = "r", label=r"$n_{e}$")
# plt.yscale("log")
# ax.set_xlabel(xlabel)
# ax.set_ylabel(r"n [$cm^{-3}$]")
# ax.set_title(f"FLASH 2D Density, t = 0.185 ns")
# ax.grid(True, alpha=0.3)
# plt.legend()
# plt.tight_layout()
# plt.show()
pcm = ax.pcolormesh(
    x_plot_edges,
    y_plot_edges,
    arr_plot.T,
    shading="auto",
    cmap="plasma",
    norm=norm,
)

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_title(field_name)
ax.set_aspect("auto")

cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label(r"$n_{i}$")

plt.tight_layout()
plt.show()

