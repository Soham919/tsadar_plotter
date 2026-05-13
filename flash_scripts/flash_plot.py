import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
import h5_helpers.h5_helper as h5h
from pathlib import Path
import yt

baseDir = Path().resolve().parent
# --------- Mac ------------
fp = baseDir/".."/"Flash"/"test_runs"/"1D_Al_1um"
file = "ks_hdf5_plt_cnt_0001"
fp = fp / file

# --------- Windows ------------
# shot = "117827"  # change this for the shot you want to plot
# fp = baseDir/"Kinshock"/"Kinshock-26A"/"data"/shot
# saveDir = baseDir/"Kinshock"/"Kinshock-26A"/"data"/shot/"plots"
# if not saveDir.exists():
#     saveDir.mkdir(parents = True, exist_ok = True)


# ---------- Using yt, works --------------- #
ds = yt.load(fp)
print("max level:", ds.index.max_level)
ad = ds.all_data()

# Get simulation time
sim_time = float(ds.current_time)      # usually in code time / seconds depending on setup
sim_time_ns = sim_time * 1e9

print("current_time =", ds.current_time)
print(f"time = {sim_time_ns:.4f} ns")

## print fields ##
print("Native fields:")
for field in ds.field_list:
    try:
        print("  ",f"{field}: {ad[field].units}")
    except Exception as e:
        print(f"{field}: failed ({e})")
    
print("\nDerived fields:")
for field in ds.derived_field_list:
    try:
        print("  ",f"{field}: {ad[field].units}")
    except Exception as e:
        print(f"{field}: failed ({type(e).__name__}: {e})")

level = ds.index.max_level
dims = ds.domain_dimensions * (2**level)

print("covering grid dims:", dims)

cg = ds.covering_grid(
    level=level,
    left_edge=ds.domain_left_edge,
    dims=dims
)

temp = cg[("flash", "tion")].to_ndarray()

print("shape from covering_grid:", temp.shape)
#print("min/max:", np.min(temp), np.max(temp))

# Pick one z-slice explicitly
temp2d = temp[:, :, 0]

#print("temp2d shape:", temp2d.shape)
#print("temp2d min/max:", np.min(temp2d), np.max(temp2d))

# Domain edges
x0, x1 = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
y0, y1 = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

nx, ny = temp2d.shape

# Cell edges for pcolormesh
x_edges = np.linspace(x0, x1, nx + 1)
y_edges = np.linspace(y0, y1, ny + 1)

X, Y = np.meshgrid(x_edges, y_edges, indexing="ij")

positive = temp2d[temp2d > 0]

fig, ax = plt.subplots(figsize=(6, 6))

pcm = ax.pcolormesh(
    X, Y, temp2d,
    shading="auto",
    cmap="plasma",
    norm=LogNorm(vmin=np.min(positive), vmax=np.max(temp2d))
)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Electron Number Density, t = {sim_time_ns:.3f} ns")
ax.set_aspect("equal")

cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label(r"$n_{e}$(cm^-3)")

plt.tight_layout()
plt.show()
# ----------------------------------- #

# ------- This is the normal h5py part, doesn't work rn -------- #
# h5h.h5_show(fp)

# with h5py.File(fp, "r") as f:
#     pres = f["pres"][:]            # shape: (nblocks, 1, ny, nx)
#     bbox = f["bounding box"][:]    # shape: (nblocks, ndim, 2)

# print("pres shape:", pres.shape)
# print("bbox shape:", bbox.shape)
# print("global min/max:", np.nanmin(pres), np.nanmax(pres))

# # For log scale, only positive values are allowed
# positive_vals = pres[pres > 0]

# if positive_vals.size == 0:
#     raise ValueError("No positive values found in 'pres'; cannot use LogNorm.")

# vmin = np.nanmin(positive_vals)
# vmax = np.nanmax(pres)

# fig, ax = plt.subplots(figsize=(6, 6))
# im = None

# for i in range(pres.shape[0]):
#     block = pres[i, 0, :, :].T   # try .T here if orientation looks wrong

#     xmin, xmax = bbox[i, 0, 0], bbox[i, 0, 1]
#     ymin, ymax = bbox[i, 1, 0], bbox[i, 1, 1]

#     im = ax.imshow(
#         block,
#         origin="lower",
#         extent=[xmin, xmax, ymin, ymax],
#         interpolation="nearest",
#         aspect="equal",
#         norm=LogNorm(vmin=vmin, vmax=vmax),
#         cmap="viridis"
#     )

# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_title("Pressure (log scale)")

# cbar = fig.colorbar(im, ax=ax)
# cbar.set_label("Pressure")

# plt.tight_layout()
# plt.show()