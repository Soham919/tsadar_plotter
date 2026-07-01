import numpy as np
import matplotlib.pyplot as plt

# ---------------- user inputs ----------------
outfile = "curved_target_table.txt"

nx = 10000
ny = 10000

x_min, x_max = -0.05, 0.05   # cm
y_min, y_max = -0.2, 0.8   # cm

y0 = 0.005        # cm, target center at x=0
thickness = 1e-4  # cm, 1 micron target thickness

curvature_amp = 5e-4  # cm, max bow/curve amount
x_half_width = 0.05   # cm, used to normalize curvature

rho_target = 0.113    # g/cm^3 or whatever FLASH expects
rho_bg = 1e-10
# ---------------------------------------------

x = np.linspace(x_min, x_max, nx)
y = np.linspace(y_min, y_max, ny)

with open(outfile, "w") as f:
    f.write("# Curved target table\n")
    f.write("# Columns: x_cm y_cm targ rho\n")
    f.write(f"# nx = {nx}, ny = {ny}\n")
    f.write(f"# bounds: x [{x_min}, {x_max}], y [{y_min}, {y_max}] cm\n")

    for j in range(ny):
        for i in range(nx):
            xx = x[i]
            yy = y[j]

            yc = y0 + curvature_amp * (xx / x_half_width)**2

            inside_target = abs(yy - yc) <= thickness / 2

            targ = 1.0 if inside_target else 0.0
            rho = rho_target if inside_target else rho_bg

            f.write(f"{xx:.8e} {yy:.8e} {targ:.1f} {rho:.8e}\n")


# Build target mask for plotting

X, Y = np.meshgrid(x, y)

Yc = y0 + curvature_amp * (X / x_half_width)**2

target_mask = np.abs(Y - Yc) <= thickness/2

plt.figure(figsize=(8,6))

plt.imshow(

    target_mask.astype(float),

    extent=[x_min, x_max, y_min, y_max],

    origin='lower',

    aspect='auto'

)

plt.xlabel("x (cm)")

plt.ylabel("y (cm)")

plt.title("Curved Target Geometry")

plt.colorbar(label="Target")

plt.tight_layout()

plt.show()