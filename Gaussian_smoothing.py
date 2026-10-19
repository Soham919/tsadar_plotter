import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def mandelbrot(h, w, maxit=20, r=2):
    """Returns an image of the Mandelbrot fractal of size (h,w)."""
    x = np.linspace(-2.5, 1.5, 4*h+1)
    y = np.linspace(-1.5, 1.5, 3*w+1)
    A, B = np.meshgrid(x, y)
    C = A + B*1j
    z = np.zeros_like(C)
    divtime = maxit + np.zeros(z.shape, dtype=np.int_)

    for i in range(maxit):
        z = z**2 + C
        diverge = abs(z) > r                    # who is diverging
        div_now = diverge & (divtime == maxit)  # who is diverging now
        divtime[div_now] = i                    # note when
        z[diverge] = r                          # avoid diverging too much

    return divtime

# Your 2D array
arr = mandelbrot(400, 400)

# Gaussian smoothing
sigma = 5.0
arr_smooth = gaussian_filter(arr, sigma=sigma)

# Use same color limits so comparison is meaningful
vmin = np.nanmin(arr)
vmax = np.nanmax(arr)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

im1 = ax[0].imshow(arr, origin="lower", cmap="viridis",
                  vmin=vmin, vmax=vmax, aspect="auto")
ax[0].set_title("Original")

im2 = ax[1].imshow(arr_smooth, origin="lower", cmap="viridis",
                  vmin=vmin, vmax=vmax, aspect="auto")
ax[1].set_title(rf"Gaussian smoothed, $\sigma={sigma}$")

fig.colorbar(im2, ax=ax, label="Value")

plt.show()