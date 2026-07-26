import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve


def gaussian(x, mean=0.0, sigma=1.0, amplitude=1.0):
    """
    Return a Gaussian function.

    Parameters
    ----------
    x : np.ndarray
        Coordinate array.
    mean : float
        Center of the Gaussian.
    sigma : float
        Standard deviation.
    amplitude : float
        Peak amplitude.
    """
    if sigma <= 0:
        raise ValueError("sigma must be positive.")

    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def convolve_functions(x, function_1, function_2):
    """
    Numerically convolve two functions sampled on the same uniform grid.

    The factor dx makes the discrete sum approximate the continuous
    convolution integral:

        h(x) = integral f(x') g(x - x') dx'
    """
    if len(x) != len(function_1) or len(x) != len(function_2):
        raise ValueError("x and both functions must have the same length.")

    dx_values = np.diff(x)

    if not np.allclose(dx_values, dx_values[0]):
        raise ValueError("x must be uniformly spaced.")

    dx = dx_values[0]

    convolution = fftconvolve(
        function_1,
        function_2,
        mode="same",
    ) * dx

    return convolution


# ---------------------------------------------------------
# Define the coordinate grid
# ---------------------------------------------------------
x = np.linspace(-15, 15, 3001)

# First Gaussian: imagine this is the true physical spectrum
f = gaussian(x, mean=-10.0, sigma=1.0, amplitude=1.0) + gaussian(x, mean=7.0, sigma=1.0, amplitude=0.5)

# Second Gaussian: imagine this is the instrument response
g = gaussian(x, mean=10.0, sigma=1.0, amplitude=1.0)

# Normalize the instrument response so its total area is 1.
# This prevents convolution from artificially changing the
# total integrated signal.
g /= np.trapezoid(g, x)

# Perform the convolution
f_convolved = convolve_functions(x, f, g)


# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------
fig, axes = plt.subplots(
    1,
    1,
    figsize=(14, 4),
    sharex=True,
)

axes.plot(x, f, label="Original function", color="blue")
axes.set_title("Original function")
axes.set_xlabel("x")
axes.set_ylabel("Amplitude")
axes.grid(alpha=0.3)

axes.plot(x, g, label="Gaussian response", color="orange")
axes.set_title("Gaussian response")
axes.set_xlabel("x")
axes.set_ylabel("Response")
axes.grid(alpha=0.3)

axes.plot(x, f_convolved, label="Convolved result", color="green")
axes.set_title("Convolved result")
axes.set_xlabel("x")
axes.set_ylabel("Amplitude")
axes.grid(alpha=0.3)

plt.tight_layout()
plt.legend()
plt.show()