import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz


def plasma_dispersion(zeta):
    """
    Plasma dispersion function

        Z(zeta) = i * sqrt(pi) * w(zeta)

    where w(zeta) is the Faddeeva function.

    Parameters
    ----------
    zeta : complex, float, or array-like
        Dimensionless phase velocity.

    Returns
    -------
    complex or np.ndarray
        Complex plasma dispersion function.
    """
    zeta = np.asarray(zeta, dtype=complex)
    return 1j * np.sqrt(np.pi) * wofz(zeta)


# -------------------------------------------------
# Evaluate Z(zeta) along the real zeta axis
# -------------------------------------------------
zeta = np.linspace(-5.0, 5.0, 2000)
Z = plasma_dispersion(zeta)
dZ_dzeta = np.gradient(Z, zeta)  # Numerical derivative of Z with respect to zeta

# -------------------------------------------------
# Plot real and imaginary parts
# -------------------------------------------------
fig, ax = plt.subplots(2,1,figsize=(8, 10))

ax[0].plot(zeta, Z.real, label=r"$\mathrm{Re}[Z(\zeta)]$")
ax[0].plot(zeta, Z.imag, label=r"$\mathrm{Im}[Z(\zeta)]$")

ax[0].axhline(0.0, linewidth=0.8)
ax[0].axvline(0.0, linewidth=0.8)

ax[0].set_xlabel(r"$\zeta$")
ax[0].set_ylabel(r"$Z(\zeta)$")
ax[0].set_title("Plasma dispersion function")
ax[0].grid(True, alpha=0.3)
ax[0].legend()

ax[1].plot(zeta, dZ_dzeta.real, label=r"$\mathrm{Re}[Z(\zeta)]$")
ax[1].plot(zeta, dZ_dzeta.imag, label=r"$\mathrm{Im}[Z(\zeta)]$")

ax[1].axhline(0.0, linewidth=0.8)
ax[1].axvline(0.0, linewidth=0.8)

ax[1].set_xlabel(r"$\zeta$")
ax[1].set_ylabel(r"$Z'(\zeta)$")
ax[1].set_title("Derivative of Plasma dispersion function")
ax[1].grid(True, alpha=0.3)
ax[1].legend()

plt.show()