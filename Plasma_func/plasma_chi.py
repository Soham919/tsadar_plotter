import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy.special import wofz
from Plasma_func import v_th, lam_db

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h


import numpy as np
from scipy import constants
from scipy.special import wofz


EPS0 = constants.epsilon_0
E_CHARGE = constants.e
M_E = constants.electron_mass


def plasma_dispersion(zeta):
    """
    Plasma dispersion function:

        Z(zeta) = i sqrt(pi) w(zeta)
    """
    zeta = np.asarray(zeta, dtype=complex)
    return 1j * np.sqrt(np.pi) * wofz(zeta)


class PlasmaSpecies:
    """
    Properties and kinetic response of one Maxwellian plasma species.

    Parameters
    ----------
    name : str
        Species name.

    Z : float
        Signed charge state:
            electron: -1
            H+:        1
            He2+:      2

    m : float
        Particle mass [kg].

    T_eV : float
        Temperature [eV].

    n_cm3 : float
        Number density [cm^-3].

    ud : float
        Drift velocity projected along k [m/s].
    """

    def __init__(self, name, Z, m, T_eV, n_cm3, ud=0.0):
        if m <= 0:
            raise ValueError("Species mass must be positive.")
        if T_eV <= 0:
            raise ValueError("Species temperature must be positive.")
        if n_cm3 < 0:
            raise ValueError("Species density cannot be negative.")
        if Z == 0:
            raise ValueError("Charge state Z cannot be zero.")

        self.name = name
        self.Z = float(Z)
        self.m = float(m)
        self.T_eV = float(T_eV)
        self.n_cm3 = float(n_cm3)
        self.ud = float(ud)

    @property
    def n_m3(self):
        """Number density [m^-3]."""
        return self.n_cm3 * 1.0e6

    def thermal_speed(self):
        """
        Thermal speed:

            v_t = sqrt(2 T_eV e / m)
        """
        return np.sqrt(
            2.0 * self.T_eV * E_CHARGE / self.m
        )

    def debye_length(self):
        """
        Species Debye length:

            lambda_D^2 = eps0 T_eV e / (n Z^2 e^2)
        """
        if self.n_m3 == 0:
            return np.inf

        charge = abs(self.Z) * E_CHARGE

        return np.sqrt(
            EPS0 * self.T_eV * E_CHARGE
            / (self.n_m3 * charge**2)
        )

    def maxwellian_1d(self, phase_velocity):
        """
        Normalized 1D Maxwellian:

            f(v) = exp[-((v-u)/v_t)^2] / (sqrt(pi) v_t)
        """
        phase_velocity = np.asarray(phase_velocity)
        vt = self.thermal_speed()

        return (
            np.exp(
                -((phase_velocity - self.ud) / vt) ** 2
            )
            / (np.sqrt(np.pi) * vt)
        )

    def susceptibility(self, omega, k):
        """
        Collisionless Maxwellian susceptibility:

            chi_s = 1/(k^2 lambda_Ds^2)
                    [1 + zeta_s Z(zeta_s)]

            zeta_s = (omega - k u_s)/(k v_ts)
        """
        if k == 0:
            raise ValueError("Wave number k cannot be zero.")

        omega = np.asarray(omega)
        vt = self.thermal_speed()
        lambda_D = self.debye_length()

        zeta = (omega - k * self.ud) / (k * vt)

        chi = (1.0 / (k**2 * lambda_D**2)) * (1.0 + zeta * plasma_dispersion(zeta))

        return chi


def thomson_structure_factor(
    k,
    phase_velocity,
    electron,
    ions,
    neutrality_rtol=1.0e-6,
):
    """
    Collisionless Maxwellian Thomson form factor for one electron
    species and an arbitrary number of ion species.

    Ion contributions are weighted using ion number fractions:

        x_i = n_i / sum_j(n_j)

    so that

        Z_i^2 n_i / n_e = Z_i^2 x_i / Z_bar

    where

        Z_bar = sum_i(Z_i x_i).

    Parameters
    ----------
    k : float
        Scattering wave number [m^-1].

    phase_velocity : array-like
        Phase velocity omega/k [m/s].

    electron : PlasmaSpecies
        Electron species.

    ions : iterable of PlasmaSpecies
        Ion species in the plasma mixture.

    neutrality_rtol : float
        Relative tolerance used when checking quasineutrality.

    Returns
    -------
    results : dict
        Dictionary containing the total form factor, individual
        species contributions, susceptibilities, dielectric function,
        ion fractions, and average ionization.
    """

    if k == 0:
        raise ValueError("Wave number k cannot be zero.")

    phase_velocity = np.asarray(
        phase_velocity,
        dtype=float,
    )

    ions = list(ions)

    if len(ions) == 0:
        raise ValueError(
            "At least one ion species is required."
        )

    if electron.n_m3 <= 0:
        raise ValueError(
            "Electron density must be positive."
        )

    for ion in ions:
        if ion.Z <= 0:
            raise ValueError(
                f"Ion species {ion.name!r} must have Z > 0."
            )

    omega = k * phase_velocity

    # ========================================================
    # Ion mixture properties
    # ========================================================

    total_ion_density = sum(ion.n_m3 for ion in ions)

    if total_ion_density <= 0:
        raise ValueError(
            "Total ion density must be positive."
        )

    ion_fractions = {ion.name: ion.n_m3 / total_ion_density for ion in ions}

    fraction_sum = sum(ion_fractions.values())

    if not np.isclose(fraction_sum, 1.0):
        raise RuntimeError(
            "Calculated ion fractions do not sum to one."
        )

    Z_bar = sum(ion.Z * ion_fractions[ion.name] for ion in ions)

    # Electron density implied by the ion mixture
    quasineutral_ne = total_ion_density * Z_bar

    if not np.isclose(electron.n_m3, quasineutral_ne,rtol=neutrality_rtol,):
        raise ValueError(
            "Electron and ion densities do not satisfy "
            "quasineutrality:\n"
            f"  supplied ne      = {electron.n_m3:.6e} m^-3\n"
            f"  sum(Z_i n_i)     = {quasineutral_ne:.6e} m^-3"
        )

    # ========================================================
    # Susceptibilities and dielectric function
    # ========================================================

    chi_e = electron.susceptibility(omega=omega,k=k,)

    chi_ions = {ion.name: ion.susceptibility(omega=omega,k=k,) for ion in ions}

    epsilon = (1.0 + chi_e + sum(chi_ions.values()))

    # ========================================================
    # Electron contribution
    # ========================================================

    f_e = electron.maxwellian_1d(phase_velocity)

    electron_screening = np.abs(1.0 - chi_e / epsilon) ** 2

    ion_screening = np.abs(chi_e / epsilon) ** 2

    prefactor = 2.0 * np.pi / abs(k)

    S_e = (prefactor * electron_screening * f_e)

    # ========================================================
    # Ion contributions
    # ========================================================

    S_ions = {}
    ion_weights = {}

    for ion in ions:
        fraction = ion_fractions[ion.name]

        # Equivalent to Z_i^2 n_i / n_e
        weight = (ion.Z**2 * fraction / Z_bar)

        ion_weights[ion.name] = weight

        f_i = ion.maxwellian_1d(phase_velocity)

        S_ions[ion.name] = (prefactor * ion_screening * weight * f_i)

    S_total = (S_e + sum(S_ions.values()))

    return {
        "S_total": S_total,
        "S_e": S_e,
        "S_ions": S_ions,
        "epsilon": epsilon,
        "chi_e": chi_e,
        "chi_ions": chi_ions,
        "omega": omega,
        "phase_velocity": phase_velocity,
        "ion_fractions": ion_fractions,
        "ion_weights": ion_weights,
        "total_ion_density": total_ion_density,
        "Z_bar": Z_bar,
        "quasineutral_ne": quasineutral_ne,
    }

if __name__ == "__main__":
    # ============================================================
    # Example H-He plasma
    # ============================================================

    # Probe geometry
    probe_wavelength_nm = 526.5
    scattering_angle_deg = 60.0

    probe_wavelength_m = probe_wavelength_nm * 1e-9
    scattering_angle_rad = np.deg2rad(scattering_angle_deg)

    # Scattering wave number:
    #
    #     k = |k_s - k_0|
    #       = (4 pi / lambda_0) sin(theta / 2)
    #
    k = (
        4.0 * np.pi / probe_wavelength_m
        * np.sin(scattering_angle_rad / 2.0)
    )

    # Plasma conditions
    ne_cm3 = 2.19e19  # cm-3
    ne = ne_cm3 * 1.0e6  # cm^-3 -> m^-3

    Te_eV = 245.0  # eV
    TH_eV = 614.0
    THe_eV = 386.0

    ud_e = 0.0  # 10^6 cm/s
    ud_H = 28.52
    ud_He = 32.3
    # Ion charge states
    ZH = 1.0    
    ZHe = 2.0

    # Ion number fractions:
    #
    #     f_H + f_He = 1
    #
    fH = 0.55
    fHe = 0.45
    Zbar = ZH * fH + ZHe * fHe

    if not np.isclose(fH + fHe, 1.0):
        raise ValueError("Ion number fractions must sum to one.")

    electron = PlasmaSpecies(name = 'electron', Z = -1, m = me, T_eV = Te_eV, n_cm3 = ne_cm3, ud = ud_e*(10**4))
    hydrogen = PlasmaSpecies(name = 'hydrogen', Z = ZH, m = mp, T_eV = TH_eV, n_cm3 = fH * (ne_cm3 / Zbar), ud = ud_H*(10**4))
    helium = PlasmaSpecies(name = 'helium', Z = ZHe, m = 4.0 * mp, T_eV = THe_eV, n_cm3 = fHe * (ne_cm3 / Zbar), ud = ud_He*(10**4))
    ions  = [hydrogen, helium]

    # ============================================================
    # Frequency axis
    # ============================================================

    # Use phase velocity as the plotting variable first.
    #
    # Since omega = k*v_phase, this makes the ion resonances easier
    # to interpret physically.
    phase_velocity = np.linspace(
        -7.0e7,
        7.0e7,
        10000,
    )

    # ============================================================
    # Calculate everything
    # ============================================================

    Results = thomson_structure_factor(k=k, phase_velocity = phase_velocity, electron = electron, ions = ions)

    S_total = Results["S_total"]
    S_e = Results["S_e"]
    S_ions = Results["S_ions"]
    epsilon = Results["epsilon"]
    chi_e = Results["chi_e"]
    chi_ions = Results["chi_ions"]
    omega = Results["omega"]
    phase_velocity = Results["phase_velocity"]
    ion_fractions = Results["ion_fractions"]
    ion_weights = Results["ion_weights"]
    total_ion_density = Results["total_ion_density"]
    Z_bar = Results["Z_bar"]

    # ============================================================
    # Print useful scales
    # ============================================================

    print(f"Scattering k = {k:.6e} m^-1")
    print()

    print(f"ne  = {electron.n_cm3:.6e} cm^-3")
    print(f"nH  = {hydrogen.n_cm3:.6e} cm^-3")
    print(f"nHe = {helium.n_cm3:.6e} cm^-3")
    print()

    print(
        f"Electron thermal speed = "
        f"{electron.thermal_speed():.6e} m/s"
    )

    print(
        f"H thermal speed        = "
        f"{hydrogen.thermal_speed():.6e} m/s"
    )

    print(
        f"He thermal speed       = "
        f"{helium.thermal_speed():.6e} m/s"
    )


    # ============================================================
    # Plot real parts
    # ============================================================

    # fig, ax = plt.subplots(figsize=(9, 5))

    # ax.plot(
    #     phase_velocity / 1e3,
    #     chi_e.real,
    #     label=r"$\mathrm{Re}(\chi_e)$",
    # )

    # ax.plot(
    #     phase_velocity / 1e3,
    #     chi_H.real,
    #     label=r"$\mathrm{Re}(\chi_H)$",
    # )

    # ax.plot(
    #     phase_velocity / 1e3,
    #     chi_He.real,
    #     label=r"$\mathrm{Re}(\chi_{\mathrm{He}})$",
    # )

    # ax.axhline(0.0, linewidth=0.8)

    # ax.set_xlabel(
    #     r"Phase velocity $\omega/k$ [km/s]"
    # )
    # ax.set_ylabel(r"$\mathrm{Re}(\chi_s)$")
    # ax.set_title("Real parts of species susceptibilities")
    # ax.grid(True, alpha=0.3)
    # ax.legend()

    # plt.show()


    # ============================================================
    # Plot imaginary parts
    # ============================================================

    # fig, ax = plt.subplots(figsize=(9, 5))

    # ax.plot(
    #     phase_velocity / 1e3,
    #     chi_e.imag,
    #     label=r"$\mathrm{Im}(\chi_e)$",
    # )

    # ax.plot(
    #     phase_velocity / 1e3,
    #     chi_H.imag,
    #     label=r"$\mathrm{Im}(\chi_H)$",
    # )

    # ax.plot(
    #     phase_velocity / 1e3,
    #     chi_He.imag,
    #     label=r"$\mathrm{Im}(\chi_{\mathrm{He}})$",
    # )

    # ax.axhline(0.0, linewidth=0.8)

    # ax.set_xlabel(
    #     r"Phase velocity $\omega/k$ [km/s]"
    # )
    # ax.set_ylabel(r"$\mathrm{Im}(\chi_s)$")
    # ax.set_title("Imaginary parts of species susceptibilities")
    # ax.grid(True, alpha=0.3)
    # ax.legend()

    # plt.show()


    # ============================================================
    # Find zero crossings of Re(epsilon)
    # ============================================================

    v_kms = phase_velocity / 1e3

    crossing_indices = np.where(
        np.diff(np.sign(epsilon.real)) != 0
    )[0]

    zero_crossings = []

    for i in crossing_indices:
        # Linear interpolation between adjacent points
        v1 = phase_velocity[i]
        v2 = phase_velocity[i + 1]

        eps1 = epsilon.real[i]
        eps2 = epsilon.real[i + 1]

        v_zero = v1 - eps1 * (v2 - v1) / (eps2 - eps1)
        zero_crossings.append(v_zero)

    print("\nZero crossings of Re(epsilon):")

    for v_zero in zero_crossings:
        nearest_index = np.argmin(
            np.abs(phase_velocity - v_zero)
        )

        print(
            f"  v_phase = {v_zero / 1e3:8.3f} km/s, "
            f"Im(epsilon) = "
            f"{epsilon.imag[nearest_index]:+.6e}"
        )


    # ============================================================
    # 2 x 2 diagnostic panel
    # ============================================================

    ion_mask = np.abs(phase_velocity) <= 1.0e6

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(15, 10),
        constrained_layout=True,
    )

    ax_eps = axes[0, 0]
    ax_full = axes[0, 1]
    ax_chi = axes[1, 0]
    ax_ion = axes[1, 1]


    # ------------------------------------------------------------
    # Helper: add Re(epsilon)=0 markers
    # ------------------------------------------------------------
    def add_mode_markers(ax, crossings):
        """
        Draw vertical lines at zero crossings of Re(epsilon).
        """
        for v_zero in crossings:
            ax.axvline(
                v_zero / 1e3,
                color="gray",
                linestyle=":",
                linewidth=1.2,
                alpha=0.8,
            )


    # ------------------------------------------------------------
    # (1,1) Dielectric function
    # ------------------------------------------------------------
    ax_eps.plot(v_kms, epsilon.real, linewidth=2, label=r"$\mathrm{Re}(\epsilon)$")

    ax_eps.plot(v_kms, epsilon.imag, linewidth=2, label=r"$\mathrm{Im}(\epsilon)$")

    ax_eps.axhline(0.0, color="k", linewidth=0.8, linestyle="--")

    add_mode_markers(ax_eps, zero_crossings)

    ax_eps.set_xlabel(
        r"Phase velocity $\omega/k$ [km/s]"
    )
    ax_eps.set_ylabel(
        r"$\epsilon(k,\omega)$"
    )
    ax_eps.set_title(
        r"Dielectric function: "
        r"$\epsilon=1+\chi_e+\chi_H+\chi_{\mathrm{He}}$"
    )

    ax_eps.grid(True, alpha=0.3)
    ax_eps.legend()


    # ------------------------------------------------------------
    # (1,2) Full Thomson spectrum
    # ------------------------------------------------------------
    ax_full.plot(v_kms,S_total,linewidth=2.2,label=r"$S_{\mathrm{total}}$",)

    ax_full.plot(v_kms,S_e, linestyle="--", color = "green", linewidth=1.8, label=r"$S_e$")

    ax_full.plot(v_kms, S_ions["hydrogen"], linestyle=":", linewidth=1.8, color = "red",label=r"$S_H$")

    ax_full.plot(v_kms, S_ions["helium"], linestyle="-.", linewidth=1.8, color = "blue", label=r"$S_{\mathrm{He}}$")

    add_mode_markers(ax_full, zero_crossings)

    ax_full.set_xlabel(
        r"Phase velocity $\omega/k$ [km/s]"
    )
    ax_full.set_ylabel(
        r"$S(k,\omega)$ [s]"
    )
    ax_full.set_title(
        "Full Thomson-scattering spectrum"
    )
    ax_full.set_ylim([0,0.1e-13])
    ax_full.grid(True, alpha=0.3)
    ax_full.legend()


    # ------------------------------------------------------------
    # (2,1) Real and imaginary susceptibilities
    # ------------------------------------------------------------

    # Electrons
    ax_chi.plot(v_kms, chi_e.real, linewidth=2, color = "green", label=r"$\mathrm{Re}(\chi_e)$")

    ax_chi.plot(v_kms, chi_e.imag, linewidth=1.7, linestyle="--", color = "green", label=r"$\mathrm{Im}(\chi_e)$")

    # Hydrogen
    ax_chi.plot(v_kms, chi_ions['hydrogen'].real, linewidth=2, color = "red", label=r"$\mathrm{Re}(\chi_H)$")

    ax_chi.plot(v_kms, chi_ions['hydrogen'].imag, linewidth=1.7, linestyle="--", color = "red", label=r"$\mathrm{Im}(\chi_H)$")

    # Helium
    ax_chi.plot(v_kms, chi_ions['helium'].real, linewidth=2, color = "blue", label=r"$\mathrm{Re}(\chi_{\mathrm{He}})$",)

    ax_chi.plot(v_kms, chi_ions['helium'].imag, linewidth=1.7, linestyle="--", color = "blue", label=r"$\mathrm{Im}(\chi_{\mathrm{He}})$",)

    ax_chi.axhline( 0.0, color="k", linewidth=0.8, linestyle="--")

    add_mode_markers(ax_chi, zero_crossings)

    ax_chi.set_xlabel(
        r"Phase velocity $\omega/k$ [km/s]"
    )
    ax_chi.set_ylabel(
        r"$\chi_s(k,\omega)$"
    )
    ax_chi.set_title(
        "Species susceptibilities"
    )

    ax_chi.grid(True, alpha=0.3)
    ax_chi.legend(
        fontsize=9,
        ncol=2,
    )


    # ------------------------------------------------------------
    # (2,2) Zoomed ion feature
    # ------------------------------------------------------------
    ax_ion.plot(v_kms[ion_mask], S_total[ion_mask],
        linewidth=2.2,
        label=r"$S_{\mathrm{total}}$",
    )

    ax_ion.plot(v_kms[ion_mask], S_ions['hydrogen'][ion_mask],
        linestyle=":",
        color = "red",
        linewidth=1.8,
        label=r"$S_H$ contribution",
    )

    ax_ion.plot(v_kms[ion_mask], S_ions['helium'][ion_mask],
        linestyle="-.",
        color = "blue",
        linewidth=1.8,
        label=r"$S_{\mathrm{He}}$ contribution",
    )

    # Only show zero crossings inside the ion window
    ion_zero_crossings = [
        v_zero
        for v_zero in zero_crossings
        if abs(v_zero) <= 4.0e5
    ]

    add_mode_markers(
        ax_ion,
        ion_zero_crossings,
    )

    ax_ion.set_xlabel(
        r"Phase velocity $\omega/k$ [km/s]"
    )
    ax_ion.set_ylabel(
        r"$S(k,\omega)$ [s]"
    )
    ax_ion.set_title(
        "Ion feature"
    )

    ax_ion.grid(True, alpha=0.3)
    ax_ion.legend()

    plt.show()

    fig2, ax2 = plt.subplots(figsize=(9, 5))
    ax2.plot(v_kms, 1/((epsilon.real)**2 + (epsilon.imag)**2) , linewidth=2, label=r"$\mathrm{Re}(\epsilon)$")
    plt.show()