import numpy as np
from dataclasses import dataclass
from scipy.constants import epsilon_0, elementary_charge, m_p
import matplotlib.pyplot as plt


@dataclass
class IonSpecies:
    name: str
    A: float          # Atomic mass in proton-mass units
    Z: float          # Ion charge state
    fraction: float   # Ion number fraction; fractions must sum to 1


def electron_debye_length(Te_eV: float, ne_m3: float) -> float:
    """
    Electron Debye length in meters.

    Te_eV : electron temperature [eV]
    ne_m3 : electron density [m^-3]
    """
    Te_J = Te_eV * elementary_charge

    return np.sqrt(
        epsilon_0 * Te_J
        / (ne_m3 * elementary_charge**2)
    )


def weighted_sum(species, quantity):
    """
    Compute <X> = sum_beta f_beta X_beta.
    """
    return sum(ion.fraction * quantity(ion) for ion in species)


def eq12_phase_speed(
    Te_eV: float,
    Ti_eV: float,
    k_m_inv: float,
    ne_m3: float,
    heavy_species,
    light_species,
    all_species,
):
    """
    Evaluate Eq. (12) of Williams et al. (1995).

    Returns
    -------
    v_phase : float
        Approximate real phase velocity omega_r/k [m/s].
    details : dict
        Intermediate quantities used in the calculation.
    """

    if Te_eV <= 0 or Ti_eV <= 0:
        raise ValueError("Te_eV and Ti_eV must be positive.")

    fraction_sum = sum(ion.fraction for ion in all_species)

    if not np.isclose(fraction_sum, 1.0):
        raise ValueError(
            f"Ion fractions must sum to 1; received {fraction_sum:.6g}."
        )

    # Average charge: Z_bar = sum f_beta Z_beta
    Z_bar = weighted_sum(all_species, lambda ion: ion.Z)

    # Heavy-species weighted sums
    Z2_over_A_h = weighted_sum(
        heavy_species,
        lambda ion: ion.Z**2 / ion.A,
    )

    Z2_over_A2_h = weighted_sum(
        heavy_species,
        lambda ion: ion.Z**2 / ion.A**2,
    )

    # Light-species contribution to Debye screening
    Z2_l = weighted_sum(
        light_species,
        lambda ion: ion.Z**2,
    )

    mu = (Z2_l / Z_bar) * (Te_eV / Ti_eV)

    lambda_D = electron_debye_length(Te_eV, ne_m3)
    k_lambda_D_sq = (k_m_inv * lambda_D)**2

    Te_J = Te_eV * elementary_charge
    Ti_J = Ti_eV * elementary_charge

    # Electron-pressure contribution
    electron_term = (
        (Z2_over_A_h / Z_bar)
        * (Te_J / m_p)
        / (1.0 + mu + k_lambda_D_sq)
    )

    # Ion-pressure correction
    ion_term = (
        3.0
        * (Z2_over_A2_h / Z2_over_A_h)
        * (Ti_J / m_p)
    )

    v_phase_sq = electron_term + ion_term

    if v_phase_sq < 0:
        raise ValueError("Equation produced a negative phase-speed squared.")

    return np.sqrt(v_phase_sq), {
        "Z_bar": Z_bar,
        "mu": mu,
        "lambda_D_m": lambda_D,
        "k_lambda_D": k_m_inv * lambda_D,
        "electron_term_m2_s2": electron_term,
        "ion_term_m2_s2": ion_term,
    }


def binary_fast_slow_modes(
    ion1: IonSpecies,
    ion2: IonSpecies,
    Te_eV: float,
    Ti_eV: float,
    k_m_inv: float,
    ne_m3: float,
):
    """
    Calculate the approximate fast and slow binary-ion modes from Eq. (12).

    Fast mode:
        Both ion species are treated as 'heavy' and mu = 0.

    Slow mode:
        The species with the smaller thermal velocity is placed in the
        heavy group. The faster-thermal-speed species is placed in the
        light group.
    """

    ions = [ion1, ion2]

    # For a common Ti, v_th,beta = sqrt(Ti / (A_beta M_p)).
    # The larger-A species therefore has the smaller thermal speed.
    thermal_speeds = {
        ion.name: np.sqrt(
            Ti_eV * elementary_charge / (ion.A * m_p)
        )
        for ion in ions
    }

    slower_species = min(
        ions,
        key=lambda ion: thermal_speeds[ion.name],
    )

    faster_species = max(
        ions,
        key=lambda ion: thermal_speeds[ion.name],
    )

    # Fast mode: all ions are in the heavy group
    v_fast, fast_details = eq12_phase_speed(
        Te_eV=Te_eV,
        Ti_eV=Ti_eV,
        k_m_inv=k_m_inv,
        ne_m3=ne_m3,
        heavy_species=ions,
        light_species=[],
        all_species=ions,
    )

    # Slow mode: slow-thermal-speed ion is heavy;
    # fast-thermal-speed ion contributes to screening.
    v_slow, slow_details = eq12_phase_speed(
        Te_eV=Te_eV,
        Ti_eV=Ti_eV,
        k_m_inv=k_m_inv,
        ne_m3=ne_m3,
        heavy_species=[slower_species],
        light_species=[faster_species],
        all_species=ions,
    )

    return {
        "fast_mode_m_s": v_fast,
        "slow_mode_m_s": v_slow,
        "fast_mode_cm_s": 100.0 * v_fast,
        "slow_mode_cm_s": 100.0 * v_slow,
        "fast_details": fast_details,
        "slow_details": slow_details,
        "slow_mode_heavy_species": slower_species.name,
        "slow_mode_light_species": faster_species.name,
        "thermal_speeds_m_s": thermal_speeds,
    }

if __name__ == "__main__":
    # ============================================================
    # Example: fully ionized H-He plasma
    # ============================================================

    Te_eV = 200.0
    Ti_eV = 200.0

    # Example values — replace with your Thomson-scattering conditions
    ne_m3 = 2.0e25
    k_m_inv = 11933875.227311652
    #fr = np.linspace(0.01, 0.99, 50)  # H fraction
    fH = 0.9
    fHe = 1 - fH
    Te = 500
    Ti = np.linspace(50,1000,1000)
    fast = []
    slow = []
    # for f in fr:
    #     H = IonSpecies(name="H+", A=1.0, Z=1.0, fraction=f,)
    #     He = IonSpecies(name="He2+",A=12.0,Z=4.0,fraction=1 - f,)
    for i in range(len(Ti)):
        H = IonSpecies(name="H+", A=1.0, Z=1.0, fraction=fH,)
        He = IonSpecies(name="He2+",A=14.0,Z=3.0,fraction=fHe,)
        result = binary_fast_slow_modes(
            ion1=H,
            ion2=He,
            Te_eV=Te,
            Ti_eV=Ti[i],
            k_m_inv=k_m_inv,
            ne_m3=ne_m3,
        )

        print(f"Fast mode: {result['fast_mode_cm_s']:.4e} cm/s")
        print(f"Slow mode: {result['slow_mode_cm_s']:.4e} cm/s")

        fast = np.append(fast, result['fast_mode_cm_s'])
        slow = np.append(slow, result['slow_mode_cm_s'])

        # print("Slow-mode groups: " f"heavy = {result['slow_mode_heavy_species']}, "f"light = {result['slow_mode_light_species']}")

        # print(f"Fast-mode mu: {result['fast_details']['mu']:.4g}")
        # print(f"Slow-mode mu: {result['slow_details']['mu']:.4g}")
        # print("k lambda_D: " f"{result['fast_details']['k_lambda_D']:.4g}")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(Ti, fast/10**6, label="Fast mode", color="red")
    ax.plot(Ti, slow/10**6, label="Slow mode", color="blue")
    ax.set_xlabel("Ti (eV)")
    ax.set_ylabel("Phase speed [10^6 cm/s]")
    ax.legend()
    plt.show()