import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import erf
from scipy import constants
from dataclasses import dataclass
import matplotlib.pyplot as plt

# -----------------------------
# Physical constants, SI units
# -----------------------------
eps0 = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c # speed of light in vacuum (m/s)
h = constants.h # Planck constant (J s)
k = 1/(4*pi*eps0) # Coulomb's constant

J_per_eV = e
J_per_keV = 1.0e3 * e
J_per_MeV = 1.0e6 * e


@dataclass
class Species:
    name: str
    Z: float          # charge state, electron should use Z = -1
    mass: float       # kg
    density: float    # m^-3
    T_eV: float       # eV


@dataclass
class Projectile:
    name: str
    Z: float          # charge state
    mass: float       # kg
    E0_MeV: float     # initial kinetic energy in MeV


def temperature_to_joule(T_eV):
    return T_eV * e


def velocity_from_energy(E_J, mass):
    """
    Nonrelativistic kinetic energy: E = 1/2 m v^2.
    Fine for MeV ions like alpha particles in ICF-ish plasmas.
    """
    return np.sqrt(2.0 * E_J / mass)


def debye_length(species_list):
    """
    Multi-species Debye length:
        1/lambda_D^2 = sum_s n_s Z_s^2 e^2 / (eps0 k_B T_s)

    Since T is passed in eV, k_B T = T_eV * e.
    """
    inv_lam2 = 0.0

    for s in species_list:
        T_J = temperature_to_joule(s.T_eV)
        inv_lam2 += s.density * (s.Z ** 2) * e**2 / (eps0 * T_J)

    return 1.0 / np.sqrt(inv_lam2)


def thermal_speed(species):
    """
    Uses v_th = sqrt(kT/m).
    Since T is in eV, kT = T_eV * e.
    """
    return np.sqrt(temperature_to_joule(species.T_eV) / species.mass)


def coulomb_log_regularized(projectile, background, E_J, all_species):
    """
    Regularized Coulomb logarithm.

    b_max ~ Debye length.
    b_min ~ max(classical 90-degree impact parameter, de Broglie length).

    This is a simple classical/semiclassical choice.
    For serious Li-Petrasso work, this is the part we may replace.
    """
    v = velocity_from_energy(E_J, projectile.mass)
    vb = thermal_speed(background)

    # crude effective relative speed
    vrel = np.sqrt(v**2 + vb**2)

    mu = projectile.mass * background.mass / (projectile.mass + background.mass)

    # Classical distance of closest approach / 90-degree impact scale
    b90 = abs(projectile.Z * background.Z) * k * e**2 / (mu * vrel**2)

    # Quantum diffraction scale
    bq = h / (2.0 * mu * vrel)

    bmin = np.maximum(b90, bq)
    bmax = debye_length(all_species)

    # Regularized form avoids negative or pathological ln Lambda
    lnL = 0.5 * np.log(1.0 + (bmax / bmin) ** 2)

    return lnL


def g_maxwellian(x):
    """
    Maxwellian velocity function:
        G(x) = erf(x) - 2x/sqrt(pi) exp(-x^2)

    This suppresses stopping by species whose thermal speed is much larger
    than the projectile speed.
    """
    return erf(x) - (2.0 * x / np.sqrt(np.pi)) * np.exp(-x**2)


def stopping_power_species(projectile, background, E_J, all_species):
    """
    Approximate classical plasma stopping power contribution from one species.

    Returns positive S(E) = -dE/dx in J/m.

    Model:
        S_b ~ 4 pi n_b (Z_p Z_b e^2 / 4pi eps0)^2 / (m_b v^2)
              * ln Lambda_b * G(v / sqrt(2) v_th,b)

    This is a useful weak-coupling Maxwellian model, not full Li-Petrasso.
    """
    v = velocity_from_energy(E_J, projectile.mass)
    vb = thermal_speed(background)

    x = v / (np.sqrt(2.0) * vb)
    G = g_maxwellian(x)

    lnL = coulomb_log_regularized(projectile, background, E_J, all_species)

    prefactor = (
        4.0
        * np.pi
        * background.density
        * (projectile.Z * background.Z * k * e**2) ** 2
        / (background.mass * v**2)
    )

    S = prefactor * lnL * G

    return np.maximum(S, 0.0)


def stopping_power_total(projectile, species_list, E_J):
    """
    Total positive stopping power:
        S(E) = -dE/dx = sum_b S_b(E)

    Returns J/m.
    """
    S = np.zeros_like(E_J, dtype=float)

    for b in species_list:
        S += stopping_power_species(projectile, b, E_J, species_list)

    return S


def stopping_length(projectile, species_list, Emin_keV=0.0, ngrid=5000):
    """
    Integrates:
        L = integral_{Emin}^{E0} dE / S(E)

    where S(E) = -dE/dx.

    Returns:
        E_MeV grid
        x_m cumulative range from E0 down to E
        S_J_per_m
    """
    E0_J = projectile.E0_MeV * J_per_MeV
    Emin_J = Emin_keV * J_per_keV

    # Use log grid because stopping changes strongly at low energy
    E_desc = np.geomspace(E0_J, Emin_J, ngrid)
    E_asc = E_desc[::-1]

    S_asc = stopping_power_total(projectile, species_list, E_asc)

    # Avoid divide-by-zero crash
    tiny = 1.0e-300
    invS = 1.0 / np.maximum(S_asc, tiny)

    # x(E): range accumulated from Emin upward
    x_from_Emin = cumulative_trapezoid(invS, E_asc, initial=0.0)

    # Total stopping length from E0 to Emin
    L_total = x_from_Emin[-1]

    # Convert to cumulative distance from initial E0 downward
    x_from_E0_desc = L_total - x_from_Emin[::-1]

    E_MeV_desc = E_desc / J_per_MeV
    S_desc = S_asc[::-1]

    return E_MeV_desc, x_from_E0_desc, S_desc, L_total


if __name__ == "__main__":

    # ------------------------------------------------------
    # Example: 3.5 MeV alpha in equimolar DT plasma
    # ------------------------------------------------------
    # Pick plasma conditions here.
    # Number densities are in m^-3.
    # Temperatures are in eV.
    #
    # Example:
    #   ne = 1e31 m^-3
    #   Te = Ti = 5 keV
    #
    # For equimolar D-T with full ionization:
    #   nD = nT = ne / 2
    # ------------------------------------------------------

    ne = 2.0e25
    Te_eV = 50.0
    Ti_eV = 1.0

    species = [
        Species(name="e", Z=-1.0, mass=me, density=ne, T_eV=Te_eV),
        Species(name="H", Z=1.0, mass=2.014 * mp, density= 0.8 * ne/1.2, T_eV=Ti_eV),
        Species(name="He", Z=2.0, mass=4.002603 * mp, density=0.2 * ne/1.2, T_eV=Ti_eV),
    ]

    Si_KeV = Projectile(
        name="Si",
        Z=14.0,
        mass= 28*mp,
        E0_MeV=0.001
    )

    Si_10KeV = Projectile(
        name="Si",
        Z=14.0,
        mass= 28*mp,
        E0_MeV=0.01
    )

    Si_100KeV = Projectile(
        name="Si",
        Z=14.0,
        mass= 28*mp,
        E0_MeV=0.1
    )

    Si_MeV = Projectile(
        name="Si",
        Z=14.0,
        mass= 28*mp,
        E0_MeV=1
    )

    E_KeV, x_m1, S_J_per_m1, L_m1 = stopping_length(
        Si_KeV,
        species,
        Emin_keV=0.001,
        ngrid=8000,
    )

    E_10KeV, x_m2, S_J_per_m2, L_m2 = stopping_length(
        Si_10KeV,
        species,
        Emin_keV=0.001,
        ngrid=8000,
    )

    E_100KeV, x_m3, S_J_per_m3, L_m3 = stopping_length(
        Si_100KeV,
        species,
        Emin_keV=0.001,
        ngrid=8000,
    )

    E_MeV, x_m, S_J_per_m, L_m = stopping_length(
        Si_MeV,
        species,
        Emin_keV=0.001,
        ngrid=8000,
    )

    
    # print("Projectile:", Si.name)
    # print(f"Initial energy: {Si.E0_MeV:.3f} MeV")
    # # print(f"Stopping length to 1 keV: {L_m:.6e} m")
    # print(f"Stopping length: {L_m * 1e3} mm")

    fig,ax = plt.subplots()

    ax.plot(x_m1*1e3, S_J_per_m1 / (10**6 * e), label="1 KeV Si")
    ax.plot(x_m2*1e3, S_J_per_m2 / (10**6 * e), label="10 KeV Si")
    ax.plot(x_m3*1e3, S_J_per_m3 / (10**6 * e), label="100 KeV Si")
    #ax.plot(x_m*1e3, S_J_per_m / (1e6 * e), label="1 MeV Si")
    
    ax.set_xlabel("Distance (mm)")
    ax.set_ylabel("Stopping power (keV/mm)")
    ax.set_title("Stopping power vs distance for Si projectiles in H(0.8) + He(0.2) plasma")
    plt.legend()
    plt.show()

    # Optional: save table
    # out = np.column_stack([
    #     E_MeV,
    #     x_m,
    #     x_m * 1e6,
    #     S_J_per_m,
    #     S_J_per_m / J_per_MeV,  # MeV/m
    #     S_J_per_m / J_per_MeV / 1e6,  # MeV/micron
    # ])

    # header = (
    #     "E_MeV x_m x_micron S_J_per_m "
    #     "S_MeV_per_m S_MeV_per_micron"
    # )

    # np.savetxt("stopping_length_output.txt", out, header=header)

    # print("Saved: stopping_length_output.txt")