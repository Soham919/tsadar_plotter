import numpy as np
from scipy.optimize import brentq

def compression_mach_error(n1, dn1, n2, dn2, gamma=5/3):
    """Density compression and ideal-shock Mach number with 1-sigma propagated errors."""
    r = n2 / n1
    dr = r * np.sqrt((dn1/n1)**2 + (dn2/n2)**2)

    if not (1 < r < (gamma + 1)/(gamma - 1)):
        raise ValueError("Compression ratio outside the physical shock range.")

    M = np.sqrt(2*r / ((gamma + 1) - (gamma - 1)*r))
    dM = M * (gamma + 1) * dr / (r * ((gamma + 1) - (gamma - 1)*r))

    return r, dr, M, dM



def mach_from_temperature(T1, dT1, T2, dT2, gamma=5/3):
    """Mach number from ideal-shock temperature jump, with propagated 1σ error."""
    R = T2/T1
    dR = R*np.hypot(dT1/T1, dT2/T2)

    def ratio(M):
        return ((2*gamma*M**2-(gamma-1))*((gamma-1)*M**2+2))/((gamma+1)**2*M**2)

    if R <= 1:
        raise ValueError("Temperature ratio must exceed 1.")

    M = brentq(lambda m: ratio(m)-R, 1, 1e4)
    h = 1e-4*M
    dR_dM = (ratio(M+h)-ratio(M-h))/(2*h)
    dM = dR/abs(dR_dM)

    return M, dM


def mach_from_velocity(u1, du1, u2, du2, us, dus, gamma=5/3):
    """Mach number from laboratory-frame velocity jump, with propagated 1σ error."""
    a, b = us-u1, us-u2
    r = a/b
    rmax = (gamma+1)/(gamma-1)

    if not (1 < r < rmax):
        raise ValueError("Velocity-derived compression is outside ideal-shock range.")

    dr = np.sqrt((du1/b)**2 + (a*du2/b**2)**2 +
                 ((u1-u2)*dus/b**2)**2)

    M = np.sqrt(2*r/((gamma+1)-(gamma-1)*r))
    dM = M*(gamma+1)*dr/(r*((gamma+1)-(gamma-1)*r))

    return M, dM


# Example values only — replace with your measurements
M_T, dM_T = mach_from_temperature(
    T1=205, dT1=9.1,
    T2=480, dT2=20.2
)

M_v, dM_v = mach_from_velocity(
    u1=0, du1=1,
    u2=29, du2=2,
    us=45, dus=3
)

print(f"Temperature Mach = {M_T:.2f} ± {dM_T:.2f}")
print(f"Velocity Mach    = {M_v:.2f} ± {dM_v:.2f}")


# Replace with your measured values
r, dr, M, dM = compression_mach_error(
    n1=0.171, dn1=0.01,
    n2=0.379, dn2=0.011
)

print(f"Compression ratio = {r:.2f} ± {dr:.2f}")
print(f"Mach number       = {M:.2f} ± {dM:.2f}")
