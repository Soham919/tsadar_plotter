
import numpy as np

def spherical_to_cartesian(r, theta, phi):
    """Theta measured from +z, phi from +x. Angles in degrees."""
    theta, phi = np.deg2rad([theta, phi])
    return r * np.array([np.sin(theta)*np.cos(phi),
                         np.sin(theta)*np.sin(phi),
                         np.cos(theta)])

def unit(v):
    v = np.asarray(v, dtype=float)
    norm = np.linalg.norm(v)
    if norm == 0:
        raise ValueError("Cannot normalize a zero-length vector.")
    return v / norm

def rotation_between_vectors(a, b):
    """Minimum-angle rotation matrix mapping vector a onto vector b."""
    a, b = unit(a), unit(b)
    dot = np.clip(np.dot(a, b), -1.0, 1.0)

    if np.isclose(dot, 1.0):
        return np.eye(3)

    if np.isclose(dot, -1.0):
        basis = np.eye(3)[np.argmin(np.abs(a))]
        axis = unit(np.cross(a, basis))
        return 2*np.outer(axis, axis) - np.eye(3)

    v = np.cross(a, b)
    K = np.array([[0, -v[2], v[1]],
                  [v[2], 0, -v[0]],
                  [-v[1], v[0], 0]])

    return np.eye(3) + K + K@K / (1.0 + dot)

def rotate_about_y(point, angle_deg, center):
    """Rotate point around a y-parallel axis passing through center."""
    phi = np.deg2rad(angle_deg)
    c, s = np.cos(phi), np.sin(phi)

    R_y = np.array([[ c, 0,  s],
                    [ 0, 1,  0],
                    [-s, 0,  c]])

    return center + R_y @ (point - center)

def transform_laser(port_sph, focus_sph, tracer_exp,
                    tracer_flash, axis_exp, axis_flash,
                    exp_to_cm=1.0):

    port_exp = spherical_to_cartesian(*port_sph) * exp_to_cm
    focus_exp = spherical_to_cartesian(*focus_sph) * exp_to_cm
    tracer_exp = np.asarray(tracer_exp, dtype=float) * exp_to_cm
    tracer_flash = np.asarray(tracer_flash, dtype=float)

    R = rotation_between_vectors(axis_exp, axis_flash)

    port_flash = tracer_flash + R @ (port_exp - tracer_exp)
    focus_flash = tracer_flash + R @ (focus_exp - tracer_exp)

    k_exp = unit(focus_exp - port_exp)
    k_flash = unit(focus_flash - port_flash)

    angle_exp = np.degrees(np.arccos(
        np.clip(abs(k_exp @ unit(axis_exp)), 0.0, 1.0)))

    angle_flash = np.degrees(np.arccos(
        np.clip(abs(k_flash @ unit(axis_flash)), 0.0, 1.0)))

    return port_flash, focus_flash, angle_exp, angle_flash, R


# ======================= INPUTS =======================

port_sph = (333.0, 74.0, 197.0)
focus_sph = (0.041, 113.558, 172.582)

axis_point1 = np.array([-1.18886, 0.0426033, -0.431458])/10
axis_point2 = np.array([-0.384439, 0.0426055, -1.02552])/10
axis_exp = axis_point1 - axis_point2

tracer_exp = np.array([-0.0786648, 0.00426033, -0.0728491])
tracer_flash = np.array([0.0, 0.0510, 0.0])

axis_flash = np.array([0.0, -1.0, 0.0])
exp_to_cm = 1.0

# Extra azimuthal rotation:
# Set to None to automatically align the laser projection with +x.
# Or enter an angle in degrees, e.g. 30.0, -45.0, 90.0.
azimuth_rotation = None

# ===================== CALCULATION ===================

lens, target, angle_exp, angle_flash, R = transform_laser(
    port_sph, focus_sph, tracer_exp, tracer_flash,
    axis_exp, axis_flash, exp_to_cm
)

# Determine extra azimuthal rotation.
# Align laser focus relative to Cu tracer along the x-axis.
d = target - tracer_flash

if azimuth_rotation is None:
    phi_deg = np.degrees(np.arctan2(d[2], d[0]))
else:
    phi_deg = azimuth_rotation
    
# Rotate both positions about the Cu tracer center.
lens = rotate_about_y(lens, phi_deg, tracer_flash)
target = rotate_about_y(target, phi_deg, tracer_flash)

# ======================== OUTPUT ======================

print(f"Incidence angle (experiment): {angle_exp:.4f} deg")
print(f"Incidence angle (FLASH):      {angle_flash:.4f} deg")
print(f"Extra azimuthal rotation:    {phi_deg:.4f} deg")

print("\nRotation matrix (experimental axis alignment):")
print(R)

print("\nFLASH laser parameters (cm):")
for name, value in zip("XYZ", lens):
    print(f"ed_lens{name}_1 = {value:.10e}")

for name, value in zip("XYZ", target):
    print(f"ed_target{name}_1 = {value:.10e}")

print("\nLens-to-target distance (cm):",
      np.linalg.norm(target - lens))

print("Focus relative to tracer (um):",
      (target - tracer_flash) * 1e4)
