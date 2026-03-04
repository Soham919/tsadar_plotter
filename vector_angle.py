import numpy as np

def angle_between(v1, v2, degrees=False):
    """
    Compute angle between two vectors in cartesian geometry.

    Parameters
    ----------
    v1, v2 : array-like
        Input vectors (1D arrays)
    degrees : bool
        If True, return angle in degrees. Otherwise radians.

    Returns
    -------
    angle : float
    """
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)

    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)

    if norm1 == 0 or norm2 == 0:
        raise ValueError("Cannot compute angle with zero-length vector")

    cos_theta = np.dot(v1, v2) / (norm1 * norm2)

    # Numerical safety: clamp to [-1, 1]
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    angle = np.arccos(cos_theta)

    if degrees:
        angle = np.degrees(angle)

    return angle

def cross_product(a, b):
    """
    Returns the cross product of two 3D vectors a and b.
    Each vector should be an iterable of length 3.
    """
    if len(a) != 3 or len(b) != 3:
        raise ValueError("Both vectors must have exactly 3 components")

    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ]


## Working in Si3N4 coordinate system

v1 = [-2.5, -2.5, 0]
v2 = [-2.5, 2.5, 0]
n = cross_product(v1,v2)
n = n/np.linalg(n)
a = 
b = 
stalk = a-b

print(f"normal vector of the Si3N4 frame = {n}")
print(f"Angle between vectors is = {angle_between(v1,v2,True)}")