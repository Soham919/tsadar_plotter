import numpy as np

def angle_between(v1, v2, degrees=False):
    """
    Compute angle between two vectors.

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