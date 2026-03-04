import numpy as np
import matplotlib.pyplot as plt

def supergaussian_fwhm(x, center=0.0, fwhm=1.0, m=2,
                       amplitude=1.0, baseline=0.0):
    """
    Super-Gaussian profile parameterized by FWHM.

    Parameters
    ----------
    x : array_like
        Points at which to evaluate.
    center : float
        Center position (x0).
    fwhm : float
        Full width at half maximum.
    m : float
        Order / flatness parameter.
        m=1 -> Gaussian
        larger m -> flatter top, steeper edges.
    amplitude : float
        Peak height above baseline.
    baseline : float
        Vertical offset.

    Returns
    -------
    y : ndarray
        Super-Gaussian evaluated at x.
    """

    x = np.asarray(x, dtype=float)

    if fwhm <= 0:
        raise ValueError("fwhm must be > 0")
    if m <= 0:
        raise ValueError("m must be > 0")

    # Convert FWHM to width parameter
    w = fwhm / (2 * (np.log(2))**(1/(2*m)))

    return baseline + amplitude * np.exp(-np.abs((x - center)/w)**(2*m))

t = np.linspace(-1,8,500)
a1 = 300/600 # db power in TW
a2 = 5/100 # 100 ps probe power in TW
a3 = 45/3700 # 3.7 ns probe power in TW
db = supergaussian_fwhm(t,0.45,0.6,3,a1)  # Drive beams
pb = supergaussian_fwhm(t,6,0.1,1,a2)  # Probe beam 100 ps
pb2 = supergaussian_fwhm(t,6,3.7,10,a3)  # Probe beam 100 ps



fig, ax = plt.subplots(figsize=(6,3))
ax.plot(t, db, color ='red',linewidth=2,label=r'Drive Beam : 3$\omega$')
ax.plot(t, pb, color ='limegreen',linewidth=2, label = r'Probe Beam : 2$\omega$')
ax.set_title(r'Timing diagram for RIDs : 103409, 103615, 103653', pad=15)
ax.set_xlabel(r"$t (ns)$")
ax.set_ylabel(r"$Approx. P (TW)$")
ax.set_ylim([0, 0.6])
ax.set_xlim([-0.2, 8])

t1 = ax.axvline(x=4, color='limegreen', linestyle='--')
t2 = ax.axvline(x=8, color='limegreen', linestyle='--')


ax.axvspan(4,8, ymin = 0, ymax = 1,
                    color='greenyellow',        # fill color
                    alpha=0.2,           # transparency
                    edgecolor='limegreen',    # outline color
                    linewidth=1.5,       # outline width
                    linestyle='-',      # outline style (optional)
                    label='Timing variable')

ax.grid(True)
plt.legend()
plt.show()

