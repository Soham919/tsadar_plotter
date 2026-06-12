import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def fit_1d_curve(
    x,
    y,
    fit_func,
    p0=None,
    bounds=(-np.inf, np.inf),
    xfit=None,
    sigma=None,
    absolute_sigma=False,
    plot=True,
    xlabel="x",
    ylabel="y",
    title="Curve Fit",
    data_label="Data",
    fit_label="Fit",
):
    """
    Simple 1D curve fitter.

    Parameters
    ----------
    x, y : array-like
        Data arrays.
    
    fit_func : callable
        Function to fit. Must have form:
            fit_func(x, p1, p2, p3, ...)
    
    p0 : list or array, optional
        Initial guesses for fit parameters.
        Example: p0=[1, 0, 2]
    
    bounds : tuple, optional
        Lower and upper bounds for parameters.
        Example: bounds=([0, -np.inf], [np.inf, np.inf])
    
    xfit : array-like, optional
        X values where best-fit curve is evaluated.
        If None, a smooth range is made automatically.
    
    sigma : array-like, optional
        Uncertainty in y values.
    
    absolute_sigma : bool
        If True, sigma is treated as absolute uncertainty.
    
    plot : bool
        If True, makes a plot.

    Returns
    -------
    result : dict
        Dictionary containing:
            popt      : best-fit parameters
            perr      : 1-sigma parameter uncertainties
            pcov      : covariance matrix
            xfit      : x values for fitted curve
            yfit      : fitted y values
            residuals : y - fit_func(x, *popt)
    """

    x = np.asarray(x)
    y = np.asarray(y)

    # Remove NaNs / infs
    good = np.isfinite(x) & np.isfinite(y)
    if sigma is not None:
        sigma = np.asarray(sigma)
        good = good & np.isfinite(sigma)

    x_clean = x[good]
    y_clean = y[good]

    if sigma is not None:
        sigma_clean = sigma[good]
    else:
        sigma_clean = None

    # Fit
    popt, pcov = curve_fit(
        fit_func,
        x_clean,
        y_clean,
        p0=p0,
        bounds=bounds,
        sigma=sigma_clean,
        absolute_sigma=absolute_sigma,
        maxfev=100000,
    )

    # Parameter errors
    perr = np.sqrt(np.diag(pcov))

    # Smooth fitted curve
    if xfit is None:
        xfit = np.linspace(np.min(x_clean), np.max(x_clean), 1000)
    else:
        xfit = np.asarray(xfit)

    yfit = fit_func(xfit, *popt)

    # Residuals on original data points
    residuals = y_clean - fit_func(x_clean, *popt)

    # Print fit result
    print("\nBest-fit parameters:")
    for i, (val, err) in enumerate(zip(popt, perr)):
        print(f"p{i} = {val:.6g} ± {err:.3g}")

    # Plot
    if plot:
        fig, ax = plt.subplots(figsize=(7, 5))

        ax.plot(x_clean, y_clean, "o", label=data_label)
        ax.plot(xfit, yfit, "-", linewidth=2, label=fit_label)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    result = {
        "popt": popt,
        "perr": perr,
        "pcov": pcov,
        "xfit": xfit,
        "yfit": yfit,
        "residuals": residuals,
        "x_clean": x_clean,
        "y_clean": y_clean,
    }

    return result