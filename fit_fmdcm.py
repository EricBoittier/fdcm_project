import os
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt

def angle(a, b, c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    a3 = np.degrees(np.arccos(cosine_angle))
    return a3


def fit_poly(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq",
    "period" and "fitfunc" '''
    tt = np.array(tt)  # x
    yy = np.array(yy)  # y = local position (x/y/z)
    def polyfunc(x, a, b, c, d): return a * x + b * x ** 2 + c * x ** 3 + d
    guess = [1, 1, 1, 1]
    popt, pcov = scipy.optimize.curve_fit(polyfunc, tt,
                                          yy, p0=guess, maxfev=10000)
    a, b, c, d = popt
    fitfunc = lambda x: a * x + b * x ** 2 + c * x ** 3 + d
    return {"a": a, "b": b, "c": c, "d": d,
            "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess, popt, pcov)}

def load_csv(path, parm=False):
    """load a csv and return dataframe (optional) with parm key added"""
    _ = pd.read_csv(path)
    if parm:
        # add the parms to fit to
        df_energy = pd.read_csv(parm)
        _ = _.merge(df_energy, on="i")

    return _


def fit_charge(chg_n, df, key="parm", plot=False):
    fit_x = fit_poly(df[key], df[f"x_c{chg_n}"])
    fit_y = fit_poly(df[key], df[f"y_c{chg_n}"])
    fit_z = fit_poly(df[key], df[f"z_c{chg_n}"])
    fits_ = [fit_x, fit_y, fit_z]

    if plot:
        """Make a plot showing the fitted function"""
        fig, ax = plt.subplots(3, 1, sharex=True)
        ax[0].scatter(df[key], df[f"x_c{chg_n}"], s=10.5, c="red")
        ax[1].scatter(df[key], df[f"y_c{chg_n}"], s=10.5, c="green")
        ax[2].scatter(df[key], df[f"z_c{chg_n}"], s=10.5, c="blue")

        ax[0].plot(df[key], fit_x["fitfunc"](df[key]), "--", c="k", linewidth=0.675)
        ax[1].plot(df[key], fit_y["fitfunc"](df[key]), "--", c="k", linewidth=0.675)
        ax[2].plot(df[key], fit_z["fitfunc"](df[key]), "--", c="k", linewidth=0.675)

        ax[0].set_ylabel("$e_x$")
        ax[1].set_ylabel("$e_y$")
        ax[2].set_ylabel("$e_z$")

        for i in range(3):
            plt.text(0.5, 0.5, f"$q_{{{chg_n + 1}}}$",
                     rotation=0, verticalalignment="center",
                     transform=ax[i].transAxes)

            covariance = fits_[i]["maxcov"]
            plt.text(1, 1, "$\sigma_{X,Y}$ =" + "{:.3e}".format(covariance),
                     rotation=0, verticalalignment="center",
                     transform=ax[i].transAxes)
        plt.savefig("{0}_c{1}.pdf".
                    format(plot, chg_n + 1))

    return fits_


def fit_fmdcm():
    """fit the entire fMDCM model"""
    pass


def update_chg_models():
    """update existing charge models from a previous fit to match the fitted function, useful for recalculating the
    errors without any approximation """
    pass




