import os
import pandas as pd
import numpy as np


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
            "fitfunc": fitfunc, "maxcov": numpy.max(pcov), "rawres": (guess, popt, pcov)}

def load_csv(path, parm=False):
    """load a csv and return dataframe (optional) with parm key added"""
    _ = pd.read_csv(path)
    if parm:
        # add the parms to fit to
        df_energy = pd.read_csv(parm)
        _ = _.merge(df_energy, on="i")

    return _


def fit_charge(chg_n, df, key="parm"):
    fit_x = fit_poly(df[key], df[f"x_c{chg_n}"])
    fit_y = fit_poly(df[key], df[f"y_c{chg_n}"])
    fit_z = fit_poly(df[key], df[f"z_c{chg_n}"])
    return [fit_x, fit_y, fit_z]


def fit_fmdcm():
    """fit the entire fMDCM model"""
    pass


def update_chg_models():
    """update existing charge models from a previous fit to match the fitted function, useful for recalculating the
    errors without any approximation """
    pass




