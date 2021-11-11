import pickle

import numpy as np
import pandas as pd
from ase import Atoms
from dscribe.descriptors import SOAP
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split


def min_max_scale(data):
    #  Scale data so it lies between 0,1
    data = data.T
    print(data.shape)
    _mins = []
    _maxs = []
    _scaled_data = []
    for i, variable in enumerate(data):
        _mins.append(min(variable))
        _maxs.append(max(variable))
        data[i] = [(x - min(variable)) / (max(variable) - min(variable)) for x in variable]
    return data.T, _mins, _maxs


def z_scale(data):
    #  Scale data so it lies between 0,1
    data = data.T
    print(data.shape)
    _means = []
    _stds = []
    _scaled_data = []
    for i, variable in enumerate(data):
        _means.append(np.mean(variable))
        _stds.append(np.std(variable))
        data[i] = [(x - np.mean(variable)) / np.std(variable) for x in variable]
    return data.T, _means, _stds


df = pd.read_pickle("../MD/25.pkl")

ONLY_HEAVY = False
MODEL = "KR"
TEST_SIZE = 0.5
CROSS_FOLDS = 2
RANDOM_SEED = 0
RCUT = 5
NMAX = 4
LMAX = 5
SCALING = "None"

#  SOAP descriptors
atom_types = ["C", "C", "H", "H", "H", "O", "O", "C", "H", "H", "H"]
atom_types_no_H = ["C", "C", "O", "O", "C"]
heavy_atoms = [1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0]
soap_desc = SOAP(species=["C", "H", "O"], rcut=RCUT, nmax=NMAX, lmax=LMAX, crossover=True)
#  Give charge indices for fitting, 1th indexed
# charge_indices = [1, 2, 5, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 22, 24]
charge_indices = range(1, 25)

indices = []
for i in range(len(df)):
    indices.append(heavy_atoms)

samples = []
for i in range(len(df)):
    t = i // 500
    if t > 3:
        t = 3
    string = list(df.iloc[i]["Z"])
    if ONLY_HEAVY:
        string = [x for i, x in enumerate(string) if heavy_atoms[i] == 1]
    if ONLY_HEAVY:
        atoms = Atoms(atom_types_no_H, string)
    else:
        atoms = Atoms(atom_types, string)
    samples.append(atoms)

der, des = soap_desc.derivatives(samples, method="analytical", return_descriptor=True, n_jobs=8)
# des = soap_desc.create(samples, n_jobs=2)
with open('des.pkl', 'wb') as f:
    pickle.dump(des, f)
with open('der.pkl', 'wb') as f:
    pickle.dump(der, f)

# print(der)
# print(der.shape)

best_validation = -np.inf
try:
    with open("best_score.txt") as f:
        best_validation = float(f.read())
        print("Best validation score: ", best_validation)
except FileNotFoundError:
    print("Best validation score NOT FOUND")
    pass

iterations = 1000
for i in range(iterations):

    keys = [[f"x_c{i - 1}", f"y_c{i - 1}", f"z_c{i - 1}"] for i in charge_indices]
    __ = []
    for _ in range(len(charge_indices)):
        for e in range(3):
            __.append(keys[_][e])
    y = df[__].values

    #  Reshape the training data
    a_, b_, c_ = des.shape
    if ONLY_HEAVY:
        X = des.reshape((len(df), c_ * sum(heavy_atoms)))
    else:
        X = des.reshape((len(df), c_ * len(heavy_atoms)))

    if SCALING == "minmax":
        y, mins, maxs = min_max_scale(y)
    elif SCALING == "z":
        y, means, stds = z_scale(y)
    else:
        pass

    #  Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)#, random_state=RANDOM_SEED)
    _X_train, _y_train = X_train, y_train

    for n in [1]:
        print(n)
        charge_n = []
        train_scores = []
        val_score = []
        val_data = []
        models = []

        #X_train, y_train = _X_train[:len(_X_train) // n], _y_train[:len(_X_train) // n]
        print("Size of X_train ", len(X_train))
        #  Choose model
        if MODEL == "GPR":
            kernel = DotProduct() + WhiteKernel()
            model = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(X_train, y_train)
        elif MODEL == "LR":
            model = LinearRegression().fit(X_train, y_train)
        elif MODEL == "RCV":
            model = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1]).fit(X_train, y_train)
        elif MODEL == "KR":
            model = KernelRidge(alpha=1.0)
            model.fit(X_train, y_train)
        else:
            print("Model not available")
            raise Exception

        #  Get scores
        train_score = model.score(X_train, y_train)
        print("training score: ", train_score)
        test_score = model.score(X_test, y_test)
        print("test set score: ", test_score, (len(y_test)))
        val_scores = cross_val_score(model, X_test, y_test, cv=CROSS_FOLDS)
        mean_val_score = np.mean(val_scores)
        print(f"Mean validation score, of {CROSS_FOLDS} fold: ", mean_val_score)

    models.append(model)

    charge_n.append(i)
    train_scores.append(train_score)
    val_score.append(mean_val_score)
    val_data.append(val_scores)

    mean_val_score = test_score

    if mean_val_score > best_validation:
        print("Saving best model with validation accuracy: ", mean_val_score)
        with open('best_model.pkl', 'wb') as f:
            pickle.dump(model, f)
        with open("best_score.txt", "w") as f:
            f.write("{}".format(mean_val_score))
        best_validation = mean_val_score

    with open('models.pkl', 'wb') as f:
        pickle.dump(models, f)
