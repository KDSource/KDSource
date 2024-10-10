#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for bandwidth selection methods

Eventually, all the content of this module should be merged into KDEpy.
"""

from KDEpy import TreeKDE

from joblib import Parallel, delayed

import matplotlib.pyplot as plt

import numpy as np

from sklearn.model_selection import KFold
from sklearn.neighbors import NearestNeighbors


def bw_silv(dim, N_eff):
    """
    Silverman's Rule.

    Estimates optimal bandwidth based on dimension and effective number
    of samples.

    Parameters
    ----------
    dim: int
        Dimension of the KDE model.
    N_eff: float
        Effective number of samples, defined as:
            N_eff = sum(weights)**2 / sum(weights**2)

    Returns
    -------
    bw_opt: float
        Bandwidth resulting from Silverman's Rule.
    """
    return (4 / (2 + dim)) ** (1 / (4 + dim)) * N_eff ** (-1 / (4 + dim))


def bw_knn(data, weights=None, K_eff=100, k=None, batch_size=10000):
    """
    K Nearest Neighbors method.

    For each sample, computes its bandwidth as the distance to the k-th
    neighbor, within each batch.

    Parameters
    ----------
    data: array-like
        Array of samples. Must have shape (obs, dim).
    weights: array-like, optional
        Array of sample weights. By default all weights are 1.
    K_eff: float, optional
        Effective k for all dataset. It is used to compute the k so that
        the estimated number of neighbors in all dataset is K_eff.
    k: int, optional
        Number of neighbors within each batch. If set, overrides the
        value estimated with K_eff.
    batch_size: int, optional
        Batch size for KNN search.

    Returns
    -------
    bw_opt: numpy.array
        Array of KNN bandwidths.
    """
    N, dim = data.shape
    batches = int(max(1, np.round(N / batch_size)))
    N_eff = (
        N
        if weights is not None
        else np.sum(weights) ** 2 / np.sum(weights ** 2)
    )
    if k is None:  # Compute k based on K_eff
        k_float = batch_size * K_eff / N_eff
        k = int(max(1, round(k_float)))
        f_k = k_float / k  # Correction factor for k
    else:  # k set by user
        f_k = 1.0
        K_eff = k * N_eff / batch_size
    print(
        "Using k = {} neighbors per batch (batch_size = {})".format(
            k, batch_size
        )
    )
    print("Correction factor: f_k = k_float / k = {}".format(f_k))
    print("Effective total neighbors: K_eff = {}".format(K_eff))
    bw_knn = np.empty(len(data))
    idx = 0
    for batch, batch_data in enumerate(np.array_split(data, batches)):
        print("batch =", batch + 1, "/", batches)
        knn = NearestNeighbors(
            n_neighbors=k + 1, n_jobs=-1
        )  # Add 1 to count self point
        knn.fit(batch_data)
        ds, idxs = knn.kneighbors(batch_data)
        bws = ds[:, -1] * (f_k) ** (
            1 / dim
        )  # Select k-th column and apply correction factor
        bw_knn[idx : idx + len(batch_data)] = bws
        idx += len(batch_data)
    return bw_knn


def _kde_cv_score(bw, data, weights=None, n_splits=10, modelKDE=TreeKDE):
    """
    Computes cross-validation likelihood score.

    Uses a K-Fold scheme to compute cross validation likelihood score,
    defined as the mean log-density of KDE model constructed with train
    samples, evaluated on test samples.

    Parameters
    ----------
    bw: float or array-like
        Bandwidth of KDE model. Can be constant or adaptive (one for
        each sample).
    data: array-like
        Array of samples. Must have shape (obs, dim).
    weights: array-like, optional
        Array of sample weights. By default all weights are 1.
    n_splits: int, optional
        Number of folds. Must be at least 2.
    modelKDE: NaiveKDE or TreeKDE, optional
        KDE model to be used for evaluating scores.
    """
    if np.ndim(bw) == 1:
        if len(bw) < len(data):
            raise ValueError("If bw is array, must have same len as data")
        if len(bw) > len(data):
            print("Warning: bw longer than data.")

    folds = KFold(n_splits=n_splits)
    scores = []
    for train_idx, test_idx in folds.split(data):
        data_train, data_test = data[train_idx], data[test_idx]
        if weights is not None:
            weights_train, weights_test = weights[train_idx], weights[test_idx]
        else:
            weights_train = weights_test = None
        if np.ndim(bw) == 1:
            bw_train = bw[train_idx]
        else:
            bw_train = bw
        kde = modelKDE(bw=bw_train)
        kde.fit(data_train, weights=weights_train)
        if weights_test is not None:
            score = np.mean(weights_test * np.log(kde.evaluate(data_test)))
        else:
            score = np.mean(np.log(kde.evaluate(data_test)))
        scores.append(score)
    return np.mean(scores)


def bw_mlcv(data, weights=None, n_splits=10, seed=None, grid=None, show=True):
    """
    Maximum Likelihood Cross-Validation bandwidth.

    Builds a grid of bandwidths, and computes the cross-validation
    likelihood score for each. Chooses the bandwidth that maximizes the
    score, or raises exception if maximum score is found in beginning or
    end of bandwidth grid. Also plots the scores over the bandwidth
    grid.

    Bandwidths grid is computed with seed and grid, so that:
        bw_grid[i] = seed * grid[i]

    Parameters
    ----------
    data: array-like
        Array of samples. Must have shape (obs, dim).
    weights: array-like, optional
        Array of sample weights. By default all weights are 1.
    n_splits: int, optional
        Number of folds for cross-validation. Must be at least 2.
    seed: float or array-like
        Seed bandwidth for computing bandwidth grid. Can be constant or
        adaptive. It is recommended to use the output of other bandwidth
        selection method. By default it is computed with Silverman's
        Rule.
    grid: array-like
        Grid of factors for computing bandwidth grid. Default:
        numpy.logspace(-1, 1, 10)
    show: bool
        Whether to show the scores plot.

    Returns
    -------
    bw_opt: float or numpy.array
        Optimal bandwidth, with same shape as seed.
    """
    if seed is None:
        N_eff = (
            len(data)
            if weights is None
            else np.sum(weights) ** 2 / np.sum(weights ** 2)
        )
        seed = bw_silv(np.shape(data)[1], N_eff)
    if grid is None:
        print("No grid specified. Using np.logspace(-0.1, 0.1, 10).")
        print("If fitting fails, change the grid.")
        grid = np.logspace(-0.1, 0.1, 10)
    bw_grid = np.reshape(grid, (-1, *np.ones(np.ndim(seed), dtype=int))) * seed
    if n_splits > len(data):
        n_splits = len(data)

    cv_scores = Parallel(n_jobs=-1, verbose=10)(
        delayed(_kde_cv_score)(bw, data, weights=weights, n_splits=n_splits)
        for bw in bw_grid
    )
    idx_best = np.argmax(cv_scores)

    if show:
        plt.plot(grid, cv_scores)
        plt.xlabel("Scaling factor")
        plt.ylabel("MLCV Figure of Merit (FoM)")
        plt.tight_layout()
        plt.show()
    if idx_best in (0, len(bw_grid) - 1):
        raise Exception(
            "Maximum not found in bw range selected. Move grid and try again."
        )
    bw_mlcv = bw_grid[idx_best]
    return bw_mlcv


def bw_auto(
    data,
    weights,
    grid=None,
    Nmlcv=1e4,
    Nbatch=1e4,
    Nknn=10,
    show=True,
    n_splits=10,
):

    """
    Automatic bandwidth optimization with kNN and MLCV methods combined.

    Parameters
    ----------
    data: array-like
        Array of samples. Must have shape (obs, dim).
    weights: array-like, optional
        Array of sample weights. By default all weights are 1.
    grid: array-like
        Grid to be used for the MLCV bandwidth optimization. Default:
        numpy.logspace(-1, 1, 10)
    Nmlcv: int
        Number of particles to use for the MLCV bandwidth optimization.
    Nbatch: int
        Number of particles per batch to use for the kNN bandwidth
        optimization.
    Nknn: int
        Number of closest neighbors used for the kNN bandiwdth
        optimization.
    show: bool
        Whether to show the scores plot of cross-validation.
    n_splits: int, optional
        Number of folds for cross-validation. Must be at least 2.
    """

    N, dim = data.shape

    # First fit with kNN to generate a seed adaptative bandwidth
    print("Fitting first with kNN method:")
    if Nbatch == -1:
        print("No batch size for kNN selected. Using all the particles list.")
        Nbatch = N
    BW_knn = bw_knn(data, weights=weights, k=Nknn, batch_size=Nbatch)
    print("")

    # Then fit with MLCV
    print("Fitting now with MLCV using the previous fitting as seed:")
    if Nmlcv == -1:
        print("No size for MLCV selected. Using all the particles list.")
        Nmlcv = N
    seed = BW_knn[:Nmlcv]
    print("If fitting takes too long, consider reducing Nmlcv.")
    BW_mlcv = bw_mlcv(
        data[:Nmlcv],
        weights=weights,
        seed=seed,
        grid=grid,
        show=show,
        n_splits=n_splits,
    )
    print("")

    print("Extending the MLCV optimization to full kNN bandwidth:")
    BW_knn_mlcv = BW_knn * BW_mlcv[0] / BW_knn[0]
    BW_knn_mlcv *= bw_silv(dim, len(BW_knn)) / bw_silv(dim, len(BW_mlcv))
    print("Done.")
    return BW_knn_mlcv


def optimize_bw(
    bw_method, vecs, ws=None, weightfun=None, maskfun=None, **kwargs
):
    """
    Optimize bandwidth with given method.

    Parameters
    ----------
    bw_method: str
        Bandwidth selection method. See bw_methods for available
        methods.
    vecs: array-like
        Array of samples. Must have shape (obs, dim).
    ws: array-like, optional
        Array of sample weights. By default all weights are 1.
    weightfun: function, optional
        Weighting function. If set, ws will be multiplied by
        weightfun(vecs).
    maskfun: function, optional
        Masking function. If set, vecs and ws will be replaced by
        vecs[maskfun(vecs)] and ws[maskfun(vecs)], respectively.
    **kwargs: optional
        Parameters to be passed to bandwidth selection method.

    Returns
    -------
    bw_opt: float or array-like
        Output of the selected bandwidth selection method.
    """
    vecs = np.array(vecs)
    if ws is not None:
        ws = np.array(ws)
        mask = ws > 0
    else:
        mask = np.ones(len(vecs), dtype=bool)
    if weightfun is not None:  # Apply weightfun
        ws = ws * weightfun(vecs)
    if maskfun is not None:  # Apply maskfun
        mask = np.logical_and(mask, maskfun(vecs))
    vecs = vecs[mask]
    ws = ws[mask]
    #
    if bw_method == "silv":  # Silverman's Rule
        dim = np.shape(vecs)[1]
        N_eff = len(vecs) if ws is None else np.sum(ws) ** 2 / np.sum(ws ** 2)
        return bw_silv(dim, N_eff)
    elif bw_method == "knn":  # K Nearest Neighbors method
        return bw_knn(vecs, weights=ws, **kwargs)
    elif bw_method == "mlcv":  # Maximum Likelihood Cross-Validation method
        return bw_mlcv(vecs, weights=ws, **kwargs)
    elif bw_method == "auto":  # kNN and MLCV combination
        return bw_auto(vecs, weights=ws, **kwargs)
    else:
        keys = list(bw_methods.keys())
        raise Exception("Invalid bw_method. Available: {}".format(keys))


bw_methods = {"auto": bw_auto, "silv": bw_silv, "knn": bw_knn, "mlcv": bw_mlcv}
