#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for statistic analysis
"""

import numpy as np
import matplotlib.pyplot as plt

def apply_weight_mask(vecs, ws, weightfun=None, maskfun=None):
    """
    Apply weighting and masking functions to particle list.

    Parameters
    ---------
    vecs: array-like
        Array of particles. Must have shape (N, dim).
    ws: array-like, optional
        Array of particle weights.
    weightfun: function, optional
        Weighting function. If set, ws will be multiplied by
        weightfun(vecs).
    maskfun: function, optional
        Masking function. If set, vecs and ws will be replaced by
        vecs[maskfun(vecs)] and ws[maskfun(vecs)], respectively.

    Returns
    -------
    [vecs, ws]: list
        Particle list and weights after applying weighting and masking
        functions.
    """
    if weightfun is not None:
        ws = ws*weightfun(vecs)
    mask = (ws > 0)
    if maskfun is not None:
        mask = np.logical_and(mask, maskfun(vecs))
    return [vecs[mask,:], ws[mask]]

def convergence(vecs, ws, param, fracmin=0.1, steps=10, plot=True):
    """
    Compute statistical parameter over subsets of particle list.

    A number of subsets of particle list are built, with linearly growing
    size and the last one being the full particle list. For each subset the
    specified statistic parameter is computed.

    Parameters
    ---------
    vecs: array-like
        Array of particles. Must have shape (N, dim).
    ws: array-like, optional
        Array of particle weights.
    param: function
        Function that computes the statistic parameter. Must have signature:
            param(vecs, ws) -> [param, err]
    fracmin: float, optional
        Size of the first subset, as a fraction of the full particle list
        size.
    steps: int, optional
        Number of subsets to build and compute statistic parameter.
    plot: bool, optional
        Whether to show the plot of statistic parameter vs. subset size.

    Returns
    -------
    [Ns, params, errs]: list
        Subsets sizes, statistic parameter values, and errors.
    """
    Nmax = len(ws)
    Nmin = int(fracmin * Nmax)
    Ns = np.rint(np.linspace(Nmin, Nmax, num=steps+1)).astype(int)
    Ns = Ns[Ns>0]
    params = []
    errs = []
    for N in Ns:
        param_err = param(vecs[:N],ws[:N])
        params.append(param_err[0])
        errs.append(param_err[1])
    params = np.array(params)
    errs = np.array(errs)
    if Ns[0]!=Ns[-1] and plot:
        plt.plot(Ns, params, 'o-')
        plt.fill_between(Ns, params-errs, params+errs, color='blue', alpha=0.3)
        plt.xlabel("Number of particles")
    return [Ns, params, errs]

def mean_weight(vecs, ws):
    """Compute mean weight and its statistic error."""
    return [ws.mean(), ws.std()/np.sqrt(len(ws))]
def mean(vecs, ws, var):
    """Compute mean particle and its statistic error."""
    mean = np.average(vecs[:,var], weights=ws)
    err = np.sqrt(np.average((vecs[:,var]-mean)**2, weights=ws)/np.sum(ws))
    return [mean, err]
def std(vecs, ws, var):
    """Compute particles standard deviation and its statistic error."""
    mean = np.average(vecs[:,var], weights=ws)
    s2 = np.average((vecs[:,var]-mean)**2, weights=ws)
    std = np.sqrt(s2)
    varerr = np.sqrt(np.average(((vecs[:,var]-mean)**2-s2)**2, weights=ws)/np.sum(ws))
    stderr = varerr * 0.5 / std
    return [std, stderr]

class Stats:
    def __init__(self, vecs, ws, weightfun=None, maskfun=None):
        """
        Object for statistic analysis of particle list.

        This class methods show the variation and convergence of some statistic
        parameters as the number of particles in the list grows. This can be
        used to determine whether the list size is enough or not, based on some
        user-defined criterion.

        Weighting and masking functions allow using an importance function based
        on particle parameters, and selecting a region of the phase-space to
        analyze.

        Parameters
        ----------
        vecs: array-like
            Array of particles. Must have shape (N, dim).
        ws: array-like, optional
            Array of particle weights.
        weightfun: function, optional
            Weighting function. If set, ws will be multiplied by
            weightfun(vecs).
        maskfun: function, optional
            Masking function. If set, vecs and ws will be replaced by
            vecs[maskfun(vecs)] and ws[maskfun(vecs)], respectively.
        """
        if len(vecs) != len(ws):
            raise ValueError("vecs and ws must have same len.")
        self.vecs, self.ws = apply_weight_mask(vecs, ws, weightfun, maskfun)
        if len(ws) == 0:
            raise Exception("Empty particle list.")
        self.dim = vecs.shape[1]
        self.N = len(ws)
    def mean_weight(self, fracmin=0.1, steps=10, plot=True):
        """
        Compute convergence of mean weight.

        Parameters
        ----------
        fracmin: float, optional
            Size of the first subset, as a fraction of the full particle list
            size.
        steps: int, optional
            Number of subsets to build and compute mean weight.
        plot: bool, optional
            Whether plot mean weight vs. subset size.

        Returns
        -------
        [Ns, params, errs]: list
            Subsets sizes, mean weight values, and errors.
        """
        Ns,params,errs = convergence(self.vecs, self.ws, mean_weight, fracmin=fracmin, steps=steps, plot=plot)
        if(plot):
            plt.ylabel("Mean weight")
        return [Ns, params, errs]
    def mean(self, var, varname=None, fracmin=0.1, steps=10, plot=True):
        """
        Compute convergence of a variable mean.

        Parameters
        ----------
        var: int
            Index of variable to compute mean.
        varname: str, optional
            Selected variable name for display in plot.
        fracmin: float, optional
            Size of the first subset, as a fraction of the full particle list
            size.
        steps: int, optional
            Number of subsets to build and compute variable mean.
        plot: bool, optional
            Whether plot variable mean vs. subset size.

        Returns
        -------
        [Ns, params, errs]: list
            Subsets sizes, variable mean values, and errors.
        """
        mean_var = lambda vecs,ws: mean(vecs,ws,var)
        Ns,params,errs = convergence(self.vecs, self.ws, mean_var, fracmin=fracmin, steps=steps, plot=plot)
        if(plot):
            if varname is None:
                varname = "v%d"%var
            plt.ylabel("{} mean".format(varname))
        return [Ns, params, errs]
    def std(self, var, varname=None, fracmin=0.1, steps=10, plot=True):
        """
        Compute convergence of a variable standard deviation.

        Parameters
        ----------
        var: int
            Index of variable to compute standard deviation.
        varname: str, optional
            Selected variable name for display in plot.
        fracmin: float, optional
            Size of the first subset, as a fraction of the full particle list
            size.
        steps: int, optional
            Number of subsets to build and compute variable standard deviation.
        plot: bool, optional
            Whether plot variable standard deviation vs. subset size.

        Returns
        -------
        [Ns, params, errs]: list
            Subsets sizes, variable standard deviation values, and errors.
        """
        std_var = lambda vecs,ws: std(vecs,ws,var)
        Ns,params,errs = convergence(self.vecs, self.ws, std_var, fracmin=fracmin, steps=steps, plot=plot)
        if(plot):
            if varname is None:
                varname = "v%d"%var
            plt.ylabel("{} standard deviation".format(varname))
        return [Ns, params, errs]
