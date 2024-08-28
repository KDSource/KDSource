#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for KDSource object
"""

import os
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, parse, tostring

from KDEpy import TreeKDE

import matplotlib.colors as col
import matplotlib.pyplot as plt

import numpy as np

from .geom import Geometry
from .kde import bw_silv
from .kde import optimize_bw
from .plist import PList

R = {'g': 1 / (2 * np.sqrt(np.pi)),
     'e': 3 / 5,
     'b': 1 / 3}   # Roughness of the kernel

STD_DEFAULT = 1

varnames = ["ekin", "x", "y", "z", "dx", "dy", "dz", "t"]
varmap = {name: idx for idx, name in enumerate(varnames)}
units = ["MeV", "cm", "cm", "cm", "", "", "", "ms"]


def load(xmlfilename, N=-1):
    """
    Load KDSource from XML parameters file. After building the KDSource
    object, a fit is performed to load the particle list,
    leaving the KDSource ready to be evaluated. The bandwidths
    given in the XML file are not modified.

    Parameters
    ----------
    xmlfilename: str
        Name of XML parameters file
    N: int
        Number of particles to use for fitting. Default is -1, meaning
        that all particles in the list will be used.

    Returns
    -------
    kdsource: KDSource object
    """
    tree = parse(xmlfilename)
    root = tree.getroot()
    J = np.double(root.find("J").text)
    kelem = root.find("kernel")
    if kelem is not None:
        kern = kelem.text
    else:
        print("No kernel specified. Using gaussian as default.")
        kern = 'g'
    if kern == 'g':
        kern = 'gaussian'
    elif kern == 'e':
        kern = 'epa'
    elif kern == 'b':
        kern = 'box'
    plist = PList.load(root.find("PList"))
    geom = Geometry.load(root.find("Geom"))
    scaling = np.array(root.find("scaling").text.split(), dtype="float64")
    bwel = root.find("BW")
    if bool(int(bwel.attrib["variable"])):
        bw = np.fromfile(bwel.text, dtype="float32").astype("float64")
    else:
        bw = np.double(bwel.text)
    s = KDSource(plist, geom, bw, J=J, kernel=kern)
    s.fit(N=N, scaling=scaling)
    return s


class KDSource:
    def __init__(self, plist, geom, bw="silv", J=1.0, kernel="gaussian"):
        """
        Object representing Kernel Density Estimation (KDE) sources.

        A KDSource object is based on a particle list in MCPL format,
        wrapped in a PList object. It also includes a Geometry object,
        which defines variables (energy, position, direction) treatment.
        Having these basic components, KDSource can fit (optimize) a KDE
        model on the particle list, creating a distributional source.

        After optimization, analysis, plotting, etc., a KDSource can be
        exported in XML format, for later use in Python or other
        KDSource tools.

        Parameters
        ----------
        plist: PList object
            The PList wrapping the MCPL file containing the particle
            list.
        geom: Geometry object
            The Geometry defining particle variables treatment.
        bw: float, str or numpy.ndarray, optional
            Bandwidth or bandwidth selection method for KDE. If a float
            is passed, it is the bandwidth for all particles. If a
            numpy.ndarray is passed, it is the bandwidth for each
            particle. If a string is passed it is the bandwidth
            selection method. See bw_methods for available methods.
        J: float, optional
            The source total current, in [1/s]. If set, the density
            plots will have the correct units.
        kernel: string, optional
            The function kernel to fit the variables. Available options
            are 'gaussian', 'box', and 'epa'
        """
        self.plist = plist
        self.geom = geom
        self.bw_method = None
        if isinstance(bw, str):
            self.bw_method = bw
            bw = 1.0
        elif isinstance(bw, np.ndarray):
            if np.ndim(bw) >= 2:
                msg = "BW dimension must be < 2. \
                    Use scaling for anisotropic KDE."
                raise ValueError(msg)
        self.kde = TreeKDE(kernel=kernel, bw=bw)
        self.scaling = None
        self.J = J
        self.kernel = kernel
        self.R = R[self.kernel[0]]
        self.fitted = False

    def fit(self, N=-1, skip=0, scaling=None, **kwargs):
        """
        Fit KDE model to particle list.

        A number of particles are loaded from MCPL file to build the KDE
        model. If the KDSource object has a bandwidth selection method,
        it will be applied to optimize the bandwidth.

        Parameters
        ----------
        N: int
            Number of particles to use for fitting. The real number of
            particles used may be lower if end of particle list is
            reached or there are particles with zero weight. -1 to use
            all particles.
        skip: int
            Number of particles to skip in the list before starting to
            read.
        scaling: array-like, optional
            Scaling to be applied to each variable. This means each
            particle variable will be divided by the corresponding
            scaling element before applying KDE. By default, the
            standard deviation of each variable is used.
        **kwargs: optional
            Parameters to be passed to bandwidth selection method. Refer
            corresponding method for docs (see bw_methods for method
            names).
        """
        parts, ws = self.plist.get(N, skip)
        N = len(parts)
        if N == 0:
            raise Exception("No particles for fitting.")
        print("Using {} particles for fit.".format(N))
        vecs = self.geom.transform(parts)
        self.N_eff = np.sum(ws) ** 2 / np.sum(ws ** 2)
        if scaling is None:
            scaling = self.geom.std(vecs=vecs, weights=ws)
        else:
            scaling = np.array(scaling)
        scaling[scaling == 0] = STD_DEFAULT
        self.scaling = scaling
        if self.bw_method is not None:
            print("Calculating bw ... ")
            bw = optimize_bw(self.bw_method, vecs / self.scaling, ws, **kwargs)
            bw_opt = np.reshape(bw, (-1, 1)) * self.scaling
            print("Done\nOptimal bw ({}) = {}".format(self.bw_method, bw_opt))
            self.kde = TreeKDE(kernel=self.kernel, bw=bw)
        self.kde.fit(vecs / self.scaling, weights=ws)
        self.fitted = True

    def evaluate(self, parts):
        """
        Evaluate density and statistic error in a set of points.

        Variables and scaling treatment is so that the evaluated density
        units will be the inverse of the product of metrics volunits
        parameter, times [s-1] (e.g.: [s-1 MeV-1 cm-2 sr-1]).

        Parameters
        ----------
        parts: array-like
            Array of particles where evaluate density. Must have shape
            (obs, 8), with columns being [ekin, x, y, z, dx, dy, dz, t].

        Returns
        -------
        [evals, errs]: list of arrays
            Array of evaluated densities and estimated statistic errors.
        """
        if self.fitted is False:
            raise Exception("Must fit before evaluating.")
        vecs = self.geom.transform(parts)
        jacs = self.geom.jac(parts)
        evals = (
            1 / np.prod(self.scaling) * self.kde.evaluate(vecs / self.scaling)
        )
        errs = np.sqrt(
            evals
            * self.R ** self.geom.dim
            / (self.N_eff * np.mean(self.kde.bw) * np.prod(self.scaling))
        )
        evals *= self.J * jacs
        errs *= self.J * jacs
        return [evals, errs]

    def save(self, xmlfilename=None, bwfile=None, adjust_N=True):
        """
        Save KDSource object to XML parameters file.

        Parameters
        ----------
        xmlfilename: str
            Name of XML parameters file. By default it is set as:
            [MCPL file name]_source.xml
        bwfile: str or file object
            File name or file object to write bandwidths in binary
            float32 format, if a variable bandwidth is used. A file
            object with append mode can be used to save several
            KDSource's in same parameters file, useful when particle
            list is so big it can't be loaded at once.
        adjust_N: bool
            If True, before saving the following factor is applied on
            bandwidth:
            bw_silv(dim, N_tot) / bw_silv(dim, N)
            where N_tot is the total number of particles in the MCPL
            file, and N is the one passed to fit method. This way the
            bandwidth optimized for a subset of the particle list can be
            adapted to the full list.

        Returns
        -------
        xmlfilename: str
            Name of the written XML file.
        """
        # Process arguments
        if xmlfilename is None:
            xmlfilename = self.plist.filename.split(".")[0] + "_source.xml"
        if bwfile is None:
            bwfile = bwfilename = self.plist.filename.split(".")[0] + "_bws"
        elif isinstance(bwfile, str):
            bwfilename = bwfile
        else:  # Assume bwfile is file object
            bwfilename = bwfile.name
        bw = self.kde.bw
        if adjust_N:  # Adjust N_eff with Silverman factor
            dim = self.geom.dim
            if not self.plist.params_set:
                self.plist.set_params()
            bw *= bw_silv(dim, self.plist.N_eff) / bw_silv(dim, self.N_eff)
        # Build XML tree
        root = Element("KDSource")
        Jel = SubElement(root, "J")
        Jel.set("units", "1/s")
        Jel.text = str(self.J)
        kernelel = SubElement(root, "kernel")
        kernelel.text = str(self.kernel[0])
        pltree = SubElement(root, "PList")
        self.plist.save(pltree)
        gtree = SubElement(root, "Geom")
        self.geom.save(gtree)
        SubElement(root, "scaling").text = np.array_str(self.scaling)[1:-1]
        bwel = SubElement(root, "BW")
        if np.isscalar(bw):  # Constant bandwidth
            bwel.set("variable", "0")
            bwel.text = str(bw)
        else:  # Variable bandwidth
            bwel.set("variable", "1")
            bw.astype("float32").tofile(bwfile, format="float32")
            print("Bandwidth file: {}".format(bwfilename))
            bwel.text = os.path.abspath(bwfilename)
        # Write XML file
        xmlstr = tostring(root, encoding="utf8", method="xml")
        xmlstr = minidom.parseString(xmlstr).toprettyxml()
        with open(xmlfilename, "w") as file:
            file.write(xmlstr)
        print("Successfully saved parameters file {}".format(xmlfilename))
        return xmlfilename

    # Plot methods

    def plot_point(self, var, grid, part0, **kwargs):
        """
        1D plot of the full multivariate distribution.

        The current densities plotted are evaluated with the evaluate
        method, along a grid of particles, varying only one variable.

        Parameters
        ----------
        var: int or str
            Variable to be plotted, or its index. Names and indices of
            variables can be found in varnames list.
        grid: array-like
            1D array with values of var variable.
        part0: array-like
            Particle defining where to evaluate the rest of variables.
            Must have shape (7,), with variables ordered as in varnames
            list. Value of var variable will be ignored.
        **kwargs: optional
            Additional parameters for plotting options:

            xscale: 'linear' or 'log'
                Scale for x axis. Default: 'linear'
            yscale: 'linear' or 'log'. Default: 'log'
                Scale for y axis.
            fact: float
                Factor to apply on all densities. Default: 1
            label: str
                Line label in plot legend.

        Returns
        -------
        [fig, [scores, errs]]:
            Figure object, and evaluated densities and statistic errors.
        """
        if self.fitted is False:
            raise Exception("Must fit before plotting.")
        if isinstance(var, str):
            var = varmap[var]
        if "xscale" not in kwargs:
            kwargs["xscale"] = "linear"
        if "yscale" not in kwargs:
            kwargs["yscale"] = "log"
        if "label" not in kwargs:
            kwargs["label"] = "KDE"
        parts = np.zeros((len(grid), 7))
        parts[:, var] = grid
        part0[var] = 0
        parts += part0
        scores, errs = self.evaluate(parts)
        if "fact" in kwargs:
            scores *= kwargs["fact"]
            errs *= kwargs["fact"]

        plt.errorbar(
            grid, scores, errs, label=kwargs["label"], capsize=1, linewidth=1
        )
        plt.xscale(kwargs["xscale"])
        plt.yscale(kwargs["yscale"])
        plt.xlabel(r"${}\ [{}]$".format(varnames[var], units[var]))
        plt.ylabel(
            r"$J\ \left[ \frac{{{}}}{{{} s}} \right]$".format(
                self.plist.pt, self.geom.volunits
            )
        )
        plt.grid()
        plt.legend()
        plt.tight_layout()

        return [plt.gcf(), [scores, errs]]

    def plot_integr(self, var, grid, vec0=None, vec1=None, **kwargs):
        """
        1D plot of the univariate distribution along one variable.

        Densities are integrated over the specified range of the other
        variables, and evaluated on a newly created 1D KDE model. Units
        are the inverse of the units of the parametrized chosen variable
        (see geometry units parameter), times [s-1] (e.g.: [s-1 cm-1],
        [s-1 deg-1]).

        Parameters
        ----------
        var: int or str
            Variable to be plotted, or its index. Names and indices of
            variables can be found in varnames parameter of geometry.
            Note that this is a parametrized variable (e.g.: lethargy,
            theta, etc.).
        grid: array-like
            1D array with values of var variable.
        vec0: array-like
            Lower limits of integration range for other variables. Must
            be of shape (geom.dim,), with variables ordered as in
            varnames parameter of the geometry. Value of var variable
            will be ignored.
        vec1: array-like
            Upper limits of integration range for other variables. Same
            considerations as for vec0.
        **kwargs: optional
            Additional parameters for plotting options:

            xscale: 'linear' or 'log'
                Scale for x axis. Default: 'linear'
            yscale: 'linear' or 'log'. Default: 'log'
                Scale for y axis.
            fact: float
                Factor to apply on all densities. Default: 1
            label: str
                Line label in plot legend.
            adjust_bw: str
                If True, a Silverman factor will be applied to adjust
                the bandwidth to the dimension and number of samples
                used in plot. Default: False.

        Returns
        -------
        [fig, [scores, errs]]:
            Figure object, and evaluated densities and statistic errors.
        """
        if self.fitted is False:
            raise Exception("Must fit before plotting.")
        if isinstance(var, str):
            var = self.geom.varmap[var]
        if "xscale" not in kwargs:
            kwargs["xscale"] = "linear"
        if "yscale" not in kwargs:
            kwargs["yscale"] = "log"
        if "label" not in kwargs:
            kwargs["label"] = "KDE"
        if "adjust_bw" not in kwargs:
            kwargs["adjust_bw"] = False
        trues = np.ones(len(self.kde.data), dtype=bool)
        if vec0 is not None:
            mask1 = np.logical_and.reduce(
                vec0 / self.scaling <= self.kde.data, axis=1
            )
        else:
            mask1 = trues
        if vec1 is not None:
            mask2 = np.logical_and.reduce(
                self.kde.data <= vec1 / self.scaling, axis=1
            )
        else:
            mask2 = trues
        mask = np.logical_and(mask1, mask2)
        vecs = self.kde.data[:, var][mask].reshape(-1, 1)
        ws = self.kde.weights[mask]
        N_eff = np.sum(ws) ** 2 / np.sum(ws ** 2)
        scaling = self.scaling[var]
        bw = self.kde.bw
        if kwargs["adjust_bw"]:
            bw *= bw_silv(1, N_eff) / bw_silv(self.geom.dim, self.N_eff)
        kde = TreeKDE(kernel=self.kernel, bw=bw)
        kde.fit(vecs, weights=ws)
        scores = 1 / scaling * kde.evaluate(grid.reshape(-1, 1) / scaling)
        errs = np.sqrt(scores * self.R / (N_eff * np.mean(bw) * scaling))
        scores *= self.J * np.sum(ws) / np.sum(self.kde.weights)
        errs *= self.J * np.sum(ws) / np.sum(self.kde.weights)
        if "fact" in kwargs:
            scores *= kwargs["fact"]
            errs *= kwargs["fact"]

        plt.errorbar(
            grid, scores, errs, label=kwargs["label"], capsize=1, linewidth=1
        )
        plt.xscale(kwargs["xscale"])
        plt.yscale(kwargs["yscale"])
        plt.xlabel(
            r"${}\ [{}]$".format(self.geom.varnames[var], self.geom.units[var])
        )
        plt.ylabel(
            r"$J\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(
                self.plist.pt, self.geom.units[var]
            )
        )
        plt.grid()
        plt.legend()
        plt.tight_layout()

        return [plt.gcf(), [scores, errs]]

    def plot_E(self, grid_E, vec0=None, vec1=None, **kwargs):
        """
        1D plot of the energy spectrum.

        Same as plot("ekin", grid_E), but applying jacobian of energy
        parametrization (if any). When using Lethargy metric, plot
        method will give plots with units [s-1] (lethargy is
        dimesionless), while this method will produce the expected
        spectrum with units [s-1 MeV-1].

        Parameters
        ----------
        grid_E: array-like
            1D array with values of ekin.
        vec0: array-like
            Lower limits of integration range for other variables. Must
            be of shape (geom.dim,), with variables ordered as in
            varnames parameter of the geometry. Value of energy variable
            will be ignored.
        vec1: array-like
            Upper limits of integration range for other variables. Same
            considerations as for vec0.
        **kwargs: optional
            Additional parameters for plotting options:

            fact: float
                Factor to apply on all densities. Default: 1
            label: str
                Line label in plot legend.
            adjust_bw: str
                If True, a Silverman factor will be applied to adjust
                the bandwidth to the dimension and number of samples
                used in plot. Default: False.

        Returns
        -------
        [fig, [scores, errs]]:
            Figure object, and evaluated densities and statistic errors.
        """
        if self.fitted is False:
            raise Exception("Must fit before plotting.")
        if "label" not in kwargs:
            kwargs["label"] = "KDE"
        if "adjust_bw" not in kwargs:
            kwargs["adjust_bw"] = False
        trues = np.ones(len(self.kde.data), dtype=bool)
        if vec0 is not None:
            mask1 = np.logical_and.reduce(
                vec0 / self.scaling <= self.kde.data, axis=1
            )
        else:
            mask1 = trues
        if vec1 is not None:
            mask2 = np.logical_and.reduce(
                self.kde.data <= vec1 / self.scaling, axis=1
            )
        else:
            mask2 = trues
        mask = np.logical_and(mask1, mask2)
        if np.sum(mask) == 0:
            raise Exception("No particles in [vec0,vec1] range.")
        # This function assumes that energy is the first defined metric
        vecs = self.kde.data[:, 0][mask].reshape(-1, 1)
        ws = self.kde.weights[mask]
        N_eff = np.sum(ws) ** 2 / np.sum(ws ** 2)
        scaling = self.scaling[0]
        bw = self.kde.bw
        if kwargs["adjust_bw"]:
            bw *= bw_silv(1, N_eff) / bw_silv(self.geom.dim, self.N_eff)
        kde = TreeKDE(kernel=self.kernel, bw=bw)
        kde.fit(vecs, weights=ws)
        grid = self.geom.ms[0].transform(grid_E)
        jacs = self.geom.ms[0].jac(grid_E)
        scores = 1 / scaling * kde.evaluate(grid.reshape(-1, 1) / scaling)
        errs = np.sqrt(scores * self.R / (N_eff * np.mean(bw) * scaling))
        scores *= self.J * np.sum(ws) / np.sum(self.kde.weights) * jacs
        errs *= self.J * np.sum(ws) / np.sum(self.kde.weights) * jacs
        if "fact" in kwargs:
            scores *= kwargs["fact"]
            errs *= kwargs["fact"]

        plt.errorbar(
            grid_E, scores, errs, label=kwargs["label"], capsize=1, linewidth=1
        )
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"$E\ [MeV]$")
        plt.ylabel(
            r"$J\ \left[ \frac{{{}}}{{MeV\ s}} \right]$".format(self.plist.pt)
        )
        plt.grid()
        plt.legend()
        plt.tight_layout()

        return [plt.gcf(), [scores, errs]]

    def plot_t(self, grid_t, vec0=None, vec1=None, **kwargs):
        """
        1D plot of the time distribution.

        Same as plot("t", grid_t), but applying jacobian of time
        parametrization (if any). When using Decade metric, plot
        method will give plots with units [s-1] (decade is
        dimesionless), while this method will produce the expected
        spectrum with units [s-2].

        Parameters
        ----------
        grid_t: array-like
            1D array with values of t.
        vec0: array-like
            Lower limits of integration range for other variables. Must
            be of shape (geom.dim,), with variables ordered as in
            varnames parameter of the geometry. Value of energy variable
            will be ignored.
        vec1: array-like
            Upper limits of integration range for other variables. Same
            considerations as for vec0.
        **kwargs: optional
            Additional parameters for plotting options:

            fact: float
                Factor to apply on all densities. Default: 1
            label: str
                Line label in plot legend.
            adjust_bw: str
                If True, a Silverman factor will be applied to adjust
                the bandwidth to the dimension and number of samples
                used in plot. Default: False.

        Returns
        -------
        [fig, [scores, errs]]:
            Figure object, and evaluated densities and statistic errors.
        """
        if self.fitted is False:
            raise Exception("Must fit before plotting.")
        if "label" not in kwargs:
            kwargs["label"] = "KDE"
        if "adjust_bw" not in kwargs:
            kwargs["adjust_bw"] = False
        trues = np.ones(len(self.kde.data), dtype=bool)
        if vec0 is not None:
            mask1 = np.logical_and.reduce(
                vec0 / self.scaling <= self.kde.data, axis=1
            )
        else:
            mask1 = trues
        if vec1 is not None:
            mask2 = np.logical_and.reduce(
                self.kde.data <= vec1 / self.scaling, axis=1
            )
        else:
            mask2 = trues
        mask = np.logical_and(mask1, mask2)
        if np.sum(mask) == 0:
            raise Exception("No particles in [vec0,vec1] range.")
        # This function assumes that time is the last defined metric
        vecs = self.kde.data[:, -1][mask].reshape(-1, 1)
        ws = self.kde.weights[mask]
        N_eff = np.sum(ws) ** 2 / np.sum(ws ** 2)
        scaling = self.scaling[-1]
        bw = self.kde.bw
        if kwargs["adjust_bw"]:
            bw *= bw_silv(1, N_eff) / bw_silv(self.geom.dim, self.N_eff)
        kde = TreeKDE(kernel=self.kernel, bw=bw)
        kde.fit(vecs, weights=ws)
        grid = self.geom.ms[-1].transform(grid_t)
        jacs = self.geom.ms[-1].jac(grid_t)
        scores = 1 / scaling * kde.evaluate(grid.reshape(-1, 1) / scaling)
        errs = np.sqrt(scores * self.R / (N_eff * np.mean(bw) * scaling))
        scores *= self.J * np.sum(ws) / np.sum(self.kde.weights) * jacs
        errs *= self.J * np.sum(ws) / np.sum(self.kde.weights) * jacs
        if "fact" in kwargs:
            scores *= kwargs["fact"]
            errs *= kwargs["fact"]

        plt.errorbar(
            grid_t, scores, errs, label=kwargs["label"], capsize=1, linewidth=1
        )
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"$t$ [ms]")
        plt.ylabel(
            r"$J\ \left[ \frac{{{}}}{{s^2}} \right]$".format(self.plist.pt)
        )
        plt.grid()
        plt.legend()
        plt.tight_layout()

        return [plt.gcf(), [scores, errs]]

    def plot2D_point(self, vrs, grids, part0, **kwargs):
        """
        2D plot of the full multivariate distribution.

        The current densities plotted are evaluated with the evaluate
        method, along a grid of particles, varying only two variable.

        Parameters
        ----------
        vrs: list
            List of the 2 variables to be plotted, or its indices. Names
            and indices of variables can be found in varnames list.
        grids: list
            List of the 2 1D arrays with values of vrs variables.
        part0: array-like
            Particle defining where to evaluate the rest of variables.
            Must have shape (7,), with variables ordered as in varnames
            list. Value of vrs variables will be ignored.
        **kwargs: optional
            Additional parameters for plotting options:

            scale: 'linear' or 'log'
                Scale for color map. Default: 'log'
            fact: float
                Factor to apply on all densities. Default: 1
            title: str
                Plot title.

        Returns
        -------
        [fig, [scores, errs]]:
            Figure object, and evaluated densities and statistic errors.
        """
        if self.fitted is False:
            raise Exception("Must fit before plotting.")
        if isinstance(vrs[0], str):
            vrs[0] = self.geom.varmap[vrs[0]]
        if isinstance(vrs[1], str):
            vrs[1] = self.geom.varmap[vrs[1]]
        if "scale" not in kwargs:
            kwargs["scale"] = "linear"
        if "title" not in kwargs:
            kwargs["title"] = "KDE"
        parts = np.zeros((len(grids[0]) * len(grids[1]), 7))
        parts[:, vrs] = np.reshape(np.meshgrid(*grids), (2, -1)).T
        part0 = np.array(part0)
        part0[vrs] = 0
        parts += part0
        scores, errs = self.J * self.evaluate(parts)
        if "fact" in kwargs:
            scores *= kwargs["fact"]
            errs *= kwargs["fact"]

        xx = np.concatenate(
            (grids[0][:1], (grids[0][1:] + grids[0][:-1]) / 2, grids[0][-1:])
        )
        yy = np.concatenate(
            (grids[1][:1], (grids[1][1:] + grids[1][:-1]) / 2, grids[1][-1:])
        )
        if kwargs["scale"] == "log":
            norm = col.LogNorm()
        else:
            norm = None
        plt.pcolormesh(
            xx,
            yy,
            scores.reshape(len(grids[1]), len(grids[0])),
            cmap="jet",
            norm=norm,
            rasterized=True,
        )
        plt.colorbar()
        title = kwargs["title"]
        title += "\n" + r"$J\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(
            self.plist.pt, self.geom.volunits
        )
        plt.title(title)
        plt.xlabel(r"${}\ [{}]$".format(varnames[vrs[0]], units[vrs[0]]))
        plt.ylabel(r"${}\ [{}]$".format(varnames[vrs[1]], units[vrs[1]]))
        plt.tight_layout()

        return [plt.gcf(), [scores, errs]]

    def plot2D_integr(self, vrs, grids, vec0=None, vec1=None, **kwargs):
        """
        2D plot of the bivariate distribution along two variables.

        Densities are integrated over the specified range of the other
        variables, and evaluated on a newly created 2D KDE model. Units
        are the inverse of the units of the parametrized chosen
        variables (see geometry units parameter), times [s-1] (e.g.:
        [s-1 cm-2], [s-1 deg-2]).

        Parameters
        ----------
        vrs: list
            List of the 2 variables to be plotted, or its indices. Names
            and indices of variables can be found in varnames parameter
            of geometry. Note that these are parametrized variables
            (e.g.: lethargy, theta, etc.).
        grids: list
            1D array with values of vrs variable.
        vec0: array-like
            Lower limits of integration range for other variables. Must
            be of shape (geom.dim,), with variables ordered as in
            varnames parameter of the geometry. Value of vrs variable
            will be ignored.
        vec1: array-like
            Upper limits of integration range for other variables. Same
            considerations as for vec0.
        **kwargs: optional
            Additional parameters for plotting options:

            scale: 'linear' or 'log'
                Scale for color map. Default: 'log'
            fact: float
                Factor to apply on all densities. Default: 1
            title: str
                Plot title.
            adjust_bw: str
                If True, a Silverman factor will be applied to adjust
                the bandwidth to the dimension and number of samples
                used in plot. Default: False.

        Returns
        -------
        [fig, [scores, errs]]:
            Figure object, and evaluated densities and statistic errors.
        """
        if self.fitted is False:
            raise Exception("Must fit before plotting.")
        if isinstance(vrs[0], str):
            vrs[0] = self.geom.varmap[vrs[0]]
        if isinstance(vrs[1], str):
            vrs[1] = self.geom.varmap[vrs[1]]
        if "scale" not in kwargs:
            kwargs["scale"] = "linear"
        if "title" not in kwargs:
            kwargs["title"] = "KDE"
        if "adjust_bw" not in kwargs:
            kwargs["adjust_bw"] = False
        trues = np.array(len(self.kde.data) * [True])
        if vec0 is not None:
            mask1 = np.logical_and.reduce(
                vec0 / self.scaling <= self.kde.data, axis=1
            )
        else:
            mask1 = trues
        if vec1 is not None:
            mask2 = np.logical_and.reduce(
                self.kde.data <= vec1 / self.scaling, axis=1
            )
        else:
            mask2 = trues
        mask = np.logical_and(mask1, mask2)
        if np.sum(mask) == 0:
            raise Exception("No particles in [vec0,vec1] range.")
        vecs = self.kde.data[:, vrs][mask]
        ws = self.kde.weights[mask]
        N_eff = np.sum(ws) ** 2 / np.sum(ws ** 2)
        scaling = self.scaling[vrs]
        bw = self.kde.bw
        if kwargs["adjust_bw"]:
            bw *= bw_silv(1, N_eff) / bw_silv(self.geom.dim, self.N_eff)
        kde = TreeKDE(kernel=self.kernel, bw=bw)
        kde.fit(vecs, weights=ws)
        grid = np.reshape(np.meshgrid(*grids), (2, -1)).T
        scores = 1 / np.prod(scaling) * kde.evaluate(grid / scaling)
        errs = np.sqrt(
            scores * self.R ** 2 / (N_eff * np.mean(bw) * np.prod(scaling))
        )
        scores *= self.J * np.sum(ws) / np.sum(self.kde.weights)
        errs *= self.J * np.sum(ws) / np.sum(self.kde.weights)
        if "fact" in kwargs:
            scores *= kwargs["fact"]
            errs *= kwargs["fact"]

        xx = np.concatenate(
            (grids[0][:1], (grids[0][1:] + grids[0][:-1]) / 2, grids[0][-1:])
        )
        yy = np.concatenate(
            (grids[1][:1], (grids[1][1:] + grids[1][:-1]) / 2, grids[1][-1:])
        )
        if kwargs["scale"] == "log":
            norm = col.LogNorm()
        else:
            norm = None
        plt.pcolormesh(
            xx,
            yy,
            scores.reshape(len(grids[1]), len(grids[0])),
            cmap="jet",
            norm=norm,
            rasterized=True,
        )
        plt.colorbar()
        if self.geom.units[vrs[0]] == self.geom.units[vrs[1]]:
            units = self.geom.units[vrs[0]] + "^2"
        else:
            units = self.geom.units[vrs[0]] + self.geom.units[vrs[1]]
        title = kwargs["title"]
        title += "\n" + r"$J\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(
            self.plist.pt, units
        )
        plt.title(title)
        plt.xlabel(
            r"${}\ [{}]$".format(
                self.geom.varnames[vrs[0]], self.geom.units[vrs[0]]
            )
        )
        plt.ylabel(
            r"${}\ [{}]$".format(
                self.geom.varnames[vrs[1]], self.geom.units[vrs[1]]
            )
        )
        plt.tight_layout()

        return [plt.gcf(), [scores, errs]]


def run(plist, geom, kernel='gaussian', skip=0, importance=None,
        N=-1, Nmlcv=-1, Nbatch=-1, Nknn=10, grid=None,
        show=True):

    """
    Run KDSource to fit model to particle list with kNN and MLCV methods.

    Parameters
    ----------
    plist: PList object
        The PList wrapping the MCPL file containing the particle
        list.
    geom: Geometry object
        The Geometry defining particle variables treatment.
    kernel: string, optional
        The function kernel to fit the variables. Available options
        are 'gaussian', 'box', and 'epa'
    skip: int
        Number of particles to skip in the list before starting to
        read.
    importance: array-like, optional
        Scaling to be applied to each variable. This means each
        particle variable standard deviation will be divided by
        the corresponding scaling element before applying KDE.
        By default, the importance is set to 1, i.e.
        the standard deviation of each variable is used.
    N: int
        Number of particles to use for fitting. The real number of
        particles used may be lower if end of particle list is
        reached or there are particles with zero weight. -1 to use
        all particles.
    Nmlcv: int
        Number of particles to use for the MLCV bandwidth optimization.
    Nbatch: int
        Number of particles per batch to use for the kNN bandwidth
        optimization.
    Nknn: int
        Number of closest neighbors used for the kNN bandiwdth
        optimization.
    grid: array-like
        Grid to be used for the MLCV bandwidth optimization.
    **kwargs: optional
        Parameters to be passed to bandwidth selection method. Refer
        corresponding method for docs (see bw_methods for method
        names).
    """

    # Create KDSource
    s = KDSource(plist, geom, kernel)

    # Get the length of the particles list
    parts, ws = s.plist.get(N=-1)
    N = len(parts)
    print('Using N={:d} particles from the original particles list.\n'.format(
        N))

    # Use the standard deviation as scaling factor
    scaling = s.geom.std(parts=parts)

    # If an importance is specified, apply to the scaling factors
    if importance is None:
        print('No importance selected. Using 1 for each variable.\n')
        importance = [1] * s.geom.dim

    scaling /= importance

    # First fit with kNN to generate a seed adaptative bandwidth
    print('Fitting first with kNN method:')
    s.bw_method = "knn"
    if Nbatch == -1:
        Nbatch = int(N / Nknn)
    s.fit(N, skip=skip, scaling=scaling, batch_size=Nbatch, k=Nknn)
    bw_knn = s.kde.bw
    print('')

    # Then fit with MLCV
    print('Fitting now with MLCV using the previous fitting as seed:')
    s.bw_method = "mlcv"
    if Nmlcv == -1:
        print("No size for MLCV selected. Using all the particles list.")
        Nmlcv = N
    if grid is None:
        print("No grid specified. Using np.logspace(-0.1, 0.1, 10).")
        print("If fitting fails, change the grid.")
        grid = np.logspace(-0.1, 0.1, 10)
    seed = bw_knn[:Nmlcv]
    print("If fitting takes too long, consider reducing Nmlcv.")
    s.fit(Nmlcv, scaling=scaling, seed=seed, grid=grid, show=show)
    bw_mlcv = s.kde.bw
    print('')

    print('Extending the MLCV optimization to full kNN bandwidth:')
    bw_knn_mlcv = bw_knn * bw_mlcv[0] / bw_knn[0]
    dim = s.geom.dim
    bw_knn_mlcv *= bw_silv(dim, len(bw_knn)) / bw_silv(dim, len(bw_mlcv))
    s = KDSource(plist, geom, bw=bw_knn_mlcv, kernel=kernel)
    s.fit(N=N, scaling=scaling)
    print('Done.')

    return s
