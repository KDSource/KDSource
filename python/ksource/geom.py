#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for Geometry and Metric objects
"""

from xml.etree import ElementTree as ET

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class Metric:
    def __init__(self, partvars, varnames, units, volunits):
        """
        Abstract object defining metrics for a subset of variables.

        The main function of a Metric is definig a transformation from a
        subset of particle variables to certain parametrized variables,
        which can be more suitable for applying KDE.

        See _metrics for available metrics.

        Parameters
        ----------
        partvars: list of int
            Indices of the particle variables to parametrize (see
            varnames global list).
        varnames: list of str
            Names of parametrized variables.
        units: list of str
            Units of each parametrized variable.
        volunits: str
            Units of the product of variables to parametrize.
        """
        self.partvars = partvars
        if len(varnames) != len(units):
            raise ValueError("varnames and units must have same len.")
        self.dim = len(varnames)
        self.varnames = varnames
        self.varmap = {name:idx for idx,name in enumerate(varnames)}
        self.units = units
        self.volunits = volunits
    def transform(self, parts):
        """Transform particle variables to parametrized variables."""
        return parts
    def inverse_transform(self, vecs):
        """Transform parametrized variables to particle variables."""
        return vecs
    def jac(self, parts):
        """Jacobian of transformation."""
        return np.ones(len(parts))
    def mean(self, parts=None, vecs=None, weights=None):
        """
        Mean of particle variables.

        Mean is computed in parametrized space, and transformed back to
        particle variables.

        Parameters
        ----------
        parts: array-like, optional
            Array of particle variables.
        vecs: array-like, optional
            Array of parametrized variables. If set, overrides parts.
        weights: array-like, optional
            Array of particle statistic weights.
        """
        if vecs is None:
            vecs = self.transform(parts)
        return self.inverse_transform(np.average(vecs, axis=0, weights=weights))
    def std(self, parts=None, vecs=None, weights=None):
        """
        Standard deviation of particle variables.

        Standard deviation is computed in parametrized space, and
        transformed back to particle variables.

        Parameters
        ----------
        parts: array-like, optional
            Array of particle variables.
        vecs: array-like, optional
            Array of parametrized variables. If set, overrides parts.
        weights: array-like, optional
            Array of particle statistic weights.
        """
        if vecs is None:
            vecs = self.transform(parts)
        mn = np.average(vecs, axis=0, weights=weights)
        return np.sqrt(np.average((vecs-mn)**2, axis=0, weights=weights))
    def save(self, mtree):
        """Save Metric parameters into XML tree."""
        ET.SubElement(mtree, "dim").text = str(self.dim)
        ET.SubElement(mtree, "params").set("nps", "0")
    @staticmethod
    def load(mtree):
        """Load parameters from XML tree and build Metric."""
        raise Exception("Load method not implemented.")

class Geometry (Metric):
    def __init__(self, metrics, trasl=None, rot=None):
        """
        Object defining particle variables treatment (metrics).

        The main function of a Geometry is definig a transformation from
        particle variables to certain parametrized variables, which can
        be more suitable for applying KDE.

        Parameters
        ----------
        metrics: list
            List of metrics for each subset of variables. They must
            cover all particle variables (energy, position and
            direction).
        trasl: array-like, optional
            Spatial traslation for source. Default is no traslation.
        rot: numpy.ndarray or scipy.spatial.transform.Rotation, optional
            Spatial rotation for source. Can be a scipy Rotation object,
            or any array-like which can be used to generate a scipy
            Rotation object (rotation vector, quaternion or rotation
            matrix). Rotation is applied after traslation in transform,
            and before traslation in inverse_transform. Default is no
            rotation.
        """
        partvars = range(7)
        varnames = sum([metric.varnames for metric in metrics], [])
        units = sum([metric.units for metric in metrics], [])
        volunits = "".join([metric.volunits+" " for metric in metrics])[:-1]
        super().__init__(partvars, varnames, units, volunits)
        self.ms = metrics
        if trasl is not None:
            trasl = np.array(trasl).reshape(-1)
            if trasl.shape != (3,):
                raise ValueError("Invalid trasl.")
        if rot is not None:
            if not isinstance(rot, st.Rotation):
                rot = np.array(rot)
                if rot.shape == (4,):
                    rot = st.Rotation.from_quat(rot)
                elif rot.shape == (3,3):
                    rot = st.Rotation.from_matrix(rot)
                elif rot.shape == (3,):
                    rot = st.Rotation.from_rotvec(rot)
                else:
                    raise ValueError("Invalid rot.")
        self.trasl = trasl
        self.rot = rot
    def transform(self, parts):
        """Transform particle variables to parametrized variables."""
        if self.trasl is not None:
            parts[:,1:4] -= self.trasl # Position
        if self.rot is not None:
            parts[:,1:4] = self.rot.apply(parts[:,1:4], inverse=True) # Position
            parts[:,4:7] = self.rot.apply(parts[:,4:7], inverse=True) # Direction
        vecss = []
        for metric in self.ms:
            vecss.append(metric.transform(parts[:,metric.partvars]))
        return np.concatenate(vecss, axis=1)
    def inverse_transform(self, vecs):
        """Transform parametrized variables to particle variables."""
        parts = np.zeros((len(vecs), 7))
        end = 0
        for metric in self.ms:
            start = end
            end = start + metric.dim
            parts[:,metric.partvars] = metric.inverse_transform(vecs[:,start:end])
        if self.trasl is not None:
            parts[:,1:4] += self.trasl # Position
        if self.rot is not None:
            parts[:,1:4] = self.rot.apply(parts[:,1:4]) # Position
            parts[:,4:7] = self.rot.apply(parts[:,4:7]) # Direction
        return parts
    def jac(self, parts):
        """Jacobian of transformation."""
        if self.trasl is not None:
            parts[:,1:4] -= self.trasl # Position
        if self.rot is not None:
            parts[:,1:4] = self.rot.apply(parts[:,1:4], inverse=True) # Position
            parts[:,4:7] = self.rot.apply(parts[:,4:7], inverse=True) # Direction
        jacs = []
        for metric in self.ms:
            jacs.append(metric.jac(parts[:,metric.partvars]))
        return np.prod(jacs, axis=1)
    def mean(self, parts=None, vecs=None, weights=None):
        """
        Mean of particle variables.

        Mean is computed in parametrized space, and transformed back to
        particle variables.

        Parameters
        ----------
        parts: array-like, optional
            Array of particle variables.
        vecs: array-like, optional
            Array of parametrized variables. If set, overrides parts.
        weights: array-like, optional
            Array of particle statistic weights.
        """
        if vecs is None:
            vecs = self.transform(parts)
        means = []
        end = 0
        for metric in self.ms:
            start = end
            end = start + metric.dim
            means.append(metric.mean(vecs=vecs[:,start:end], weights=weights))
        return np.concatenate(means)
    def std(self, parts=None, vecs=None, weights=None):
        """
        Standard deviation of particle variables.

        Standard deviation is computed in parametrized space, and
        transformed back to particle variables.

        Parameters
        ----------
        parts: array-like, optional
            Array of particle variables.
        vecs: array-like, optional
            Array of parametrized variables. If set, overrides parts.
        weights: array-like, optional
            Array of particle statistic weights.
        """
        if vecs is None:
            vecs = self.transform(parts)
        stds = []
        end = 0
        for metric in self.ms:
            start = end
            end = start + metric.dim
            stds.append(metric.std(vecs=vecs[:,start:end], weights=weights))
        return np.concatenate(stds)
    def save(self, gtree):
        """Save Geometry parameters into XML tree."""
        gtree.set("order", str(len(self.ms)))
        for metric in self.ms:
            mtree = ET.SubElement(gtree, metric.__class__.__name__)
            metric.save(mtree)
        trasl = np.array_str(self.trasl)[1:-1] if self.trasl is not None else ""
        ET.SubElement(gtree, "trasl").text = trasl
        rot = np.array_str(self.rot.as_rotvec())[1:-1] if self.rot is not None else ""
        ET.SubElement(gtree, "rot").text = rot
    @staticmethod
    def load(gtree):
        """Load parameters from XML tree and build Geometry."""
        order = int(gtree.attrib["order"])
        metrics = []
        for i in range(order):
            metricname = gtree[i].tag
            if metricname not in _metrics:
                raise Exception("Invalid metricname {}".format(metricname))
            metrics.append(_metrics[metricname].load(gtree[i]))
        if gtree[-2].text: trasl = np.array(gtree[-2].text.split(), dtype="float64")
        else: trasl = None
        if gtree[-1].text: rot = np.array(gtree[-1].text.split(), dtype="float64")
        else: rot = None
        return Geometry(metrics, trasl=trasl, rot=rot)

# Inherited Metrics

class Energy (Metric):
    def __init__(self):
        """Simple metric for energy, with no transformation."""
        super().__init__([0], ["ekin"], ["MeV"], "MeV")
    def load(mtree):
        """Build Energy."""
        return Energy()

class Lethargy (Metric):
    def __init__(self, E0=10):
        """
        Lethargy metric for energy.
        
        Lethargy is defined as:
            u = log(E0 / ekin)

        Parameters
        ----------
        E0: float
            Reference energy. Typically, it is the highest energy in the
            system.
        """
        super().__init__([0], ["u"], ["[let]"], "MeV")
        self.E0 = E0
    def transform(self, ekins):
        """Transform energy to lethargy."""
        return np.log(self.E0/ekins)
    def inverse_transform(self, us):
        """Transform lethargy to energy."""
        return self.E0 * np.exp(-us)
    def jac(self, ekins):
        """Jacobian of lethargy transformation."""
        return 1/ekins.reshape(-1)
    def save(self, mtree):
        """Save Lethargy parameters into XML tree."""
        ET.SubElement(mtree, "dim").text = str(self.dim)
        paramsel = ET.SubElement(mtree, "params")
        paramsel.set("nps", "1")
        paramsel.text = "{}".format(self.E0)
    @staticmethod
    def load(mtree):
        """Load parameters from XML tree and build Lethargy."""
        dim = int(mtree[0].text)
        params = np.array(mtree[1].text.split(), dtype="float64")
        if dim!=1 or len(params)!=1 or int(mtree[1].attrib["nps"])!=1:
            raise Exception("Invalid metric tree.")
        return Lethargy(*params)

class Vol (Metric):
    def __init__(self, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf,
        zmin=-np.inf, zmax=np.inf):
        """
        Simple metric for 3D position, with no transformation.

        Each spatial variable (x, y, z) is delimited between a min and
        max value. By default these are -infinity and infinity,
        respectively. All positions in the particle list should be
        inside these limits.
        """
        super().__init__([1,2,3], ["x","y","z"], ["cm","cm","cm"], "cm^3")
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
    def save(self, mtree):
        """Save Vol parameters into XML tree."""
        ET.SubElement(mtree, "dim").text = str(self.dim)
        paramsel = ET.SubElement(mtree, "params")
        paramsel.set("nps", "6")
        paramsel.text = "{} {} {} {} {} {}".format(self.xmin, self.xmax,
            self.ymin, self.ymax, self.zmin, self.zmax)
    @staticmethod
    def load(mtree):
        """Load parameters from XML tree and build Vol."""
        dim = int(mtree[0].text)
        params = np.array(mtree[1].text.split(), dtype="float64")
        if dim!=3 or len(params)!=6 or int(mtree[1].attrib["nps"])!=6:
            raise Exception("Invalid metric tree.")
        return Vol(*params)

class SurfXY (Metric):
    def __init__(self, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, z=0):
        """
        Simple metric for 2D position, with no transformation.

        Spatial variables x and y are delimited between a min and max
        value. By default these are -infinity and infinity,
        respectively. All positions in the particle list should be
        inside these limits.

        z has the fixed value given as argument.
        """
        super().__init__([1,2,3], ["x","y"], ["cm","cm"], "cm^2")
        self.z = z
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
    def transform(self, poss):
        """Transform volume position (x,y,z) to flat position (x,y)."""
        return poss[:,:2]
    def inverse_transform(self, poss):
        """Transform flat position (x,y) to volume position (x,y,z)."""
        z_col = np.broadcast_to(self.z, (*poss.shape[:-1],1))
        return np.concatenate((poss, z_col), axis=1)
    def save(self, mtree):
        """Save SurfXY parameters into XML tree."""
        ET.SubElement(mtree, "dim").text = str(self.dim)
        paramsel = ET.SubElement(mtree, "params")
        paramsel.set("nps", "5")
        paramsel.text = "{} {} {} {} {}".format(self.xmin, self.xmax, self.ymin, self.ymax, self.z)
    @staticmethod
    def load(mtree):
        """Load parameters from XML tree and build SurfXY."""
        dim = int(mtree[0].text)
        params = np.array(mtree[1].text.split(), dtype="float64")
        if dim!=2 or len(params)!=5 or int(mtree[1].attrib["nps"])!=5:
            raise Exception("Invalid metric tree.")
        return SurfXY(*params)

class Guide (Metric):
    def __init__(self, xwidth, yheight, zmax=np.inf, rcurv=None):
        """
        Guide metric for position and direction.

        Position is parametrized with following variables:
            z: distance along guide, following curvature (if any).
            t: transversal direction along mirrors, starting at (x+,y-)
               corner, towards (x+,y+) corner.
        Direction is parametrized with following variables:
            mu: cosine of angle between particle direction and mirror
                normal.
            phi: azimuthal angle, starting from z direction, in [deg].

        Parameters
        ----------
        xwidth: float
            Guide width.
        yheight: float
            Guide height.
        zmax: float
            Guide length.
        rcurv: float
            Curvature radius, defined as follows. Default is no
            curvature.
                rcurv > 0 for curvature towards negative x
                rcurv < 0 for curvature towards negative x
                rcurv = 0 or rcurv = infinity for no curvature
        """
        super().__init__([1,2,3,4,5,6], ["z","t","theta","phi"], ["cm","cm","deg","deg"], "cm^2 sr")
        self.xwidth = xwidth
        self.yheight = yheight
        self.zmax = zmax
        self.rcurv = rcurv
    def transform(self, posdirs):
        """Transform position and direction to guide variables."""
        xs,ys,zs,dxs,dys,dzs = posdirs.T
        if self.rcurv is not None:
            rs = np.sqrt((self.rcurv+xs)**2 + zs**2)
            xs = np.sign(self.rcurv) * rs - self.rcurv; zs = np.abs(self.rcurv) * np.arcsin(zs / rs)
            dxs2 = dxs; dzs2 = dzs; angs = zs/self.rcurv
            dxs = dxs2*np.cos(angs) + dzs2*np.sin(angs); dzs = -dxs2*np.sin(angs) + dzs2*np.cos(angs)
        mask0 = np.logical_and((ys/self.yheight > -xs/self.xwidth), (ys/self.yheight <  xs/self.xwidth)) # mirror x pos
        mask1 = np.logical_and((ys/self.yheight >  xs/self.xwidth), (ys/self.yheight > -xs/self.xwidth)) # mirror y pos
        mask2 = np.logical_and((ys/self.yheight < -xs/self.xwidth), (ys/self.yheight >  xs/self.xwidth)) # mirror x neg
        mask3 = np.logical_and((ys/self.yheight <  xs/self.xwidth), (ys/self.yheight < -xs/self.xwidth)) # mirror y neg
        ts = np.empty_like(xs)
        ts[mask0] = 0.5*self.yheight +                   ys[mask0]
        ts[mask1] = 1.0*self.yheight + 0.5*self.xwidth - xs[mask1]
        ts[mask2] = 1.5*self.yheight + 1.0*self.xwidth - ys[mask2]
        ts[mask3] = 2.0*self.yheight + 1.5*self.xwidth + xs[mask3]
        mus = np.empty_like(dxs); phis = np.empty_like(dxs)
        mus[mask0] =  dxs[mask0]; phis[mask0] = np.arctan2(-dys[mask0], dzs[mask0])
        mus[mask1] =  dys[mask1]; phis[mask1] = np.arctan2( dxs[mask1], dzs[mask1])
        mus[mask2] = -dxs[mask2]; phis[mask2] = np.arctan2( dys[mask2], dzs[mask2])
        mus[mask3] = -dys[mask3]; phis[mask3] = np.arctan2(-dxs[mask3], dzs[mask3])
        phis *= 180/np.pi
        return np.stack((zs,ts,mus,phis), axis=1)
    def inverse_transform(self, posdirs):
        """Transform guide variables to position and direction."""
        zs,ts,mus,phis = posdirs.T
        phis *= np.pi/180
        mask0 =                                                   (ts <   self.yheight)              # mirror x pos
        mask1 = np.logical_and((ts >   self.yheight)            , (ts <   self.yheight+self.xwidth)) # mirror y pos
        mask2 = np.logical_and((ts >   self.yheight+self.xwidth), (ts < 2*self.yheight+self.xwidth)) # mirror x neg
        mask3 =                (ts > 2*self.yheight+self.xwidth)                                     # mirror y neg
        xs = np.empty_like(ts); ys = np.empty_like(ts)
        xs[mask0] =  self.xwidth /2; ys[mask0] =  ts[mask0] - 0.5*self.yheight
        ys[mask1] =  self.yheight/2; xs[mask1] = -ts[mask1] + 1.0*self.yheight + 0.5*self.xwidth
        xs[mask2] = -self.xwidth /2; ys[mask2] = -ts[mask2] + 1.5*self.yheight + 1.0*self.xwidth
        ys[mask3] = -self.yheight/2; xs[mask3] =  ts[mask3] - 2.0*self.yheight - 1.5*self.xwidth
        dxs = np.empty_like(mus); dys = np.empty_like(mus); dzs = np.empty_like(mus)
        dxs[mask0] =  mus[mask0]; dzs[mask0] = np.sqrt(1-mus[mask0]**2)*np.cos(phis[mask0]); dys[mask0] = -np.sqrt(1-mus[mask0]**2)*np.sin(phis[mask0])
        dys[mask1] =  mus[mask1]; dzs[mask1] = np.sqrt(1-mus[mask1]**2)*np.cos(phis[mask1]); dxs[mask1] =  np.sqrt(1-mus[mask1]**2)*np.sin(phis[mask1])
        dxs[mask2] = -mus[mask2]; dzs[mask2] = np.sqrt(1-mus[mask2]**2)*np.cos(phis[mask2]); dys[mask2] =  np.sqrt(1-mus[mask2]**2)*np.sin(phis[mask2])
        dys[mask3] = -mus[mask3]; dzs[mask3] = np.sqrt(1-mus[mask3]**2)*np.cos(phis[mask3]); dxs[mask3] = -np.sqrt(1-mus[mask3]**2)*np.sin(phis[mask3])
        if self.rcurv is not None:
            rs = (self.rcurv + xs) * np.sign(self.rcurv); angs = zs/self.rcurv
            xs = np.sign(self.rcurv) * rs * np.cos(zs/self.rcurv) - self.rcurv; zs = rs * np.sin(np.abs(zs/self.rcurv))
            dxs2 = dxs; dzs2 = dzs
            dxs = dxs2*np.cos(angs) - dzs2*np.sin(angs); dzs = dxs2*np.sin(angs) + dzs2*np.cos(angs)
        return np.stack((xs,ys,zs,dxs,dys,dzs), axis=1)
    def save(self, mtree):
        """Save Guide parameters into XML tree."""
        ET.SubElement(mtree, "dim").text = str(self.dim)
        paramsel = ET.SubElement(mtree, "params")
        paramsel.set("nps", "4")
        paramsel.text = "{} {} {} {}".format(self.xwidth, self.yheight, self.zmax, self.rcurv)
    @staticmethod
    def load(mtree):
        """Load parameters from XML tree and build Guide."""
        dim = int(mtree[0].text)
        params = np.array(mtree[1].text.split(), dtype="float64")
        if dim!=6 or len(params)!=4 or int(mtree[1].attrib["nps"])!=4:
            raise Exception("Invalid metric tree.")
        return Guide(*params)

class Isotrop (Metric):
    def __init__(self, keep_xdir=False, keep_ydir=False, keep_zdir=False):
        """
        Simple metric for direction, with no transformation.

        Distance is measured as the euclidean distance between 3D
        unitary direction vectors.

        Parameters
        ----------
        keep_xdir: bool
            If True, when using the source for sampling new particles,
            perturbation will not change dirx sign.
        keep_ydir: bool
            If True, when using the source for sampling new particles,
            perturbation will not change diry sign.
        keep_zdir: bool
            If True, when using the source for sampling new particles,
            perturbation will not change dirz sign.
        """
        super().__init__([4,5,6], ["dx","dy","dz"], ["[dir]","[dir]","[dir]"], "[dir]^3")
        self.keep_xdir = keep_xdir
        self.keep_ydir = keep_ydir
        self.keep_zdir = keep_zdir
    def mean(self, dirs=None, vecs=None, weights=None):
        """
        Mean of directions.

        Mean is computed as the euclidean mean of direction vectors,
        normalized to 1.

        Parameters
        ----------
        dirs: array-like, optional
            Array of directions.
        vecs: array-like, optional
            Array of parametrized directions. If set, overrides dirs.
        weights: array-like, optional
            Array of particle statistic weights.
        """
        if vecs is None:
            vecs = self.transform(dirs)
        mn = np.average(vecs, axis=0, weights=weights)
        mn_norm = np.linalg.norm(mn)
        if mn_norm == 0: mn_norm = 1
        return mn / mn_norm
    def std(self, dirs=None, vecs=None, weights=None):
        """
        Standard deviation of directions.

        Standard deviation is computed as the euclidean standard
        deviation of direction vectors minus its mean (computed with
        mean method).

        Parameters
        ----------
        dirs: array-like, optional
            Array of directions.
        vecs: array-like, optional
            Array of parametrized directions. If set, overrides dirs.
        weights: array-like, optional
            Array of particle statistic weights.
        """
        if vecs is None:
            vecs = self.transform(dirs)
        mn = self.mean(vecs=vecs, weights=weights)
        std = np.sqrt(np.mean(np.average((vecs-mn)**2, axis=0, weights=weights)))
        return np.array(3*[std])
    def save(self, mtree):
        """Save Guide parameters into XML tree."""
        ET.SubElement(mtree, "dim").text = str(self.dim)
        paramsel = ET.SubElement(mtree, "params")
        paramsel.set("nps", "3")
        paramsel.text = "{:d} {:d} {:d}".format(self.keep_xdir, self.keep_ydir, self.keep_zdir)
    @staticmethod
    def load(mtree):
        """Load parameters from XML tree and build Isotrop."""
        dim = int(mtree[0].text)
        params = np.array(mtree[1].text.split(), dtype="float64")
        if dim!=3 or len(params)!=3 or int(mtree[1].attrib["nps"])!=3:
            raise Exception("Invalid metric tree.")
        return Isotrop(*params)

class Polar (Metric):
    def __init__(self):
        """
        Polar parametrization for direction.

        Polar angles are defined as follows:
            theta: angle between direction and z, in [deg].
            phi: azimuthal angle, starting from x direction, in [deg].
        """
        super().__init__([4,5,6], ["theta","phi"], ["deg","deg"], "sr^2")
    def transform(self, dirs):
        """Transform directions to polar angles."""
        thetas = np.arccos(dirs[:,2]) * 180/np.pi
        phis = np.arctan2(dirs[:,1], dirs[:,0]) * 180/np.pi
        return np.stack((thetas, phis), axis=1)
    def inverse_transform(self, tps):
        """Transform polar angles to directions."""
        thetas,phis = tps.T
        dxs = np.sin(thetas*np.pi/180) * np.cos(phis*np.pi/180)
        dys = np.sin(thetas*np.pi/180) * np.sin(phis*np.pi/180)
        dzs = np.cos(thetas*np.pi/180)
        return np.stack((dxs, dys, dzs), axis=1)
    def jac(self, dirs):
        """Jacobian of polar transformation."""
        thetas = np.arccos(dirs[:,2])
        return (180/np.pi)**2 / np.sin(thetas)
    @staticmethod
    def load(mtree):
        """Build Polar."""
        return Polar()

class PolarMu (Metric):
    def __init__(self):
        """
        Polar parametrization for direction, with mu = cos(theta).

        Polar parameters are defined as follows:
            mu: cosine of angle between direction and z, in [deg].
            phi: azimuthal angle, starting from x direction, in [deg].
        """
        super().__init__([4,5,6], ["mu","phi"], ["","deg"], "sr")
    def transform(self, dirs):
        """Transform directions to polar parameters."""
        mus = dirs[:,2]
        phis = np.arctan2(dirs[:,1], dirs[:,0]) * 180/np.pi
        return np.stack((mus, phis), axis=1)
    def inverse_transform(self, tps):
        """Transform polar parameters to directions."""
        mus,phis = tps.T
        dxs = np.sqrt(1-mus**2) * np.cos(phis*np.pi/180)
        dys = np.sqrt(1-mus**2) * np.sin(phis*np.pi/180)
        dzs = mus
        return np.stack((dxs, dys, dzs), axis=1)
    @staticmethod
    def load(mtree):
        """Build PolarMu."""
        return PolarMu()

_metrics = {
    "Energy":Energy,
    "Lethargy":Lethargy,
    "Vol":Vol,
    "SurfXY":SurfXY,
    "Guide":Guide,
    "Isotrop":Isotrop,
    "Polar":Polar,
    "PolarMu":PolarMu
}

# Aliases for usual geometries

def GeomFlat(xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, z=0, E0=10, keep_zdir=True, trasl=None, rot=None):
    """
    Build flat neutron source.

    Energy metric is Lethargy, position metric is SurfXY, and direction 
    metric is Isotrop.

    See Metric's and Geometry constructors for parameters docs.
    """
    return Geometry([Lethargy(E0), SurfXY(xmin,xmax,ymin,ymax,z), Isotrop(keep_zdir=keep_zdir)], trasl=trasl, rot=rot)

def GeomGuide(xwidth, yheight, zmax=np.inf, rcurv=None, E0=10, trasl=None, rot=None):
    """
    Build neutron source for leaks thru guide mirrors.

    Energy metric is Lethargy, position and direction metric is Guide.

    See Metric's and Geometry constructors for parameters docs.
    """
    return Geometry([Lethargy(E0), Guide(xwidth,yheight,zmax,rcurv)], trasl=trasl, rot=rot)

def GeomActiv(xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=-np.inf, zmax=np.inf, trasl=None, rot=None):
    """
    Build photon volumetric activation source.

    Energy metric is Energy, position metric is Vol, and direction
    metric is Isotrop.

    See Metric's and Geometry constructors for parameters docs.
    """
    return Geometry([Energy(), Vol(xmin,xmax,ymin,ymax,zmin,zmax), Isotrop()], trasl=trasl, rot=rot)
