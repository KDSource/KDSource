# -*- coding: utf-8 -*-

from xml.etree import ElementTree as ET

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class Metric:
	def __init__(self, partvars, varnames, units, volunits):
		self.partvars = partvars
		if len(varnames) != len(units):
			raise ValueError("Longitudes de varnames y units deben ser iguales")
		self.dim = len(varnames)
		self.varnames = varnames
		self.varmap = {name:idx for idx,name in enumerate(varnames)}
		self.units = units
		self.volunits = volunits
	def transform(self, parts):
		return parts
	def inverse_transform(self, vecs):
		return vecs
	def jac(self, parts):
		return np.ones(len(parts))
	def mean(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		return self.inverse_transform(np.mean(vecs, axis=0))
	def std(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		mn = np.mean(vecs, axis=0)
		return np.std(vecs-mn, axis=0)
	def save(self, mtree):
		ET.SubElement(mtree, "dim").text = str(self.dim)
		ET.SubElement(mtree, "params").set("nps", "0")
	@staticmethod
	def load(mtree):
		raise Exception("Load method not implemented")

class Geometry (Metric):
	def __init__(self, metrics, trasl=None, rot=None):
		partvars = range(7)
		varnames = sum([metric.varnames for metric in metrics], [])
		units = sum([metric.units for metric in metrics], [])
		volunits = "".join([metric.volunits+" " for metric in metrics])[:-1]
		super().__init__(partvars, varnames, units, volunits)
		self.ms = metrics
		if trasl is not None:
			trasl = np.array(trasl).reshape(-1)
			if trasl.shape != (3,):
				raise ValueError("trasl invalido")
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
					raise ValueError("rot invalido")
		self.trasl = trasl
		self.rot = rot
	def transform(self, parts):
		if self.trasl is not None:
			parts[:,1:4] -= self.trasl # Posicion
		if self.rot is not None:
			parts[:,1:4] = self.rot.apply(parts[:,1:4], inverse=True) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7], inverse=True) # Direccion
		vecss = []
		for metric in self.ms:
			vecss.append(metric.transform(parts[:,metric.partvars]))
		return np.concatenate(vecss, axis=1)
	def inverse_transform(self, vecs):
		parts = np.zeros((len(vecs), 7))
		end = 0
		for metric in self.ms:
			start = end
			end = start + metric.dim
			parts[:,metric.partvars] = metric.inverse_transform(vecs[:,start:end])
		if self.trasl is not None:
			parts[:,1:4] += self.trasl # Posicion
		if self.rot is not None:
			parts[:,1:4] = self.rot.apply(parts[:,1:4]) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7]) # Direccion
		return parts
	def jac(self, parts):
		if self.trasl is not None:
			parts[:,1:4] -= self.trasl # Posicion
		if self.rot is not None:
			parts[:,1:4] = self.rot.apply(parts[:,1:4], inverse=True) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7], inverse=True) # Direccion
		jacs = []
		for metric in self.ms:
			jacs.append(metric.jac(parts[:,metric.partvars]))
		return np.prod(jacs, axis=1)
	def mean(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		means = []
		end = 0
		for metric in self.ms:
			start = end
			end = start + metric.dim
			means.append(metric.mean(vecs=vecs[:,start:end]))
		return np.concatenate(means)
	def std(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		stds = []
		end = 0
		for metric in self.ms:
			start = end
			end = start + metric.dim
			stds.append(metric.std(vecs=vecs[:,start:end]))
		return np.concatenate(stds)
	def save(self, gtree):
		gtree.set("order", str(len(self.ms)))
		for metric in self.ms:
			mtree = ET.SubElement(gtree, metric.__class__.__name__)
			metric.save(mtree)
		trasl = np.array_str(self.trasl)[2:-2] if self.trasl is not None else ""
		ET.SubElement(gtree, "trasl").text = trasl
		rot = np.array_str(self.rot.as_rotvec())[2:-2] if self.rot is not None else ""
		ET.SubElement(gtree, "rot").text = rot
	@staticmethod
	def load(gtree):
		order = int(gtree.attrib["order"])
		metrics = []
		for i in range(order):
			metricname = gtree[i].tag
			if metricname not in _metrics:
				raise Exception("Invalid metric name {}".format(metricname))
			metrics.append(_metrics[metricname].load(gtree[i]))
		if gtree[-2].text: trasl = np.array(gtree[-2].text.split(), dtype="float64")
		else: trasl = None
		if gtree[-1].text: rot = np.array(gtree[-1].text.split(), dtype="float64")
		else: rot = None
		return Geometry(metrics, trasl=trasl, rot=rot)

# Clases heredadas

class Energy (Metric):
	def __init__(self):
		super().__init__([0], ["E"], ["MeV"], "MeV")
	def load(mtree):
		return Energy()

class Lethargy (Metric):
	def __init__(self, E0=10):
		super().__init__([0], ["u"], ["[let]"], "MeV")
		self.E0 = E0
	def transform(self, Es):
		return np.log(self.E0/Es)
	def inverse_transform(self, us):
		return self.E0 * np.exp(-us)
	def jac(self, Es):
		return 1/Es.reshape(-1)
	def save(self, mtree):
		ET.SubElement(mtree, "dim").text = str(self.dim)
		paramsel = ET.SubElement(mtree, "params")
		paramsel.set("nps", "1")
		paramsel.text = "{}".format(self.E0)
	@staticmethod
	def load(mtree):
		dim = int(mtree[0].text)
		params = np.array(mtree[1].text.split(), dtype="float64")
		if dim!=1 or len(params)!=1 or int(mtree[1].attrib["nps"])!=1:
			raise Exception("Invalid metric tree")
		return Lethargy(*params)

class Vol (Metric):
	def __init__(self, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=-np.inf, zmax=np.inf):
		super().__init__([1,2,3], ["x","y","z"], ["cm","cm","cm"], "cm^3")
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
		self.zmin = zmin
		self.zmax = zmax
	def save(self, mtree):
		ET.SubElement(mtree, "dim").text = str(self.dim)
		paramsel = ET.SubElement(mtree, "params")
		paramsel.set("nps", "6")
		paramsel.text = "{} {} {} {} {} {}".format(self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax)
	@staticmethod
	def load(mtree):
		dim = int(mtree[0].text)
		params = np.array(mtree[1].text.split(), dtype="float64")
		if dim!=3 or len(params)!=6 or int(mtree[1].attrib["nps"])!=6:
			raise Exception("Invalid metric tree")
		return Vol(*params)

class SurfXY (Metric):
	def __init__(self, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, z=0):
		super().__init__([1,2,3], ["x","y"], ["cm","cm"], "cm^2")
		self.z = z
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
	def transform(self, poss):
		return poss[:,:2]
	def inverse_transform(self, poss):
		z_col = np.broadcast_to(self.z, (*poss.shape[:-1],1))
		return np.concatenate((poss, z_col), axis=1)
	def save(self, mtree):
		ET.SubElement(mtree, "dim").text = str(self.dim)
		paramsel = ET.SubElement(mtree, "params")
		paramsel.set("nps", "5")
		paramsel.text = "{} {} {} {} {}".format(self.xmin, self.xmax, self.ymin, self.ymax, self.z)
	@staticmethod
	def load(mtree):
		dim = int(mtree[0].text)
		params = np.array(mtree[1].text.split(), dtype="float64")
		if dim!=2 or len(params)!=5 or int(mtree[1].attrib["nps"])!=5:
			raise Exception("Invalid metric tree")
		return SurfXY(*params)

class Guide (Metric):
	def __init__(self, xwidth, yheight, zmax=np.inf, rcurv=None):
		super().__init__([1,2,3,4,5,6], ["z","t","theta","phi"], ["cm","cm","deg","deg"], "cm^2 sr")
		self.xwidth = xwidth
		self.yheight = yheight
		self.zmax = zmax
		self.rcurv = rcurv
	def transform(self, posdirs):
		xs,ys,zs,dxs,dys,dzs = posdirs.T
		if self.rcurv is not None:
			rs = np.sqrt((self.rcurv+xs)**2 + zs**2)
			xs = np.sign(self.rcurv) * rs - self.rcurv; zs = np.abs(self.rcurv) * np.arcsin(zs / rs)
			dxs2 = dxs; dzs2 = dzs; angs = zs/self.rcurv
			dxs = dxs2*np.cos(angs) + dzs2*np.sin(angs); dzs = -dxs2*np.sin(angs) + dzs2*np.cos(angs)
		mask0 = np.logical_and((ys/self.yheight > -xs/self.xwidth), (ys/self.yheight <  xs/self.xwidth)) # espejo x pos
		mask1 = np.logical_and((ys/self.yheight >  xs/self.xwidth), (ys/self.yheight > -xs/self.xwidth)) # espejo y pos
		mask2 = np.logical_and((ys/self.yheight < -xs/self.xwidth), (ys/self.yheight >  xs/self.xwidth)) # espejo x neg
		mask3 = np.logical_and((ys/self.yheight <  xs/self.xwidth), (ys/self.yheight < -xs/self.xwidth)) # espejo y neg
		ts = np.zeros_like(xs)
		ts[mask0] = 0.5*self.yheight +                   ys[mask0]
		ts[mask1] = 1.0*self.yheight + 0.5*self.xwidth - xs[mask1]
		ts[mask2] = 1.5*self.yheight + 1.0*self.xwidth - ys[mask2]
		ts[mask3] = 2.0*self.yheight + 1.5*self.xwidth + xs[mask3]
		thetas = np.zeros_like(dxs); phis = np.zeros_like(dxs)
		thetas[mask0] = np.arccos( dxs[mask0]); phis[mask0] = np.arctan2(-dys[mask0], dzs[mask0])
		thetas[mask1] = np.arccos( dys[mask1]); phis[mask1] = np.arctan2( dxs[mask1], dzs[mask1])
		thetas[mask2] = np.arccos(-dxs[mask2]); phis[mask2] = np.arctan2( dys[mask2], dzs[mask2])
		thetas[mask3] = np.arccos(-dys[mask3]); phis[mask3] = np.arctan2(-dxs[mask3], dzs[mask3])
		thetas *= 180/np.pi; phis *= 180/np.pi
		return np.stack((zs,ts,thetas,phis), axis=1)
	def inverse_transform(self, posdirs):
		zs,ts,thetas,phis = posdirs.T
		thetas *= np.pi/180; phis *= np.pi/180
		mask0 =                                                   (ts <   self.yheight)              # espejo x pos
		mask1 = np.logical_and((ts >   self.yheight)            , (ts <   self.yheight+self.xwidth)) # espejo y pos
		mask2 = np.logical_and((ts >   self.yheight+self.xwidth), (ts < 2*self.yheight+self.xwidth)) # espejo x neg
		mask3 =                (ts > 2*self.yheight+self.xwidth)                                     # espejo y neg
		xs = np.zeros_like(ts); ys = np.zeros_like(ts)
		xs[mask0] =  self.xwidth /2; ys[mask0] =  ts[mask0] - 0.5*self.yheight
		ys[mask1] =  self.yheight/2; xs[mask1] = -ts[mask1] + 1.0*self.yheight + 0.5*self.xwidth
		xs[mask2] = -self.xwidth /2; ys[mask2] = -ts[mask2] + 1.5*self.yheight + 1.0*self.xwidth
		ys[mask3] = -self.yheight/2; xs[mask3] =  ts[mask3] - 2.0*self.yheight - 1.5*self.xwidth
		dxs = np.zeros_like(thetas); dys = np.zeros_like(thetas); dzs = np.zeros_like(thetas)
		dxs[mask0] =  np.cos(thetas[mask0]); dzs[mask0] = np.sin(thetas[mask0])*np.cos(phis[mask0]); dys[mask0] = -np.sin(thetas[mask0])*np.sin(phis[mask0])
		dys[mask1] =  np.cos(thetas[mask1]); dzs[mask1] = np.sin(thetas[mask1])*np.cos(phis[mask1]); dxs[mask1] =  np.sin(thetas[mask1])*np.sin(phis[mask1])
		dxs[mask2] = -np.cos(thetas[mask2]); dzs[mask2] = np.sin(thetas[mask2])*np.cos(phis[mask2]); dys[mask2] =  np.sin(thetas[mask2])*np.sin(phis[mask2])
		dys[mask3] = -np.cos(thetas[mask3]); dzs[mask3] = np.sin(thetas[mask3])*np.cos(phis[mask3]); dxs[mask3] = -np.sin(thetas[mask3])*np.sin(phis[mask3])
		if self.rcurv is not None:
			rs = (self.rcurv + xs) * np.sign(self.rcurv); angs = zs/self.rcurv
			xs = np.sign(self.rcurv) * rs * np.cos(zs/self.rcurv) - self.rcurv; zs = rs * np.sin(np.abs(zs/self.rcurv))
			dxs2 = dxs; dzs2 = dzs
			dxs = dxs2*np.cos(angs) - dzs2*np.sin(angs); dzs = dxs2*np.sin(angs) + dzs2*np.cos(angs)
		return np.stack((xs,ys,zs,dxs,dys,dzs), axis=1)
	def jac(self, posdirs):
		xs,ys,zs,dxs,dys,dzs = posdirs.T
		if self.rcurv is not None:
			rs = np.sqrt((self.rcurv+xs)**2 + zs**2)
			xs = np.sign(self.rcurv) * rs - self.rcurv; zs = np.abs(self.rcurv) * np.arcsin(zs / rs)
			dxs2 = dxs; dzs2 = dzs; angs = zs/self.rcurv
			dxs = dxs2*np.cos(angs) + dzs2*np.sin(angs); dzs = -dxs2*np.sin(angs) + dzs2*np.cos(angs)
		mask0 = np.logical_and((ys/self.yheight > -xs/self.xwidth), (ys/self.yheight <  xs/self.xwidth)) # espejo x pos
		mask1 = np.logical_and((ys/self.yheight >  xs/self.xwidth), (ys/self.yheight > -xs/self.xwidth)) # espejo y pos
		mask2 = np.logical_and((ys/self.yheight < -xs/self.xwidth), (ys/self.yheight >  xs/self.xwidth)) # espejo x neg
		mask3 = np.logical_and((ys/self.yheight <  xs/self.xwidth), (ys/self.yheight < -xs/self.xwidth)) # espejo y neg
		thetas = np.zeros_like(dxs)
		thetas[mask0] = np.arccos( dxs[mask0])
		thetas[mask1] = np.arccos( dys[mask1])
		thetas[mask2] = np.arccos(-dxs[mask2])
		thetas[mask3] = np.arccos(-dys[mask3])
		return (180/np.pi)**2 / np.sin(thetas)
	def save(self, mtree):
		ET.SubElement(mtree, "dim").text = str(self.dim)
		paramsel = ET.SubElement(mtree, "params")
		paramsel.set("nps", "4")
		paramsel.text = "{} {} {} {}".format(self.xwidth, self.yheight, self.zmax, rcurv)
	@staticmethod
	def load(mtree):
		dim = int(mtree[0].text)
		params = np.array(mtree[1].text.split(), dtype="float64")
		if dim!=6 or len(params)!=4 or int(mtree[1].attrib["nps"])!=4:
			raise Exception("Invalid metric tree")
		return Guide(*params)

class Isotrop (Metric):
	def __init__(self):
		super().__init__([4,5,6], ["dx","dy","dz"], ["[dir]","[dir]","[dir]"], "[dir]^3")
	def mean(self, dirs, transform=False):
		if transform:
			vecs = self.transform(dirs)
		else:
			vecs = dirs
		mn = np.mean(dirs, axis=0)
		return mn / np.linalg.norm(mn)
	def std(self, dirs=None, vecs=None):
		if vecs is None:
			vecs = self.transform(dirs)
		mn = self.mean(vecs, transform=False)
		std = np.std(vecs-mn)
		return np.array(3*[std])
	@staticmethod
	def load(mtree):
		return Isotrop()

class Polar (Metric):
	def __init__(self):
		super().__init__([4,5,6], ["theta","phi"], ["deg","deg"], "sr^2")
	def transform(self, dirs):
		thetas = np.arccos(dirs[:,2]) * 180/np.pi
		phis = np.arctan2(dirs[:,1], dirs[:,0]) * 180/np.pi
		return np.stack((thetas, phis), axis=1)
	def inverse_transform(self, tps):
		thetas,phis = tps.T
		dxs = np.sin(thetas*np.pi/180) * np.cos(phis*np.pi/180)
		dys = np.sin(thetas*np.pi/180) * np.sin(phis*np.pi/180)
		dzs = np.cos(thetas*np.pi/180)
		return np.stack((dxs, dys, dzs), axis=1)
	def jac(self, dirs):
		thetas = np.arccos(dirs[:,2])
		return (180/np.pi)**2 / np.sin(thetas)
	@staticmethod
	def load(mtree):
		return Polar()

_metrics = {
	"Energy":Energy,
	"Lethargy":Lethargy,
	"Vol":Vol,
	"SurfXY":SurfXY,
	"Guide":Guide,
	"Isotrop":Isotrop,
	"Polar":Polar
}

# Alias para geometrias mas usuales

def GeomFlat(xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, z=0, E0=10, trasl=None, rot=None):
	return Geometry([Lethargy(E0), SurfXY(xmin,xmax,ymin,ymax,z), Polar()], trasl=trasl, rot=rot)

def GeomGuide(xwidth, yheight, zmax=np.inf, rcurv=None, E0=10, trasl=None, rot=None):
	return Geometry([Lethargy(E0), Guide(xwidth,yheight,zmax,rcurv)], trasl=trasl, rot=rot)

def GeomActiv(xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=-np.inf, zmax=np.inf, E0=10, trasl=None, rot=None):
	return Geometry([Lethargy(E0), Vol(xmin,xmax,ymin,ymax,zmin,zmax), Polar()], trasl=trasl, rot=rot)