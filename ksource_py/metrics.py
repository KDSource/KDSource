# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class Metric:
	def __init__(self, partdim, varnames, units, volunits):
		self.partdim = partdim
		self.varnames = varnames
		self.varmap = {name:idx for idx,name in enumerate(varnames)}
		self.units = units
		self.volunits = volunits
		self.dim = dim = len(varnames)
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
	def save(self, file):
		file.write(self.__class__.__name__+'\n')
		file.write("{}\n".format(self.dim))
		file.write('0\n')

class Geometry (Metric):
	def __init__(self, metrics, trasl=None, rot=None):
		partdim = sum([metric.partdim for metric in metrics])
		varnames = sum([metric.varnames for metric in metrics], [])
		units = sum([metric.units for metric in metrics], [])
		volunits = "".join([metric.volunits+" " for metric in metrics])[:-1]
		super().__init__(partdim, varnames, units, volunits)
		self.ms = metrics
		if trasl is not None:
			trasl = np.array(trasl)
		if rot is not None:
			rot = st.Rotation.from_rotvec(rot)
		self.trasl = trasl
		self.rot = rot
	def transform(self, parts):
		if self.trasl is not None:
			parts[:,1:4] -= self.trasl
		if self.rot is not None:
			parts[:,1:4] = self.rot.apply(parts[:,1:4], inverse=True) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7], inverse=True) # Direccion
		vecss = []
		end = 0
		for metric in self.ms:
			start = end
			end = start + metric.partdim
			vecss.append(metric.transform(parts[:,start:end]))
		return np.concatenate(vecss, axis=1)
	def inverse_transform(self, vecs):
		partss = []
		end = 0
		for metric in self.ms:
			start = end
			end = start + metric.dim
			partss.append(metric.inverse_transform(vecs[:,start:end]))
		parts = np.concatenate(partss, axis=1)
		if self.trasl is not None:
			parts[:,1:4] += self.trasl
		if self.rot is not None:
			parts[:,1:4] = self.rot.apply(parts[:,1:4]) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7]) # Direccion
		return parts
	def jac(self, parts):
		if self.trasl is not None:
			parts[:,1:4] -= self.trasl
		if self.rot is not None:
			parts[:,1:4] = self.rot.apply(parts[:,1:4], inverse=True) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7], inverse=True) # Direccion
		jacs = []
		end = 0
		for metric in self.ms:
			start = end
			end = start + metric.partdim
			jacs.append(metric.jac(parts[:,start:end]))
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
	def save(self, file):
		file.write("{}\n".format(len(self.ms)))
		for metric in self.ms:
			metric.save(file)
		if self.trasl is not None: np.savetxt(file, self.trasl[np.newaxis,:])
		else: file.write('\n')
		if self.rot is not None: np.savetxt(file, self.rot[np.newaxis,:])
		else: file.write('\n')

class Energy (Metric):
	def __init__(self):
		super().__init__(1, ["E"], ["MeV"], "MeV")

class Lethargy (Metric):
	def __init__(self, E0):
		super().__init__(1, ["u"], ["[let]"], "MeV")
		self.E0 = E0
	def transform(self, Es):
		return np.log(self.E0/Es)
	def inverse_transform(self, us):
		return self.E0 * np.exp(-us)
	def jac(self, Es):
		return 1/Es.reshape(-1)
	def save(self, file):
		file.write(self.__class__.__name__+'\n')
		file.write("{}\n".format(self.dim))
		file.write("1 {}\n".format(self.E0))

class Vol (Metric):
	def __init__(self, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, zmin=-np.inf, zmax=np.inf):
		super().__init__(3, ["x","y","z"], ["cm","cm","cm"], "cm^3")
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
		self.zmin = zmin
		self.zmax = zmax
	def save(self, file):
		file.write(self.__class__.__name__+'\n')
		file.write("{}\n".format(self.dim))
		file.write("6 {} {} {} {} {} {}\n".format(self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax))

class SurfXY (Metric):
	def __init__(self, z, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf):
		super().__init__(3, ["x","y"], ["cm","cm"], "cm^2")
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
	def save(self, file):
		file.write(self.__class__.__name__+'\n')
		file.write("{}\n".format(self.dim))
		file.write("4 {} {} {} {}\n".format(self.xmin, self.xmax, self.ymin, self.ymax))

class Guide (Metric):
	def __init__(self, xwidth, yheight, zmax=np.inf, rcurv=None):
		super().__init__(6, ["z","t","theta","phi"], ["cm","cm","deg","deg"], "cm^2 sr")
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
	def save(self, file):
		file.write(self.__class__.__name__+'\n')
		file.write("{}\n".format(self.dim))
		rcurv = self.rcurv
		if rcurv is None: rcurv = 0
		file.write("4 {} {} {} {}\n".format(self.xwidth, self.yheight, self.zmax, rcurv))

class Isotrop (Metric):
	def __init__(self):
		super().__init__(3, ["dx","dy","dz"], ["[dir]","[dir]","[dir]"], "[dir]^3")
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

class Polar (Metric):
	def __init__(self):
		super().__init__(3, ["mu","phi"], ["deg","deg"], "sr^2")
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