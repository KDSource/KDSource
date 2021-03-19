# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class Metric:
	def __init__(self, varnames, units, volunits):
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

class SepVarMetric (Metric):
	def __init__(self, metric_E, metric_pos, metric_dir, trasl=None, rot=None):
		varnames = metric_E.varnames+metric_pos.varnames+metric_dir.varnames
		units = metric_E.units+metric_pos.units+metric_dir.units
		volunits = metric_E.volunits+" "+metric_pos.volunits+" "+metric_dir.volunits
		super().__init__(varnames, units, volunits)
		self.E = metric_E
		self.pos = metric_pos
		self.dir = metric_dir
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
		transf_Es = self.E.transform(parts[:,0:1])
		transf_poss = self.pos.transform(parts[:,1:4])
		transf_dirs = self.dir.transform(parts[:,4:7])
		return np.concatenate((transf_Es, transf_poss, transf_dirs), axis=1)
	def inverse_transform(self, vecs):
		Es = self.E.inverse_transform(vecs[:,0:self.E.dim])
		poss = self.pos.inverse_transform(vecs[:,self.E.dim:self.E.dim+self.pos.dim])
		dirs = self.dir.inverse_transform(vecs[:,self.E.dim+self.pos.dim:])
		if self.trasl is not None:
			poss += self.trasl
		if self.rot is not None:
			poss = self.rot.apply(poss) # Posicion
			dirs = self.rot.apply(dirs) # Direccion
		return np.concatenate((Es, poss, dirs), axis=1)
	def jac(self, parts):
		if self.trasl is not None:
			parts[:,1:4] -= self.trasl
		if self.rot is not None:
			parts[:,1:4] = self.rot.apply(parts[:,1:4], inverse=True) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7], inverse=True) # Direccion
		jac_E = self.E.jac(parts[:,0:1])
		jac_pos = self.pos.jac(parts[:,1:4])
		jac_dir = self.dir.jac(parts[:,4:7])
		return jac_E*jac_pos*jac_dir
	def mean(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		mean_Es = self.E.mean(vecs=vecs[:,0:self.E.dim])
		mean_poss = self.pos.mean(vecs=vecs[:,self.E.dim:self.E.dim+self.pos.dim])
		mean_dirs = self.dir.mean(vecs=vecs[:,self.E.dim+self.pos.dim:])
		return np.concatenate((mean_Es, mean_poss, mean_dirs))
	def std(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		std_Es = self.E.std(vecs=vecs[:,0:self.E.dim])
		std_poss = self.pos.std(vecs=vecs[:,self.E.dim:self.E.dim+self.pos.dim])
		std_dirs = self.dir.std(vecs=vecs[:,self.E.dim+self.pos.dim:])
		return np.concatenate((std_Es, std_poss, std_dirs))
	def save(self, file):
		file.write("# E:\n")
		self.E.save(file)
		file.write("# pos:\n")
		self.pos.save(file)
		file.write("# dir:\n")
		self.dir.save(file)
		if self.trasl is not None: np.savetxt(file, self.trasl)
		else: file.write('\n')
		if self.rot is not None: np.savetxt(file, self.rot)
		else: file.write('\n')

class Energy (Metric):
	def __init__(self):
		super().__init__(["E"], ["MeV"], "MeV")

class Lethargy (Metric):
	def __init__(self, E0):
		super().__init__(["u"], ["[let]"], "MeV")
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
		super().__init__(["x","y","z"], ["cm","cm","cm"], "cm^3")
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
		super().__init__(["x","y"], ["cm","cm"], "cm^2")
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
		super().__init__(["z","t"], ["cm","cm"], "cm^2")
		self.xwidth = xwidth
		self.yheight = yheight
		self.zmax = zmax
		self.rcurv = rcurv
	def transform(self, poss):
		xs,ys,zs = poss.T
		if self.rcurv is not None:
			rs = np.sqrt((self.rcurv+xs)**2 + zs**2)
			xs = np.sign(self.rcurv) * rs - self.rcurv
			zs = np.abs(self.rcurv) * np.arcsin(zs / rs)
		mask0 = np.logical_and((ys/self.yheight > -xs/self.xwidth), (ys/self.yheight <  xs/self.xwidth)) # espejo x pos
		mask1 = np.logical_and((ys/self.yheight >  xs/self.xwidth), (ys/self.yheight > -xs/self.xwidth)) # espejo y pos
		mask2 = np.logical_and((ys/self.yheight < -xs/self.xwidth), (ys/self.yheight >  xs/self.xwidth)) # espejo x neg
		mask3 = np.logical_and((ys/self.yheight <  xs/self.xwidth), (ys/self.yheight < -xs/self.xwidth)) # espejo y neg
		ts = np.zeros_like(xs)
		ts[mask0] = 0.5*self.yheight +                   ys[mask0]
		ts[mask1] = 1.0*self.yheight + 0.5*self.xwidth - xs[mask1]
		ts[mask2] = 1.5*self.yheight + 1.0*self.xwidth - ys[mask2]
		ts[mask3] = 2.0*self.yheight + 1.5*self.xwidth + xs[mask3]
		return np.stack((zs,ts), axis=1)
	def inverse_transform(self, poss):
		zs,ts = poss.T
		mask0 =                                                   (ts <   self.yheight)              # espejo x pos
		mask1 = np.logical_and((ts >   self.yheight)            , (ts <   self.yheight+self.xwidth)) # espejo y pos
		mask2 = np.logical_and((ts >   self.yheight+self.xwidth), (ts < 2*self.yheight+self.xwidth)) # espejo x neg
		mask3 =                (ts > 2*self.yheight+self.xwidth)                                     # espejo y neg
		xs = np.zeros_like(ts)
		ys = np.zeros_like(ts)
		xs[mask0] =  self.xwidth /2; ys[mask0] =  ts[mask0] - 0.5*self.yheight
		ys[mask1] =  self.yheight/2; xs[mask1] = -ts[mask1] + 1.0*self.yheight + 0.5*self.xwidth
		xs[mask2] = -self.xwidth /2; ys[mask2] = -ts[mask2] + 1.5*self.yheight + 1.0*self.xwidth
		ys[mask3] = -self.yheight/2; xs[mask3] =  ts[mask3] - 2.0*self.yheight - 1.5*self.xwidth
		if self.rcurv is not None:
			rs = (self.rcurv + xs) * np.sign(self.rcurv)
			angs = zs / np.abs(self.rcurv)
			xs = np.sign(self.rcurv) * rs * np.cos(angs) - self.rcurv
			zs = rs * np.sin(angs)
		return np.stack((xs,ys,zs), axis=1)
	def save(self, file):
		file.write(self.__class__.__name__+'\n')
		file.write("{}\n".format(self.dim))
		rcurv = self.rcurv
		if rcurv is None: rcurv = 0
		file.write("4 {} {} {} {}\n".format(self.xwidth, self.yheight, self.zmax, rcurv))

class Isotrop (Metric):
	def __init__(self):
		super().__init__(["dx","dy","dz"], ["[dir]","[dir]","[dir]"], "[dir]^3")
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
		super().__init__(["mu","phi"], ["[mu]","deg"], "[mu]deg")
	def transform(self, dirs):
		mus = dirs[:,2]
		phis = np.arctan2(dirs[:,1], dirs[:,0]) * 180/np.pi
		return np.stack((mus, phis), axis=1)
	def inverse_transform(self, tps):
		mus = tps[:,0]
		phis = tps[:,1]
		dxs = np.sqrt(1-mus**2) * np.cos(phis*np.pi/180)
		dys = np.sqrt(1-mus**2) * np.sin(phis*np.pi/180)
		dzs = mus
		return np.stack((dxs, dys, dzs), axis=1)