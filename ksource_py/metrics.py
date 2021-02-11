# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class Metric:
	def __init__(self, dim):
		self.dim = dim
	def transform(self, parts):
		return parts
	def inverse_transform(self, vecs):
		return parts
	def jac(self, parts):
		return np.ones(len(parts))
	def mean(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		return self.inverse_transform(np.mean(vecs, axis=0))
	def std(self, parts=None, vecs=None):
		if vecs is None:
			vecs = self.transform(parts)
		return np.std(vecs, axis=0)

class SepVarMetric (Metric):
	def __init__(self, metric_E, metric_pos, metric_dir):
		super().__init__(metric_E.dim+metric_pos.dim+metric_dir.dim)
		self.E = metric_E
		self.pos = metric_pos
		self.dir = metric_dir
	def transform(self, parts):
		transf_Es = self.E.transform(parts[:,0:1])
		transf_poss = self.pos.transform(parts[:,1:4])
		transf_dirs = self.dir.transform(parts[:,4:7])
		return np.concatenate((transf_Es, transf_poss, transf_dirs), axis=1)
	def inverse_transform(self, vecs):
		Es = self.E.inverse_transform(vecs[:,0:self.E.dim])
		poss = self.pos.inverse_transform(vecs[:,self.E.dim:self.E.dim+self.pos.dim])
		dirs = self.dir.inverse_transform(vecs[:,self.E.dim+self.pos.dim:])
		return np.concatenate((Es, poss, dirs), axis=1)
	def jac(self, parts):
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

class Energy (Metric):
	def __init__(self):
		super().__init__(1)

class Lethargy (Metric):
	def __init__(self, E0):
		super().__init__(1)
		self.E0 = E0
	def transform(self, Es):
		return np.log(self.E0/Es)
	def inverse_transform(self, us):
		return self.E0 * np.exp(-us)
	def jac(self, Es):
		return 1/Es.reshape(-1)

class Vol (Metric):
	def __init__(self):
		super().__init__(3)

class SurfXY (Metric):
	def __init__(self, z):
		super().__init__(2)
		self.z = z
	def transform(self, poss):
		return poss[:,:2]
	def inverse_transform(self, poss):
		return np.hstack((poss, self.z))

class Isotrop (Metric):
	def __init__(self):
		super().__init__(3)
	def jac(self, parts):
		return np.ones(len(parts)) # ??? Calcular factor de cambio de unidades (volum a sup)
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
		super().__init__(2)
	def transform(self, dirs):
		thetas = np.arccos(dirs[:,2])
		phis = np.arctan2(dirs[:,1], dirs[:,0])
		return np.hstack((thetas, phis))
	def inverse_transform(self, tps):
		thetas = tps[:,0]
		phis = tps[:,1]
		xs = np.sin(thetas) * np.cos(phis)
		ys = np.sin(thetas) * np.sin(phis)
		zs = np.cos(thetas)
		return np.hstack((xs, ys, zs))
	def jac(self, dirs):
		return 1/np.sqrt(1-dirs[:,2]**2)