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
	def mean(self, parts, transform=True):
		if transform:
			vecs = self.transform(parts)
		else:
			vecs = parts
		return np.mean(parts, axis=0)
	def std(self, parts, transform=True):
		if transform:
			vecs = self.transform(parts)
		else:
			vecs = parts
		return np.std(parts, axis=0)

class SepVarMetric (Metric):
	def __init__(self, metric_E, metric_pos, metric_dir):
		super().__init__(metric_E.dim+metric_pos.dim+metric_dir.dim)
		self.E = metric_E
		self.pos = metric_pos
		self.dir = metric_dir
	def transform(self, parts):
		Es = self.E.transform(parts[:,0:0])
		poss = self.pos.transform(parts[:,1:3])
		dirs = self.dir.transform(parts[:,4:6])
		return np.hstack((Es, poss, dirs))
	def inverse_transform(self, vecs):
		Es = self.E.inverse_transform(vecs[:,0:self.E.dim-1])
		poss = self.pos.inverse_transform(vecs[:,self.E.dim:self.E.dim+self.pos.dim-1])
		dirs = self.dir.inverse_transform(vecs[:,self.E.dim+self.pos.dim:])
		return np.hstack((Es, poss, dirs))
	def jac(self, parts):
		jac_E = self.E.jac(parts[:,0:0])
		jac_pos = self.pos.jac(parts[:,1:3])
		jac_dir = self.dir.jac(parts[:,4:6])
		return jac_E*jac_pos*jac_dir

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

class Dir (Metric):
	def __init__(self):
		super().__init__(3)
	def mean(self, dirs, transform=False):
		if transform:
			vecs = self.transform(dirs)
		else:
			vecs = dirs
		mn = np.mean(dirs, axis=0)
		return mn / np.linalg.norm(mn)
	def std(self, dirs, transform=False):
		if transform:
			vecs = self.transform(dirs)
		else:
			vecs = dirs
		mn = self.mean(vecs, transform=False)
		return (np.sum((vecs-mn)**2)/(len(vecs)-1))**0.5