# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV, LeaveOneOut

class KSource:
	def __init__(self, metric, bw=1):
		self.metric = metric
		if isinstance(bw, str):
			self.bw = None
			self.bw_method = bw
		elif isinstance(bw, (list, tuple, np.ndarray)):
			self.bw = np.array(bw)
			self.bw_method = None
		elif np.isscalar(bw):
			self.bw = bw
			self.bw_method = None
		else:
			print("Error: Invalid bandwidth")
		self.kde = KernelDensity(bandwidth=1.0)
		self.std = None

	def fit(self, parts, weights, N_tot=None):
		self.vecs = self.metric.transform(parts)
		self.ws = weights
		if self.bw_method is not None:
			print("Calculating bw ...")
			self.bw = self.optimize_bw(parts, N_tot=N_tot, method=self.bw_method)
			print("Done. Optimal bw =", self.bw)
		self.kde.fit(self.vecs/self.bw, sample_weight=self.ws)

	def score(self, parts):
		vecs = self.metric.transform(parts)
		jacs = self.metric.jac(parts)
		scores = self.kde.score_samples(vecs/self.bw)
		return jacs * np.exp(scores)
		
	def optimize_bw(self, parts, N_tot=None, method='silv'):
		vecs = self.metric.transform(parts)
		#
		if method == 'silv': # Metodo de Silverman
			C_silv = 0.9397 # Cte de Silverman (revisar)
			std = self.metric.std(vecs=vecs)
			if N_tot is None:
				N_tot = len(parts)
			bw_silv = C_silv * std * N_tot**(-1/(4+self.metric.dim))
			return bw_silv
		#
		elif method == 'mlcv': # Metodo Maximum Likelihood Cross Validation
			C_silv = 0.9397 # Cte de Silverman (revisar)
			std = self.metric.std(vecs=vecs)
			N = len(parts)
			bw_silv = C_silv * std * N**(-1/(4+self.metric.dim)) # Regla del dedo de Silverman
			#
			nsteps = 20 # Cantidad de pasos para bw
			max_fact = 1.5 # Rango de bw entre bw_silv/max_fact y bw_silv*max_fact
			max_log = np.log10(max_fact)
			bw_grid = np.logspace(-max_log,max_log,nsteps) # Grilla de bandwidths
			#
			cv = 10
			grid = GridSearchCV(KernelDensity(kernel='gaussian'),
			                                  {'bandwidth': bw_grid},
			                                  cv=cv,
			                                  verbose=10,
			                                  n_jobs=8)
			grid.fit(vecs/bw_silv)
			bw_mlcv = bw_silv * grid.best_params_['bandwidth']
			#
			if N_tot is not None:
				bw_mlcv *= (N_tot/N)**(-1/(4+self.metric.dim))
			#
			return bw_mlcv
		else:
			print("Error: Invalid method")

	def plot_point(self, grid, idx, part0):
		parts = np.zeros((len(grid), 7))
		parts[:,idx] = grid
		part0[idx] = 0
		parts += part0
		scores = self.score(parts)
		#
		lbl = "part = "+str(part0)
		plt.plot(grid, scores, '-s', label=lbl)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel("??? [???]")
		plt.ylabel(r"$\Phi\ \left[ \frac{1}{???} \right]$")
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), scores]

	def plot_integr(self, grid, idx, vec0, vec1, jacs=None):
		mask1 = np.logical_and.reduce(vec0 < self.vecs, axis=1)
		mask2 = np.logical_and.reduce(self.vecs < vec1, axis=1)
		mask = np.logical_and(mask1, mask2)
		vecs = self.vecs[:,idx][mask].reshape(-1,1)
		ws = self.ws[mask]
		bw = self.bw[idx]
		kde = KernelDensity(bandwidth=1.0)
		kde.fit(vecs/bw, sample_weight=ws)
		scores = np.exp(kde.score_samples(grid.reshape(-1,1)/bw))
		if jacs is not None:
			scores *= jacs
		#
		lbl = str(vec0)+" < vec < "+str(vec1)
		plt.plot(grid, scores, '-s', label=lbl)
		plt.yscale('log')
		plt.xlabel("??? [???]")
		plt.ylabel(r"$\Phi\ \left[ \frac{1}{???} \right]$")
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), scores]

	def plot2D_point(self, grids, idxs, part0):
		parts = np.zeros((len(grid[0])*len(grid[1]), 7))
		parts[:,idxs] = np.meshgrid(grids).reshape(-1,2).T
		part0[idx] = 0
		parts += part0
		scores = self.score(parts)
		#
		plt.scatter(*parts[:,idxs].T, c=scores, edgecolors='k')
		title = r"$\Phi \left[ \frac{1}{???} \right]$"
		title += "\npart = "+str(part0)
		plt.title(title)
		plt.xlabel("??? [???]")
		plt.ylabel("??? [???]")
		#
		return [plt.gcf(), scores]

	def plot2D_integr(self, grids, idxs, vec0, vec1, jacs=None):
		mask1 = np.logical_and.reduce(vec0 < self.vecs, axis=1)
		mask1 = np.logical_and.reduce(self.vecs < vec1, axis=1)
		mask = np.logical_and(mask1, mask2)
		vecs = self.vecs[:,idxs][mask]
		ws = self.ws[masks]
		bw = self.bw[idx]
		kde = KernelDensity(bandwidth=1.0)
		kde.fit(vecs/bw, sample_weight=ws)
		grid = np.meshgrid(grids).reshape(-1,2).T
		scores = np.exp(kde.score_samples(grid/bw))
		if jacs is not None:
			scores *= jacs
		#
		plt.scatter(*grid.T, c=scores, edgecolors='k')
		title = r"$\Phi \left[ \frac{1}{???} \right]$"
		title += "\npart = "+str(part0)
		plt.title(title)
		plt.xlabel("??? [???]")
		plt.ylabel("??? [???]")
		#
		return [plt.gcf(), scores]