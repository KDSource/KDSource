# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV, LeaveOneOut

varnames = ["E", "x","y","z", "dx","dy","dz"]
varmap = {name:idx for idx,name in enumerate(varnames)}
units = ["MeV", "cm","cm","cm", "","",""]

R_gaussian = 1 / (2*np.pi) # Roughness of gaussian kernel

class KSource:
	def __init__(self, metric, bw="silv", J=1):
		self.metric = metric
		if isinstance(bw, str):
			self.bw = None
			self.bw_method = bw
		elif isinstance(bw, (list, tuple, np.ndarray)):
			self.bw = np.array(bw)
			self.bw_method = None
		elif np.isscalar(bw):
			self.bw = bw * np.ones((metric.dim))
			self.bw_method = None
		else:
			print("Error: Invalid bandwidth")
		self.kde = KernelDensity(bandwidth=1.0)
		self.J = J

	def fit(self, plist, N, skip=0, **kwargs):
		self.plist = plist
		parts,ws = plist.get(N=N, skip=skip)
		if len(parts) == 0:
			print("Error: No se pudieron obtener particulas para entrenamiento")
			return
		print("Usando {} particulas para entrenamiento".format(len(parts)))
		self.vecs = self.metric.transform(parts)
		self.ws = ws
		if self.bw_method is not None:
			print("Calculando bw ... ")
			self.optimize_bw(**kwargs)
			print("Hecho\nOptimal bw ({}) = {}".format(self.bw_method, self.bw))
		self.kde.fit(self.vecs/self.bw, sample_weight=self.ws)

	def score(self, parts):
		vecs = self.metric.transform(parts)
		jacs = self.metric.jac(parts, bw=self.bw)
		scores = 1/np.prod(self.bw) * np.exp(self.kde.score_samples(vecs/self.bw))
		errs = np.sqrt(scores * R_gaussian**self.metric.dim / (len(parts) * np.prod(self.bw)))
		scores *= self.J * jacs
		errs *= self.J * jacs
		return np.array([scores, errs])

	def save_bw(self, bwfilename, append=False):
		if append:
			bwfilename = open(bwfilename, "a")
		np.savetxt(bwfilename, self.bw.reshape(-1,self.metric.dim))
		
	def optimize_bw(self, **kwargs):
		vecs = self.vecs.copy()
		N = len(vecs)
		if "N_tot" in kwargs:
			N_tot = kwargs["N_tot"]
		else:
			N_tot = N
		std = self.metric.std(vecs=vecs)
		#
		if self.bw_method == 'silv': # Metodo de Silverman
			C_silv = 0.9397 # Cte de Silverman (revisar)
			bw_silv = C_silv * std * N_tot**(-1/(4+self.metric.dim))
			self.bw = bw_silv
		#
		elif self.bw_method == 'mlcv': # Metodo Maximum Likelihood Cross Validation
			C_silv = 0.9397 # Cte de Silverman (revisar)
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
			plt.plot(bw_grid, np.exp(grid.cv_results_['mean_test_score']*cv/N))
			plt.xlabel("ancho de banda normalizado")
			plt.ylabel("mean CV score")
			plt.show()
			bw_mlcv = bw_silv * grid.best_params_['bandwidth']
			bw_mlcv *= (N_tot/N)**(-1/(4+self.metric.dim))
			#
			self.bw = bw_mlcv
		#
		elif self.bw_method == 'knn': # Metodo K Nearest Neighbours
			K = 10
			batch_size = 10000
			k = round(K * batch_size/N_tot)
			if k == 0:
				print("Warning: k = K*batch_size/N_tot = 0. Se cambiara a k=1 <=> K={}".format(N_tot/batch_size))
				k = 1
			batches = int(N / batch_size)
			vecs /= std
			bw_knn = np.zeros((0,self.metric.dim))
			for batch in range(batches):
				print("batch =", batch+1, "/", batches)
				if batch < batches-1:
					vs = vecs[batch*batch_size:(batch+1)*batch_size]
				else:
					vs = vecs[batch*batch_size:]
				bws = []
				for v in vs:
					dists2 = np.sum((vs - v)**2, axis=1)
					bws.append(np.sqrt(np.partition(dists2, k)[k]))
				bws = std * np.array(bws)[:,np.newaxis] / (N_tot/len(vs))**(1/len(std))
				bw_knn = np.concatenate((bw_knn, bws), axis=0)
			self.bw = bw_knn
		#
		else:
			print("Error: Invalid method")

	def plot_point(self, grid, idx, part0, **kwargs):
		if isinstance(idx, str):
			idx = varmap[idx]
		if not "xscale" in kwargs: kwargs["xscale"] = "linear"
		if not "yscale" in kwargs: kwargs["yscale"] = "log"
		parts = np.zeros((len(grid), 7))
		parts[:,idx] = grid
		part0[idx] = 0
		parts += part0
		scores,errs = self.score(parts)
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = "part = "+str(part0)
		plt.errorbar(grid, scores, errs, fmt='-s', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(varnames[idx], units[idx]))
		plt.ylabel(r"$\Phi\ \left[ \frac{{{}}}{{{} s}} \right]$".format(self.plist.pt, self.metric.volunits))
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), [scores,errs]]

	def plot_integr(self, grid, idx, vec0=None, vec1=None, **kwargs):
		if isinstance(idx, str):
			idx = self.metric.varnames[idx]
		if not "xscale" in kwargs: kwargs["xscale"] = "linear"
		if not "yscale" in kwargs: kwargs["yscale"] = "log"
		trues = np.array(len(self.vecs)*[True])
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0 < self.vecs, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.vecs < vec1, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		vecs = self.vecs[:,idx][mask].reshape(-1,1)
		ws = self.ws[mask]
		bw = self.bw[idx]
		kde = KernelDensity(bandwidth=1.0)
		kde.fit(vecs/bw, sample_weight=ws)
		scores = 1/bw * np.exp(kde.score_samples(grid.reshape(-1,1)/bw))
		errs = np.sqrt(scores * R_gaussian / (len(vecs) * bw))
		scores *= self.J
		errs *= self.J
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = str(vec0)+" < vec < "+str(vec1)
		plt.errorbar(grid, scores, errs, fmt='-s', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(self.metric.varnames[idx], self.metric.units[idx]))
		plt.ylabel(r"$\Phi\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt, self.metric.units[idx]))
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), [scores,errs]]

	def plot_E(self, grid_E, vec0=None, vec1=None, **kwargs):
		trues = np.array(len(self.vecs)*[True])
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0 < self.vecs, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.vecs < vec1, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		if sum(mask) == 0:
			print("Error: No hay tracks en el rango pedido")
			return
		vecs = self.vecs[:,0][mask].reshape(-1,1)
		ws = self.ws[mask]
		bw = self.bw[0]
		kde = KernelDensity(bandwidth=1.0)
		kde.fit(vecs/bw, sample_weight=ws)
		grid = self.metric.E.transform(grid_E)
		jacs = self.metric.E.jac(grid_E)
		scores = 1/bw * np.exp(kde.score_samples(grid.reshape(-1,1)/bw))
		errs = np.sqrt(scores * R_gaussian / (len(vecs) * bw))
		scores *= self.J * jacs
		errs *= self.J * jacs
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = str(vec0)+" < vec < "+str(vec1)
		plt.errorbar(grid_E, scores, errs, fmt='-s', label=lbl)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel(r"$E\ [MeV]$")
		plt.ylabel(r"$\Phi\ \left[ \frac{{{}}}{{MeV\ s}} \right]$".format(self.plist.pt))
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), [scores,errs]]

	def plot2D_point(self, grids, idxs, part0, **kwargs):
		if isinstance(idxs[0], str):
			idxs = [varnames[idx] for idx in idxs]
		if not "scale" in kwargs: kwargs["scale"] = "linear"
		parts = np.zeros((len(grids[0])*len(grids[1]), 7))
		parts[:,idxs] = np.reshape(np.meshgrid(*grids), (2,-1)).T
		part0 = np.array(part0)
		part0[idxs] = 0
		parts += part0
		scores,errs = self.J * self.score(parts)
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		xx = np.concatenate((grids[0][:1], (grids[0][1:]+grids[0][:-1])/2, grids[0][-1:]))
		yy = np.concatenate((grids[1][:1], (grids[1][1:]+grids[1][:-1])/2, grids[1][-1:]))
		if kwargs["scale"] == "log": norm = col.LogNorm()
		else: norm = None
		plt.pcolormesh(xx, yy, scores.reshape(len(grids[0]), len(grids[1])), cmap="jet", norm=norm)
		plt.colorbar()
		title = r"$\Phi\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt,self.metric.volunits)
		title += "\npart = "+str(part0)
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(varnames[idxs[0]], units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(varnames[idxs[1]], units[idxs[1]]))
		#
		return [plt.gcf(), [scores,errs]]

	def plot2D_integr(self, grids, idxs, vec0=None, vec1=None, **kwargs):
		if isinstance(idxs[0], str):
			idxs = [self.metric.varnames[idx] for idx in idxs]
		if not "scale" in kwargs: kwargs["scale"] = "linear"
		trues = np.array(len(self.vecs)*[True])
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0 < self.vecs, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.vecs < vec1, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		vecs = self.vecs[:,idxs][mask]
		ws = self.ws[mask]
		bw = self.bw[idxs]
		kde = KernelDensity(bandwidth=1.0)
		kde.fit(vecs/bw, sample_weight=ws)
		grid = np.reshape(np.meshgrid(*grids),(2,-1)).T
		scores = 1/np.prod(bw) * np.exp(kde.score_samples(grid/bw))
		errs = np.sqrt(scores * R_gaussian**2 / (len(vecs) * bw[0]*bw[1]))
		scores *= self.J
		errs *= self.J
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		xx = np.concatenate((grids[0][:1], (grids[0][1:]+grids[0][:-1])/2, grids[0][-1:]))
		yy = np.concatenate((grids[1][:1], (grids[1][1:]+grids[1][:-1])/2, grids[1][-1:]))
		if kwargs["scale"] == "log": norm = col.LogNorm()
		else: norm = None
		plt.pcolormesh(xx, yy, scores.reshape(len(grids[0]), len(grids[1])), cmap="jet", norm=norm)
		plt.colorbar()
		if self.metric.units[idxs[0]] == self.metric.units[idxs[1]]:
			units = self.metric.units[idxs[0]]+"^2"
		else:
			units = self.metric.units[idxs[0]] + self.metric.units[idxs[1]]
		title = r"$\Phi\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt,units)
		title += "\n"+str(vec0)+" < vec < "+str(vec1)
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(self.metric.varnames[idxs[0]], self.metric.units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(self.metric.varnames[idxs[1]], self.metric.units[idxs[1]]))
		#
		return [plt.gcf(), [scores,errs]]