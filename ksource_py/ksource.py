# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import os
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV, LeaveOneOut

from .aux import R_gaussian,C_gaussian

np.set_printoptions(precision=3)

BW_DEFAULT = 1
STD_DEFAULT = 1


varnames = ["E", "x","y","z", "dx","dy","dz"]
varmap = {name:idx for idx,name in enumerate(varnames)}
units = ["MeV", "cm","cm","cm", "","",""]

class KSource:
	def __init__(self, plist, geom, bw="silv", J=1):
		self.plist = plist
		self.geom = geom
		if isinstance(bw, str):
			self.bw = None
			self.bw_method = bw
		elif isinstance(bw, (list, tuple, np.ndarray)):
			self.bw = np.array(bw)
			self.bw_method = None
		elif np.isscalar(bw):
			self.bw = bw * np.ones((geom.dim))
			self.bw_method = None
		else:
			raise Exceprion("Invalid bandwidth")
		self.kde = KernelDensity(bandwidth=1.0)
		self.J = J
		self.fitted = False

	def fit(self, N=-1, skip=0, **kwargs):
		parts,ws = self.plist.get(N, skip)
		if len(parts) == 0:
			raise Exceprion("No hay particulas para entrenamiento")
		print("Usando {} particulas para entrenamiento".format(len(parts)))
		self.vecs = self.geom.transform(parts)
		self.ws = ws
		if self.bw_method is not None:
			print("Calculando bw ... ")
			self.optimize_bw(**kwargs)
			print("Hecho\nOptimal bw ({}) = {}".format(self.bw_method, self.bw))
		self.bw[self.bw == 0] = BW_DEFAULT
		self.kde.fit(self.vecs/self.bw, sample_weight=self.ws)
		self.fitted = True

	def score(self, parts):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		vecs = self.geom.transform(parts)
		jacs = self.geom.jac(parts, bw=self.bw)
		scores = 1/np.prod(self.bw) * np.exp(self.kde.score_samples(vecs/self.bw))
		N_eff = np.sum(self.ws)**2 / np.sum(self.ws**2)
		errs = np.sqrt(scores * R_gaussian**self.geom.dim / (N_eff * np.prod(self.bw)))
		scores *= self.J * jacs
		errs *= self.J * jacs
		return np.array([scores, errs])

	def save(self, sourcefilename=None, bwfilename=None):
		if sourcefilename is None:
			sourcefilename = self.plist.filename.split('.')[0]+"_source.txt"
		if bwfilename is None:
			bwfilename = self.plist.filename.split('.')[0]+"_bws"
		print("Archivo de definicion de fuente: {}".format(sourcefilename))
		file = open(sourcefilename, "w")
		file.write("# J [1/s]:\n")
		file.write("%le\n" % (self.J))
		file.write("# PList:\n")
		self.plist.save(file)
		file.write("# Metric:\n")
		self.geom.save(file)
		if self.bw.ndim == 2: # Ancho de banda variable
			file.write("1\n")
			self.bw.astype("float64").tofile(bwfilename, format="float64")
			file.write(os.path.abspath(bwfilename)+"\n")
			print("Archivo de anchos de banda: {}".format(bwfilename))
		else:
			file.write("0\n")
			np.savetxt(file, self.bw[np.newaxis,:])
		file.close()
		return sourcefilename

	def optimize_bw(self, weightfun=None, maskfun=None, **kwargs):
		vecs = self.vecs.copy()
		ws = self.ws.copy()
		if weightfun is not None: ws *= weightfun(vecs)
		mask = (ws > 0)
		if maskfun is not None: mask = np.logical_and(mask, maskfun(vecs))
		vecs = vecs[mask]
		ws = ws[mask]
		N = len(vecs)
		if "N_tot" in kwargs:
			N_tot = kwargs["N_tot"]
		else:
			N_tot = N
		std = self.geom.std(vecs=vecs)
		#
		if self.bw_method == 'silv': # Metodo de Silverman
			bw_silv = C_gaussian(self.geom.dim) * std * N_tot**(-1/(4+self.geom.dim))
			self.bw = bw_silv
		#
		elif self.bw_method == 'mlcv': # Metodo Maximum Likelihood Cross Validation
			bw_silv = C_gaussian(self.geom.dim) * std * N**(-1/(4+self.geom.dim)) # Regla del dedo de Silverman
			if "shift_fact" in kwargs: bw_silv *= kwargs["shift_fact"]
			#
			if "nsteps" in kwargs: nsteps = kwargs["nsteps"] # Cantidad de pasos para bw
			else:  nsteps = 20
			if "max_fact" in kwargs: max_fact = kwargs["max_fact"] # Rango de bw entre bw_silv/max_fact y bw_silv*max_fact
			else: max_fact = 1.5
			max_log = np.log10(max_fact)
			bw_grid = np.logspace(-max_log,max_log,nsteps) # Grilla de bandwidths
			#
			if "cv" in kwargs: cv = kwargs["cv"] # CV folds
			else:  cv = 10
			grid = GridSearchCV(KernelDensity(kernel='gaussian'),
											  {'bandwidth': bw_grid},
											  cv=cv,
											  verbose=10,
											  n_jobs=-1)
			bw_silv[bw_silv == 0] = BW_DEFAULT
			grid.fit(vecs/bw_silv)
			plt.plot(bw_grid, np.exp(grid.cv_results_['mean_test_score']*cv/N))
			plt.xlabel("ancho de banda normalizado")
			plt.ylabel("mean CV score")
			plt.show()
			print("Verificar que el grafico presenta un maximo")
			bw_mlcv = bw_silv * grid.best_params_['bandwidth']
			bw_mlcv *= (N_tot/N)**(-1/(4+self.geom.dim))
			#
			self.bw = bw_mlcv
		#
		elif self.bw_method == 'knn': # Metodo K Nearest Neighbours
			if "batch_size" in kwargs: batch_size = kwargs["batch_size"] # Tamaño de batch
			else: batch_size = 10000
			if "K" in kwargs: K = kwargs["K"] # Cantidad de vecinos
			elif "seed" in kwargs:
				vs = vecs[:batch_size] / kwargs["seed"]
				ks = []
				for v in vs:
					dists2 = np.sum((vs - v)**2, axis=1)
					ks.append(np.sum(dists2 < 1))
				K = round((np.mean(ks)-1)*N_tot/batch_size)
			else: K = 10
			print("Usando K = %d"%K)
			k = round(K * batch_size/N_tot)
			if k == 0:
				print("Warning: k = K*batch_size/N_tot = 0. Se cambiara a k=1 <=> K={}".format(int(N_tot/batch_size)))	
				k = 1
			batches = np.ceil(N / batch_size).astype(int)
			batch_size = np.round(N / batches).astype(int)
			std[std == 0] = STD_DEFAULT
			vecs /= std
			bw_knn = np.zeros((0,self.geom.dim))
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
			raise Exceprion("Invalid method")

	def plot_point(self, grid, idx, part0, **kwargs):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
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
		lbl = "part = "+np.array_str(np.array(part0), precision=2)
		plt.errorbar(grid, scores, errs, fmt='-s', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(varnames[idx], units[idx]))
		plt.ylabel(r"$\Phi\ \left[ \frac{{{}}}{{{} s}} \right]$".format(self.plist.pt, self.geom.volunits))
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), [scores,errs]]

	def plot_integr(self, grid, idx, vec0=None, vec1=None, **kwargs):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		if isinstance(idx, str):
			idx = self.geom.varmap[idx]
		if not "xscale" in kwargs: kwargs["xscale"] = "linear"
		if not "yscale" in kwargs: kwargs["yscale"] = "log"
		trues = np.ones(len(self.vecs), dtype=bool)
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0 <= self.vecs, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.vecs <= vec1, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		vecs = self.vecs[:,idx][mask].reshape(-1,1)
		ws = self.ws[mask]
		bw = self.bw[idx]
		kde = KernelDensity(bandwidth=1.0)
		kde.fit(vecs/bw, sample_weight=ws)
		scores = 1/bw * np.exp(kde.score_samples(grid.reshape(-1,1)/bw))
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		errs = np.sqrt(scores * R_gaussian / (N_eff * bw))
		scores *= self.J
		errs *= self.J
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = np.array_str(np.array(vec0), precision=2)+" <= vec <= "+np.array_str(np.array(vec1), precision=2)
		plt.errorbar(grid, scores, errs, fmt='-s', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(self.geom.varnames[idx], self.geom.units[idx]))
		plt.ylabel(r"$\Phi\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt, self.geom.units[idx]))
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), [scores,errs]]

	def plot_E(self, grid_E, vec0=None, vec1=None, **kwargs):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		trues = np.ones(len(self.vecs), dtype=bool)
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0 <= self.vecs, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.vecs <= vec1, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		if sum(mask) == 0:
			raise Exception("No hay tracks en el rango pedido")
		vecs = self.vecs[:,0][mask].reshape(-1,1)
		ws = self.ws[mask]
		bw = self.bw[0]
		kde = KernelDensity(bandwidth=1.0)
		kde.fit(vecs/bw, sample_weight=ws)
		grid = self.geom.ms[0].transform(grid_E)
		jacs = self.geom.ms[0].jac(grid_E)
		scores = 1/bw * np.exp(kde.score_samples(grid.reshape(-1,1)/bw))
		errs = np.sqrt(scores * R_gaussian / (len(vecs) * bw))
		scores *= self.J * jacs
		errs *= self.J * jacs
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = np.array_str(np.array(vec0), precision=2)+" <= vec <= "+np.array_str(np.array(vec1), precision=2)
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
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		if isinstance(idxs[0], str):
			idxs = [varmap[idx] for idx in idxs]
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
		title = r"$\Phi\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt,self.geom.volunits)
		title += "\npart = "+np.array_str(np.array(part0), precision=2)
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(varnames[idxs[0]], units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(varnames[idxs[1]], units[idxs[1]]))
		#
		return [plt.gcf(), [scores,errs]]

	def plot2D_integr(self, grids, idxs, vec0=None, vec1=None, **kwargs):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		if isinstance(idxs[0], str):
			idxs = [self.geom.varmap[idx] for idx in idxs]
		if not "scale" in kwargs: kwargs["scale"] = "linear"
		trues = np.array(len(self.vecs)*[True])
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0 <= self.vecs, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.vecs <= vec1, axis=1)
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
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		errs = np.sqrt(scores * R_gaussian**2 / (N_eff * bw[0]*bw[1]))
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
		if self.geom.units[idxs[0]] == self.geom.units[idxs[1]]:
			units = self.geom.units[idxs[0]]+"^2"
		else:
			units = self.geom.units[idxs[0]] + self.geom.units[idxs[1]]
		title = r"$\Phi\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt,units)
		title += "\n"+np.array_str(np.array(vec0), precision=2)+" <= vec <= "+np.array_str(np.array(vec1), precision=2)
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(self.geom.varnames[idxs[0]], self.geom.units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(self.geom.varnames[idxs[1]], self.geom.units[idxs[1]]))
		#
		return [plt.gcf(), [scores,errs]]