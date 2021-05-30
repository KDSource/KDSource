# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import os

np.set_printoptions(precision=3)

from sklearn.model_selection import GridSearchCV
from sklearn.base import BaseEstimator
from sklearn.neighbors import NearestNeighbors
from KDEpy import TreeKDE, FFTKDE

# Clase necesaria para compatibilidad entre KDEpy y Scikit-Learn
class CompatKDE (TreeKDE, BaseEstimator):
    def __init__(self, kernel='gaussian', bw=1, norm=2.0):
        TreeKDE.__init__(self, kernel, bw, norm)
    def score(self, X, y=None):
        return np.sum(np.log(self.evaluate(X)))


# Nombres y unidades de las variables de particula
varnames = ["E", "x","y","z", "dx","dy","dz"]
varmap = {name:idx for idx,name in enumerate(varnames)}
units = ["MeV", "cm","cm","cm", "","",""]

STD_DEFAULT = 1

R_gaussian = 1 / (2*np.sqrt(np.pi)) # Roughness of gaussian kernel

def optimize_bw(bw_method, vecs=None, ws=None, **kwargs):
	if bw_method != "silv":
		vecs = np.array(vecs).copy()
		ws = np.array(ws).copy()
		if "weightfun" in kwargs: ws *= kwargs["weightfun"](vecs) # Aplico weightfun
		mask = (ws > 0)
		if "maskfun" in kwargs: mask = np.logical_and(mask, kwargs["maskfun"](vecs)) # Aplico maskfun
		vecs = vecs[mask]
		ws = ws[mask]
		N = len(vecs)
	if ws is not None: N_eff_data = np.sum(ws)**2 / np.sum(ws**2)
	if 'N_eff' in kwargs: N_eff = kwargs['N_eff']
	else: N_eff = N_eff_data
	if 'dim' in kwargs: dim = kwargs['dim']
	else: dim = vecs.shape[1]
	#
	if bw_method == 'silv': # Metodo de Silverman
		bw_silv = (4/(2+dim))**(1/(4+dim)) * N_eff**(-1/(4+dim))
		return bw_silv
	#
	elif bw_method == 'mlcv': # Metodo Maximum Likelihood Cross Validation
		if "bw_grid" in kwargs: bw_grid = kwargs["bw_grid"]
		else: # Crear grilla log-equiespaciada entre bw_seed/max_fact y bw_seed*max_fact
			if "seed" in kwargs: bw_seed = kwargs["seed"]
			else: bw_seed = optimize_bw("silv", N_eff=N_eff_data, dim=dim)
			if "nsteps" in kwargs: nsteps = kwargs["nsteps"] # Cantidad de pasos para bw
			else:  nsteps = 20
			if "max_fact" in kwargs: max_fact = kwargs["max_fact"]
			else: max_fact = 2
			max_log = np.log10(max_fact)
			bw_grid = np.logspace(-max_log,max_log,nsteps).reshape(-1,*np.ones(np.ndim(bw_seed),dtype=int)) * bw_seed # Grilla de bandwidths
		#
		if "cv" in kwargs: cv = kwargs["cv"] # CV folds
		else: cv = 10
		grid = GridSearchCV(CompatKDE(),
		                    {'bw': list(bw_grid)},
		                    cv=cv,
		                    verbose=10,
		                    n_jobs=-1)
		grid.fit(vecs, weights=ws)
		plt.plot(grid.cv_results_['mean_test_score'])
		plt.xlabel("grilla de anchos de banda")
		plt.ylabel("mean test score")
		plt.tight_layout()
		if "show" in kwargs: show = kwargs["show"]
		else: show = True
		if(show): plt.show()
		if grid.best_index_ in (0, len(grid.param_grid["bw"])-1):
			raise Exception("No se encotro maximo en el rango de bw seleccionado. Mueva la grilla e intente nuevamente.")
		bw_mlcv = grid.best_params_['bw']
		if N_eff_data != N_eff:
			bw_mlcv *= optimize_bw("silv", N_eff=N_eff, dim=dim) / optimize_bw("silv", N_eff=N_eff_data, dim=dim) # Reajusto N_eff con factor de Silverman
		return bw_mlcv
	#
	elif bw_method == 'knn': # Metodo K Nearest Neighbours
		if "batch_size" in kwargs: batch_size = kwargs["batch_size"] # Tamaño de batch
		else: batch_size = 10000
		batches = np.max(1, np.round(N / batch_size)).astype(int)
		if 'k' in kwargs: # Cantidad de vecinos en un batch
			k = kwargs['k']
			f_k = 1.0
			K = k * N_eff / batch_size
		else:
			if 'K' in kwargs: K = kwargs['K'] # Cantidad total de vecinos
			else: K = 100 # Valor default
			k_float = batch_size * K / N_eff
			k = int(round(k_float) if round(k_float)>0 else 2)
			f_k = k_float / k # Factor de correccion para k
		print("Usando: k = {} / batch_size = {}, f_k = {} ===> K_eff = {}".format(k, batch_size, f_k, K))
		bw_knn = np.array([])
		for batch,vs in enumerate(np.array_split(vecs, batches)):
			print("batch =", batch+1, "/", batches)
			knn = NearestNeighbors(n_neighbors=k, n_jobs=-1)
			knn.fit(vs)
			ds,idxs = knn.kneighbors(vs)
			bws = ds[:,-1] * (f_k)**(1/dim) # Selecciono k-esima columna y aplico factor de correccion
			bw_knn = np.concatenate((bw_knn, bws))
		return bw_knn
	#
	else:
		raise Exception("bw_method invalido. Validos: 'silv', 'mlcv', 'knn'")

class KSource:
	def __init__(self, plist, geom, bw="silv", J=1):
		self.plist = plist
		self.geom = geom
		self.bw_method = None
		if isinstance(bw, str):
			self.bw_method = bw
			bw = 1.0
		self.kde = TreeKDE(bw=bw)
		self.J = J
		self.fitted = False

	def fit(self, N=-1, skip=0, std=None, **kwargs):
		parts,ws = self.plist.get(N, skip)
		if len(parts) == 0:
			raise Exception("No hay particulas para ajuste")
		N = len(parts)
		print("Usando {} particulas para ajuste".format(N))
		vecs = self.geom.transform(parts)
		if std is None: std = self.geom.std(vecs=vecs)
		else: std = np.array(std)
		assert len(std) == self.geom.dim
		std[std == 0] = STD_DEFAULT
		self.std = std
		self.N_eff = np.sum(ws)**2 / np.sum(ws**2)
		if self.bw_method is not None:
			print("Calculando bw ... ")
			bw = optimize_bw(self.bw_method, vecs/self.std, ws, **kwargs)
			print("Hecho\nOptimal bw ({}) = {}".format(self.bw_method, np.reshape(bw, (-1,1)) * self.std))
			self.kde = TreeKDE(bw=bw)
		self.kde.fit(vecs/self.std, weights=ws)
		self.fitted = True

	def evaluate(self, parts):
		if self.fitted == False:
			raise Exception("Se debe ajustar (fit) antes de evaluar")
		vecs = self.geom.transform(parts)
		jacs = self.geom.jac(parts)
		evals = 1/np.prod(self.std) * self.kde.evaluate(vecs/self.std)
		errs = np.sqrt(evals * R_gaussian**self.geom.dim / (self.N_eff * np.mean(self.kde.bw) * np.prod(self.std)))
		evals *= self.J * jacs
		errs *= self.J * jacs
		return [evals, errs]

	def save(self, sourcefilename=None, bwfile=None, adjust_N=True):
		if sourcefilename is None:
			sourcefilename = self.plist.filename.split('.')[0]+"_source.txt"
		if bwfile is None:
			bwfile = bwfilename = self.plist.filename.split('.')[0]+"_bws"
		elif isinstance(bwfile, str):
			bwfilename = bwfile
		else: # Asumo que es file object
			bwfilename = bwfile.name
		print("Archivo de definicion de fuente: {}".format(sourcefilename))
		with open(sourcefilename, "w") as file:
			file.write("# J [1/s]:\n")
			file.write("%le\n" % (self.J))
			file.write("# PList:\n")
			self.plist.save(file)
			file.write("# Metric:\n")
			self.geom.save(file)
			bw = np.reshape(self.kde.bw, (-1,1)) * self.std
			if adjust_N:
				N_eff_tot = self.N_eff * self.plist.N / len(self.kde.data)
				dim = self.geom.dim
				bw *= optimize_bw("silv", N_eff=N_eff_tot, dim=dim) / optimize_bw("silv", N_eff=self.N_eff, dim=dim) # Reajusto N_eff con factor de Silverman
			if len(bw) == 1: # Ancho de banda constante
				file.write("0\n")
				np.savetxt(file, bw)
			else: # Ancho de banda variable
				file.write("1\n")
				bw.astype("float32").tofile(bwfile, format="float32")
				file.write(os.path.abspath(bwfilename)+"\n")
				print("Archivo de anchos de banda: {}".format(bwfilename))
		return sourcefilename

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
		scores,errs = self.evaluate(parts)
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = "part = "+np.array_str(np.array(part0), precision=2)
		plt.errorbar(grid, scores, errs, fmt='-', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(varnames[idx], units[idx]))
		plt.ylabel(r"$J\ \left[ \frac{{{}}}{{{} s}} \right]$".format(self.plist.pt, self.geom.volunits))
		plt.grid()
		plt.legend()
		plt.tight_layout()
		#
		return [plt.gcf(), [scores,errs]]

	def plot_integr(self, grid, idx, vec0=None, vec1=None, **kwargs):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		if isinstance(idx, str):
			idx = self.geom.varmap[idx]
		if not "xscale" in kwargs: kwargs["xscale"] = "linear"
		if not "yscale" in kwargs: kwargs["yscale"] = "log"
		trues = np.ones(len(self.kde.data), dtype=bool)
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0/self.std <= self.kde.data, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.kde.data <= vec1/self.std, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		vecs = self.kde.data[:,idx][mask].reshape(-1,1)
		ws = self.kde.weights[mask]
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		std = self.std[idx]
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / optimize_bw("silv", N_eff=self.N_eff, dim=self.geom.dim)
		kde = TreeKDE(bw=bw)
		kde.fit(vecs, weights=ws)
		scores = 1/std * kde.evaluate(grid.reshape(-1,1)/std)
		errs = np.sqrt(scores * R_gaussian / (N_eff * np.mean(bw) * std))
		scores *= self.J * np.sum(ws)/np.sum(self.kde.weights)
		errs *= self.J * np.sum(ws)/np.sum(self.kde.weights)
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = np.array_str(np.array(vec0), precision=2)+" <= vec <= "+np.array_str(np.array(vec1), precision=2)
		plt.errorbar(grid, scores, errs, fmt='-', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(self.geom.varnames[idx], self.geom.units[idx]))
		plt.ylabel(r"$J\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt, self.geom.units[idx]))
		plt.grid()
		plt.legend()
		plt.tight_layout()
		#
		return [plt.gcf(), [scores,errs]]

	def plot_E(self, grid_E, vec0=None, vec1=None, **kwargs):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		trues = np.ones(len(self.kde.data), dtype=bool)
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0/self.std <= self.kde.data, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.kde.data <= vec1/self.std, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		if sum(mask) == 0:
			raise Exception("No hay tracks en el rango pedido")
		vecs = self.kde.data[:,0][mask].reshape(-1,1)
		ws = self.kde.weights[mask]
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		std = self.std[0]
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / optimize_bw("silv", N_eff=self.N_eff, dim=self.geom.dim)
		kde = TreeKDE(bw=bw)
		kde.fit(vecs, weights=ws)
		grid = self.geom.ms[0].transform(grid_E)
		jacs = self.geom.ms[0].jac(grid_E)
		scores = 1/std * kde.evaluate(grid.reshape(-1,1)/std)
		errs = np.sqrt(scores * R_gaussian / (N_eff * np.mean(bw) * std))
		scores *= self.J * np.sum(ws)/np.sum(self.kde.weights) * jacs
		errs *= self.J * np.sum(ws)/np.sum(self.kde.weights) * jacs
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		lbl = np.array_str(np.array(vec0), precision=2)+" <= vec <= "+np.array_str(np.array(vec1), precision=2)
		plt.errorbar(grid_E, scores, errs, fmt='-', label=lbl)
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel(r"$E\ [MeV]$")
		plt.ylabel(r"$J\ \left[ \frac{{{}}}{{MeV\ s}} \right]$".format(self.plist.pt))
		plt.grid()
		plt.legend()
		plt.tight_layout()
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
		scores,errs = self.J * self.evaluate(parts)
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		xx = np.concatenate((grids[0][:1], (grids[0][1:]+grids[0][:-1])/2, grids[0][-1:]))
		yy = np.concatenate((grids[1][:1], (grids[1][1:]+grids[1][:-1])/2, grids[1][-1:]))
		if kwargs["scale"] == "log": norm = col.LogNorm()
		else: norm = None
		plt.pcolormesh(xx, yy, scores.reshape(len(grids[1]), len(grids[0])), cmap="jet", norm=norm)
		plt.colorbar()
		title = r"$J\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt,self.geom.volunits)
		title += "\npart = "+np.array_str(np.array(part0), precision=2)
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(varnames[idxs[0]], units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(varnames[idxs[1]], units[idxs[1]]))
		plt.tight_layout()
		#
		return [plt.gcf(), [scores,errs]]

	def plot2D_integr(self, grids, idxs, vec0=None, vec1=None, **kwargs):
		if self.fitted == False:
			raise Exception("Se debe fittear antes de evaluar")
		if isinstance(idxs[0], str):
			idxs = [self.geom.varmap[idx] for idx in idxs]
		if not "scale" in kwargs: kwargs["scale"] = "linear"
		trues = np.array(len(self.kde.data)*[True])
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0/self.std <= self.kde.data, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.kde.data <= vec1/self.std, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		vecs = self.kde.data[:,idxs][mask]
		ws = self.kde.weights[mask]
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		std = self.std[idxs]
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / optimize_bw("silv", N_eff=self.N_eff, dim=self.geom.dim)
		kde = TreeKDE(bw=bw)
		kde.fit(vecs, weights=ws)
		grid = np.reshape(np.meshgrid(*grids),(2,-1)).T
		scores = 1/np.prod(std) * kde.evaluate(grid/std)
		errs = np.sqrt(scores * R_gaussian**2 / (N_eff * np.mean(bw) * np.prod(std)))
		scores *= self.J * np.sum(ws)/np.sum(self.kde.weights)
		errs *= self.J * np.sum(ws)/np.sum(self.kde.weights)
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		xx = np.concatenate((grids[0][:1], (grids[0][1:]+grids[0][:-1])/2, grids[0][-1:]))
		yy = np.concatenate((grids[1][:1], (grids[1][1:]+grids[1][:-1])/2, grids[1][-1:]))
		if kwargs["scale"] == "log": norm = col.LogNorm()
		else: norm = None
		plt.pcolormesh(xx, yy, scores.reshape(len(grids[1]), len(grids[0])), cmap="jet", norm=norm)
		plt.colorbar()
		if self.geom.units[idxs[0]] == self.geom.units[idxs[1]]:
			units = self.geom.units[idxs[0]]+"^2"
		else:
			units = self.geom.units[idxs[0]] + self.geom.units[idxs[1]]
		title = r"$J\ \left[ \frac{{{}}}{{{}\ s}} \right]$".format(self.plist.pt,units)
		title += "\n"+np.array_str(np.array(vec0), precision=2)+" <= vec <= "+np.array_str(np.array(vec1), precision=2)
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(self.geom.varnames[idxs[0]], self.geom.units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(self.geom.varnames[idxs[1]], self.geom.units[idxs[1]]))
		plt.tight_layout()
		#
		return [plt.gcf(), [scores,errs]]