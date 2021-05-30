# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import os
from KDEpy import TreeKDE

np.set_printoptions(precision=3)

from .kde import optimize_bw,bw_silv

R_gaussian = 1 / (2*np.sqrt(np.pi)) # Roughness of gaussian kernel

STD_DEFAULT = 1


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
			if adjust_N: # Reajusto N_eff con factor de Silverman
				dim = self.geom.dim
				bw *= bw_silv(dim, self.plist.N) / bw_silv(dim, len(self.kde.data))
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
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / bw_silv(self.geom.dim, self.N_eff)
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
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / bw_silv(self.geom.dim, self.N_eff)
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
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / bw_silv(self.geom.dim, self.N_eff)
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