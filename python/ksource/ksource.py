# -*- coding: utf-8 -*-

from xml.etree import ElementTree as ET
from xml.dom import minidom
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from KDEpy import TreeKDE

from .plist import PList
from .geom import Geometry

np.set_printoptions(precision=3)

from .kde import optimize_bw,bw_silv

R_gaussian = 1 / (2*np.sqrt(np.pi)) # Roughness of gaussian kernel

STD_DEFAULT = 1


class KSource:
	def __init__(self, plist, geom, bw="silv", scaling=None, J=1):
		self.plist = plist
		self.geom = geom
		self.bw_method = None
		if isinstance(bw, str):
			self.bw_method = bw
			bw = 1.0
		elif isinstance(bw, np.ndarray):
			if np.ndim(bw) >= 2:
				raise ValueError("BW dimension must be < 2. Use scaling for anisotropic KDE.")
		self.kde = TreeKDE(bw=bw)
		if scaling is not None:
			scaling = np.array(scaling).ravel()
			if len(scaling)!=1 and len(scaling)!=geom.dim:
				raise ValueError("scaling len must be equal to geom dim.")
		else:
			scaling = np.ones(geom.dim)
		self.scaling = scaling
		self.J = J
		self.fitted = False

	def fit(self, N=-1, skip=0, **kwargs):
		parts,ws = self.plist.get(N, skip)
		if len(parts) == 0:
			raise Exception("No hay particulas para ajuste")
		N = len(parts)
		print("Usando {} particulas para ajuste".format(N))
		vecs = self.geom.transform(parts)
		self.N_eff = np.sum(ws)**2 / np.sum(ws**2)
		if self.bw_method is not None:
			print("Calculando bw ... ")
			scaling = self.geom.std(vecs=vecs)
			scaling[scaling == 0] = STD_DEFAULT
			self.scaling = scaling
			bw = optimize_bw(self.bw_method, vecs/self.scaling, ws, **kwargs)
			print("Hecho\nOptimal bw ({}) = {}".format(self.bw_method, np.reshape(bw, (-1,1)) * self.scaling))
			self.kde = TreeKDE(bw=bw)
		self.kde.fit(vecs/self.scaling, weights=ws)
		self.fitted = True

	def evaluate(self, parts):
		if self.fitted == False:
			raise Exception("Se debe ajustar (fit) antes de evaluar")
		vecs = self.geom.transform(parts)
		jacs = self.geom.jac(parts)
		evals = 1/np.prod(self.scaling) * self.kde.evaluate(vecs/self.scaling)
		errs = np.sqrt(evals * R_gaussian**self.geom.dim / (self.N_eff * np.mean(self.kde.bw) * np.prod(self.scaling)))
		evals *= self.J * jacs
		errs *= self.J * jacs
		return [evals, errs]

	def save(self, xmlfilename=None, bwfile=None, adjust_N=True):
		# Procesar argumentos
		if xmlfilename is None:
			xmlfilename = self.plist.filename.split('.')[0]+"_source.xml"
		if bwfile is None:
			bwfile = bwfilename = self.plist.filename.split('.')[0]+"_bws"
		elif isinstance(bwfile, str):
			bwfilename = bwfile
		else: # Asumo que es file object
			bwfilename = bwfile.name
		print("Archivo de definicion de fuente: {}".format(xmlfilename))
		bw = self.kde.bw
		if adjust_N: # Reajusto N_eff con factor de Silverman
			dim = self.geom.dim
			bw *= bw_silv(dim, self.plist.N) / bw_silv(dim, len(self.kde.data))
		# Construir arbol XML
		root = ET.Element("KSource")
		Jel = ET.SubElement(root, "J")
		Jel.set('units', '1/s')
		Jel.text = str(self.J)
		pltree = ET.SubElement(root, "PList")
		self.plist.save(pltree)
		gtree = ET.SubElement(root, "Geom")
		self.geom.save(gtree)
		ET.SubElement(root, "scaling").text = np.array_str(self.scaling)[1:-1]
		bwel = ET.SubElement(root, "BW")
		if np.isscalar(bw): # Ancho de banda constante
			bwel.set('variable', '0')
			bwel.text = str(bw)
		else: # Ancho de banda variable
			bwel.set('variable', '1')
			bw.astype("float32").tofile(bwfile, format="float32")
			print("Archivo de anchos de banda: {}".format(bwfilename))
			bwel.text = os.path.abspath(bwfilename)
		# Escribir archivo XML
		xmlstr = ET.tostring(root, encoding='utf8', method='xml')
		xmlstr = minidom.parseString(xmlstr).toprettyxml()
		with open(xmlfilename, "w") as file:
			file.write(xmlstr)
		return xmlfilename

	@staticmethod
	def load(xmlfilename):
		tree = ET.parse(xmlfilename)
		root = tree.getroot()
		J = np.double(root[0].text)
		plist = PList.load(root[1])
		geom = Geometry.load(root[2])
		scaling = np.array(root[3].text.split(), dtype="float64")
		bwel = root[4]
		if bool(int(bwel.attrib["variable"])):
			bw = np.fromfile(bwel.text, dtype="float32").astype("float64")
		else:
			bw = np.double(bwel.text)
		return KSource(plist, geom, bw, scaling, J=1)

	# Plot methods

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
			mask1 = np.logical_and.reduce(vec0/self.scaling <= self.kde.data, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.kde.data <= vec1/self.scaling, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		vecs = self.kde.data[:,idx][mask].reshape(-1,1)
		ws = self.kde.weights[mask]
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		scaling = self.scaling[idx]
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / bw_silv(self.geom.dim, self.N_eff)
		kde = TreeKDE(bw=bw)
		kde.fit(vecs, weights=ws)
		scores = 1/scaling * kde.evaluate(grid.reshape(-1,1)/scaling)
		errs = np.sqrt(scores * R_gaussian / (N_eff * np.mean(bw) * scaling))
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
			mask1 = np.logical_and.reduce(vec0/self.scaling <= self.kde.data, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.kde.data <= vec1/self.scaling, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		if sum(mask) == 0:
			raise Exception("No hay tracks en el rango pedido")
		vecs = self.kde.data[:,0][mask].reshape(-1,1)
		ws = self.kde.weights[mask]
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		scaling = self.scaling[0]
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / bw_silv(self.geom.dim, self.N_eff)
		kde = TreeKDE(bw=bw)
		kde.fit(vecs, weights=ws)
		grid = self.geom.ms[0].transform(grid_E)
		jacs = self.geom.ms[0].jac(grid_E)
		scores = 1/scaling * kde.evaluate(grid.reshape(-1,1)/scaling)
		errs = np.sqrt(scores * R_gaussian / (N_eff * np.mean(bw) * scaling))
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
		if isinstance(idxs[0], str): idxs[0] = self.geom.varmap[idxs[0]]
		if isinstance(idxs[1], str): idxs[1] = self.geom.varmap[idxs[1]]
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
		if isinstance(idxs[0], str): idxs[0] = self.geom.varmap[idxs[0]]
		if isinstance(idxs[1], str): idxs[1] = self.geom.varmap[idxs[1]]
		if not "scale" in kwargs: kwargs["scale"] = "linear"
		trues = np.array(len(self.kde.data)*[True])
		if vec0 is not None:
			mask1 = np.logical_and.reduce(vec0/self.scaling <= self.kde.data, axis=1)
		else:
			mask1 = trues
		if vec1 is not None:
			mask2 = np.logical_and.reduce(self.kde.data <= vec1/self.scaling, axis=1)
		else:
			mask2 = trues
		mask = np.logical_and(mask1, mask2)
		vecs = self.kde.data[:,idxs][mask]
		ws = self.kde.weights[mask]
		N_eff = np.sum(ws)**2 / np.sum(ws**2)
		scaling = self.scaling[idxs]
		bw = self.kde.bw * optimize_bw("silv", vecs, ws) / bw_silv(self.geom.dim, self.N_eff)
		kde = TreeKDE(bw=bw)
		kde.fit(vecs, weights=ws)
		grid = np.reshape(np.meshgrid(*grids),(2,-1)).T
		scores = 1/np.prod(scaling) * kde.evaluate(grid/scaling)
		errs = np.sqrt(scores * R_gaussian**2 / (N_eff * np.mean(bw) * np.prod(scaling)))
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