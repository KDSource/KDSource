# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import os
from PIL import Image as Im

from .plists import convert2mcpl,savessv,appendssv


def read_spectrum(spectrum=None):
	Es = []
	ws = []
	if spectrum is not None:
		file = open(spectrum, "r")
		for line in file:
		    try:
		        line = line.split(sep=',')
		        Es.append(np.double(line[0])/1000.)
		        ws.append(np.double(line[2]))
		    except:
		        pass
		file.close()
	else:
		Es.append(1)
		ws.append(1)
	Es = np.array(Es)
	ws = np.array(ws)
	if len(Es) != len(ws):
		raise Exception("Formato de espectro de decaimiento invalido")
	if len(Es) == 0:
		raise Exception("Espectro de decaimiento vacio")
	return [Es, ws]

class T4Tally:
	varnames = ["x","y","z"]
	varmap = {name:idx for idx,name in enumerate(varnames)}
	units = ["cm","cm","cm"]

	def __init__(self, outputfile, tallyname, spectrum=None, geomplot=None, J=1):
		self.J = J
		self.folder = os.path.dirname(outputfile)
		self.outputfile = outputfile
		self.tallyname = tallyname
		# Leer espectro
		self.Es,self.Ews = read_spectrum(spectrum)
		# Leer grafico de geometria
		if geomplot is not None:
			geomplot = np.array(Im.open(geomplot).convert('L').crop((25,26,514,514)))
		self.geomplot = geomplot
		# Leer tallies
		file = open(outputfile, "r")
		# Buscar bloque de SCOREs
		for line in file:
			if "SCORE" in line:
				break
		else:
			file.close()
			raise Exception("Archivo no tiene SCORE")
		# Buscar tally deseado
		for line in file:
			if tallyname in line:
				break
			if "END_SCORE" in line:
				file.close()
				raise Exception("No se encontro score {}".format(tallyname))
		# Buscar grilla
		for line in file:
			if "EXTENDED_MESH" in line:
				break
			if "NAME" in line or "END_SCORE" in line:
				file.close()
				raise Exception("No se encontro grilla EXTENDED_MESH")
		# Leer grillas
		buf = []
		idx = line.split().index("EXTENDED_MESH")
		buf.extend(line.split()[idx+1:]) # Acumular datos luego de EXTENDED_MESH
		for line in file:
			if "FRAME" in line:
				break
			buf.extend(line.split())
		idx = line.split().index("FRAME")
		buf.extend(line.split()[:idx]) # Acumular datos antes de FRAME
		if len(buf) != 10:
			file.close()
			raise Exception("No se pudo leer EXTENDED_MESH")
		mins = np.double(buf[1:4])
		maxs = np.double(buf[4:7])
		Ns = list(map(int, buf[7:10]))
		grid1 = np.linspace(mins[0], maxs[0], Ns[0]+1)
		grid2 = np.linspace(mins[1], maxs[1], Ns[1]+1)
		grid3 = np.linspace(mins[2], maxs[2], Ns[2]+1)
		self.grids = [grid1, grid2, grid3]
		# Leer coordenadas
		if not "FRAME CARTESIAN" in line:
			file.close()
			raise Exception("Se debe tener FRAME CARTESIAN")
		buf = []
		idx = line.split().index("CARTESIAN")
		buf.extend(line.split()[idx+1:]) # Acumular datos luego de CARTESIAN
		for line in file:
			found = False
			for search in ["NAME", "END_SCORE", "//", "/*"]:
				if search in line:
					if not found:
						idx = line.split().index(search)
						found = True
					else:
						idx2 = line.split().index(search)
						if idx2 < idx: idx = idx2
			if found: break
			buf.extend(line.split())
		buf.extend(line.split()[:idx]) # Acumular datos antes de NAME o END_SCORE
		if len(buf) != 12:
			file.close()
			raise Exception("No se pudo leer FRAME CARTESIAN")
		self.origin = np.double(buf[:3])
		self.dx1 = np.double(buf[3:6])
		self.dx2 = np.double(buf[6:9])
		self.dx3 = np.double(buf[9:12])
		# Buscar tally
		I = []
		err = []
		for line in file:
			if "SCORE NAME : "+tallyname in line:
				break
		else:
			file.close()
			raise Exception("No se encontro tally {}".format(tallyname))
		for line in file:
			if "Energy range" in line:
				break
		for line in file:
			line = line.split()
			if len(line) == 0:
				break
			I.append(np.double(line[1]))
			err.append(np.double(line[2]))
		if len(I) == np.prod(Ns):
			print("Tally {} leido exitosamente".format(tallyname))
		else:
			print("Lectura de tally {} incompleta".format(tallyname))
		self.I = np.reshape(I, Ns)
		self.err = np.reshape(err, Ns)

	def save_tracks(self, tracksfile=None):
		pt = "p"
		# Preparar lista de posiciones
		if tracksfile is None:
			tracksfile = self.folder + "/" + self.tallyname + ".ssv"
		grids = [(grid[:-1]+grid[1:])/2 for grid in self.grids]
		poss = np.reshape(np.meshgrid(*grids, indexing='ij'),(3,-1)).T
		ws = self.I.reshape(-1)
		poss = poss[ws>0]
		ws = ws[ws>0]
		ws /= ws.mean()
		# Preparar lista de particulas
		parts = np.zeros((len(poss), 7))
		parts[:,1:4] = poss
		parts[:,6] = 1
		np.random.shuffle(parts)
		# Generar lista de energias (ciclica)
		N = len(parts)
		nloops = int(np.ceil(N/len(self.Es)) + 1)
		Es = np.tile(self.Es,nloops)
		Ews = np.tile(self.Ews,nloops)
		# Guardar en archivo
		savessv(pt, [], [], tracksfile)
		for i in range(len(self.Es)):
			parts[:,0] = Es[i:i+N] # Energias
			mus = -1 + 2*np.random.rand(N)
			dxys = np.sqrt(1-mus**2)
			phis = -np.pi + 2*np.pi*np.random.rand(N)
			dirs = np.array([dxys*np.cos(phis), dxys*np.sin(phis), mus]).T
			parts[:,4:7] = dirs # Direcciones
			appendssv(pt, parts, ws*Ews[i:i+N], tracksfile)
		tracksfile = convert2mcpl(tracksfile, "ssv")
		print("Lista de tracks guardada exitosamente en {}".format(tracksfile))
		return tracksfile

	def plot(self, idx, cells=None, **kwargs):
		if isinstance(idx, str):
			idx = self.varmap[idx]
		if not "xscale" in kwargs: kwargs["xscale"] = "linear"
		if not "yscale" in kwargs: kwargs["yscale"] = "log"
		idxs = [0,1,2]
		idxs.remove(idx)
		if cells is None: # Promediar en las otras variables
			scores = np.mean(self.I, axis=idxs)
			errs = np.sqrt(np.sum(self.err**2, axis=idxs)) / (self.err.shape[idxs[0]]*self.err.shape[idxs[1]])
		else: # Graficar para las celdas indicada
			slc = cells.copy()
			slc.insert(idx, slice(None))
			scores = self.J * self.I[tuple(slc)]
			errs = self.J * self.err[tuple(slc)]
		if np.sum(scores) == 0:
			print("Tally nulo en la region a graficar")
			return [None, [scores,errs]]
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		if cells is None:
			lbl = str(self.grids[idxs[0]][0])+" <= "+self.varnames[idxs[0]]+" <= "+str(self.grids[idxs[0]][-1])
			lbl += str(self.grids[idxs[1]][0])+" <= "+self.varnames[idxs[1]]+" <= "+str(self.grids[idxs[1]][-1])
		else:
			lbl = str(self.grids[idxs[0]][cells[0]])+" <= "+self.varnames[idxs[0]]+" <= "+str(self.grids[idxs[0]][cells[0]+1])
			lbl += str(self.grids[idxs[1]][cells[1]])+" <= "+self.varnames[idxs[1]]+" <= "+str(self.grids[idxs[1]][cells[1]+1])
		plt.errorbar(grid, scores, errs, fmt='-s', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(self.varnames[idx], self.units[idx]))
		plt.ylabel("Tally")
		plt.grid()
		plt.legend()
		return [plt.gcf(), [scores,errs]]

	def plot2D(self, idxs, cell=None, geomplot=False, levelcurves=None, **kwargs):
		if isinstance(idxs[0], str):
			idxs = [self.varmap[idx] for idx in idxs]
		if not "scale" in kwargs: kwargs["scale"] = "log"
		idx = [0,1,2]
		idx.remove(idxs[0])
		idx.remove(idxs[1])
		idx = idx[0]
		if cell is None: # Promediar en la variable idx
			scores = np.mean(self.I, axis=idx)
			errs = np.sqrt(np.sum(self.err**2, axis=idx)) / self.err.shape[idx]
		else: # Graficar para la celda indicada
			slc = 2 * [slice(None)]
			slc.insert(idx, cell)
			scores = self.J * self.I[tuple(slc)]
			errs = self.J * self.err[tuple(slc)]
		if np.sum(scores) == 0:
			print("Tally nulo en la region a graficar")
			return [None, [scores,errs]]
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		if idxs[0] > idxs[1]:
			scores = np.transpose(scores)
			errs = np.transpose(errs)
		scores = np.rot90(scores)
		errs = np.rot90(errs)
		#
		if kwargs["scale"] == "log": norm = col.LogNorm()
		else: norm = None
		extent = (self.grids[idxs[0]][0], self.grids[idxs[0]][-1], self.grids[idxs[1]][0], self.grids[idxs[1]][-1])
		plt.imshow(scores, extent=extent, cmap="jet", norm=norm, aspect='auto')
		plt.colorbar()
		title = "Tally"
		if cell is None:
			title += "\n"+str(self.grids[idx][0])+" <= "+self.varnames[idx]+" <= "+str(self.grids[idx][-1])
		else:
			title += "\n"+str(self.grids[idx][cell])+" <= "+self.varnames[idx]+" <= "+str(self.grids[idx][cell+1])
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(self.varnames[idxs[0]], self.units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(self.varnames[idxs[1]], self.units[idxs[1]]))
		plt.tight_layout()
		#
		if levelcurves is not None:
			plt.contour(scores, levelcurves, extent=extent, linewidths=0.5)
		#
		if self.geomplot is not None and geomplot:
			ext = (extent[0], extent[1], extent[3], extent[2])
			for val in np.unique(self.geomplot):
				plt.contour(self.geomplot==val, [0.5], colors='black', extent=ext, linewidths=0.25)
		#
		return [plt.gcf(), [scores,errs]]