# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col


class T4Tally:
	varnames = ["x","y","z"]
	varmap = {name:idx for idx,name in enumerate(varnames)}
	units = ["cm","cm","cm"]

	def __init__(self, filename, tallyname, J=1):
		self.J = J
		try:
			file = open(filename, "r")
		except FileNotFoundError:
			print("Error: {} no encontrado".format(filename))
		# Buscar bloque de SCOREs
		for line in file:
			if "SCORE" in line:
				break
		else:
			print("Error: Archivo no tiene SCORE")
			file.close()
			return
		# Buscar tally deseado
		for line in file:
			if tallyname in line:
				break
			if "END_SCORE" in line:
				print("Error: No se encontro score {}".format(tallyname))
				file.close()
				return
		file.readline() # Response
		file.readline() # Estimator
		file.readline() # Energy grid
		# Leer grilla
		if not "EXTENDED_MESH" in file.readline():
			print("Error: La grilla debe ser EXTENDED_MESH")
			file.close()
			return
		file.readline() # WINDOW
		mins = np.double(file.readline().split())
		maxs = np.double(file.readline().split())
		Ns = list(map(int, file.readline().split()))
		grid1 = np.linspace(mins[0], maxs[0], Ns[0]+1)
		grid2 = np.linspace(mins[1], maxs[1], Ns[1]+1)
		grid3 = np.linspace(mins[2], maxs[2], Ns[2]+1)
		self.grids = [grid1, grid2, grid3]
		# Leer coordenadas
		if not "FRAME CARTESIAN" in file.readline():
			print("Error: Se debe tener FRAME CARTESIAN")
			file.close()
			return
		self.origin = np.double(file.readline().split())
		self.dx1 = np.double(file.readline().split())
		self.dx2 = np.double(file.readline().split())
		self.dx3 = np.double(file.readline().split())
		# Crear array 3D para almacenar tally
		I = []
		err = []
		# Buscar tally
		for line in file:
			if "SCORE NAME : "+tallyname in line:
				break
		else:
			print("Error: No se encontro tally {}".format(tallyname))
			file.close()
			return
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

	def plot(self, idx, cells=[0,0], **kwargs):
		if isinstance(idx, str):
			idx = self.varmap[idx]
		if not "xscale" in kwargs: kwargs["xscale"] = "linear"
		if not "yscale" in kwargs: kwargs["yscale"] = "log"
		slc = cells.copy()
		slc.insert(idx, slice(None))
		scores = self.J * self.I[tuple(slc)]
		errs = self.J * self.err[tuple(slc)]
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		#
		idxs = [0,1,2]
		idxs.remove(idx)
		lbl = str(self.grids[idxs[0]][cells[0]])+" <= "+self.varnames[idxs[0]]+" <= "+str(self.grids[idxs[0]][cells[0]+1])
		lbl += str(self.grids[idxs[1]][cells[1]])+" <= "+self.varnames[idxs[1]]+" <= "+str(self.grids[idxs[1]][cells[1]+1])
		plt.errorbar(grid, scores, errs, fmt='-s', label=lbl)
		plt.xscale(kwargs["xscale"])
		plt.yscale(kwargs["yscale"])
		plt.xlabel(r"${}\ [{}]$".format(self.varnames[idx], self.units[idx]))
		plt.ylabel("Tally")
		plt.grid()
		plt.legend()
		#
		return [plt.gcf(), [scores,errs]]

	def plot2D(self, idxs, cell=0, **kwargs):
		if isinstance(idxs[0], str):
			idxs = [self.varmap[idx] for idx in idxs]
		if not "scale" in kwargs: kwargs["scale"] = "log"
		idx = [0,1,2]
		idx.remove(idxs[0])
		idx.remove(idxs[1])
		idx = idx[0]
		slc = 2 * [slice(None)]
		slc.insert(idx, cell)
		scores = self.J * self.I[tuple(slc)]
		errs = self.J * self.err[tuple(slc)]
		if "fact" in kwargs:
			scores *= kwargs["fact"]
			errs *= kwargs["fact"]
		if idxs[1] > idxs[0]:
			scores = np.transpose(scores)
			errs = np.transpose(errs)
		#
		if kwargs["scale"] == "log": norm = col.LogNorm()
		else: norm = None
		extent = (self.grids[idxs[0]][0], self.grids[idxs[0]][-1], self.grids[idxs[1]][0], self.grids[idxs[1]][-1])
		plt.imshow(scores, extent=extent, cmap="jet", norm=norm, aspect='auto')
		plt.colorbar()
		title = "Tally"
		title += "\n"+str(self.grids[idx][cell])+" <= "+self.varnames[idx]+" <= "+str(self.grids[idx][cell+1])
		plt.title(title)
		plt.xlabel(r"${}\ [{}]$".format(self.varnames[idxs[0]], self.units[idxs[0]]))
		plt.ylabel(r"${}\ [{}]$".format(self.varnames[idxs[1]], self.units[idxs[1]]))
		#
		return [plt.gcf(), [scores,errs]]