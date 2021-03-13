# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class PList:
	def __init__(self, readformat, filename, pt='n', trasl=None, rot=None, switch_x2z=False, set_params=True):
		self.readformat = readformat
		exec("self.readfun = "+readformat+"_read")
		self.filename = filename
		self.pt = pt
		if trasl is not None:
			trasl = np.array(trasl)
		if rot is not None:
			rot = st.Rotation.from_rotvec(rot)
		self.trasl = trasl
		self.rot = rot
		self.x2z = switch_x2z
		self.params_set = False
		if set_params:
			self.set_params()
		else:
			self.I = self.p2 = 1.0
			self.N = 1

	def set_params(self):
		file = open(self.filename, "r")
		I = p2 = N = 0
		print("Reading file", self.filename, "...")
		for line in file:
			part_w = self.readfun(line)
			if part_w is not None: # Linea de texto es particula
				w = part_w[1]
				I += w
				p2 += w*w
				N += 1
		file.close()
		print("Done")
		print("I = {}\np2 = {}\nN = {}".format(I, p2, N))
		self.I = I
		self.p2 = p2
		self.N = N
		self.params_set = True

	def get(self, N=-1, skip=0):
		file = open(self.filename, "r")
		parts = []
		ws = []
		cont = 0
		for line in file:
			part_w = self.readfun(line)
			if part_w is not None: # Linea de texto es particula
				cont += 1
				if cont == skip: break
		cont = 0
		for line in file:
			part_w = self.readfun(line)
			if part_w is not None: # Linea de texto es particula
				part,w = part_w
				parts.append(part)
				ws.append(w)
				cont += 1
				if cont==N: break
		file.close()
		parts = np.array(parts)
		ws = np.array(ws)
		if self.trasl is not None: # Aplico traslacion
			if parts.shape[1] == 7:
				parts[:,1:4] += self.trasl
			if parts.shape[1] == 3:
				parts += self.trasl
		if self.rot is not None: # Aplico rotacion
			if parts.shape[1] == 7:
				parts[:,1:4] = self.rot.apply(parts[:,1:4]) # Posicion
				parts[:,4:7] = self.rot.apply(parts[:,4:7]) # Direccion
			if parts.shape[1] == 3:
				parts = self.rot.apply(parts)
		if self.x2z: # Aplico permutacion (x,y,z) -> (y,z,x)
			if parts.shape[1] == 7:
				E,x,y,z,dx,dy,dz = parts.T
				parts = np.stack((E,y,z,x,dy,dz,dx), axis=1)
			if parts.shape[1] == 3:
				x,y,z = parts.T
				parts = np.stack((y,z,x), axis=1)
		parts = parts[ws>0]
		ws = ws[ws>0]
		return [parts, ws]

	def save(self, file):
		file.write(self.plist.pt+'\n')
		file.write(self.plist.readformat+'\n')
		file.write(self.plist.filename+'\n')
		if self.trasl is not None: np.savetxt(file, self.trasl)
		else: file.write('\n')
		if self.rot is not None: np.savetxt(file, self.rot)
		else: file.write('\n')
		file.write("%d" % (self.switch_x2z))


def PTRAC_read(line):
	line = line.split()
	if len(line) == 9: # Linea de texto es particula
		[x,y,z,dx,dy,dz,E,w,t] = np.double(line)
		part = [E,x,y,z,dx,dy,dz]
		return [part,w]
	return None

def T4stock_read(line):
	line = line.split()
	if line[0] == "NEUTRON" or line[0] == "PHOTON": # Linea de texto es particula
		[E,x,y,z,dx,dy,dz,w] = np.double(line[1:])
		part = [E,x,y,z,dx,dy,dz]
		return [part,w]
	return None

def SSV_read(line):
	line = line.split()
	if len(line) == 8: # Linea de texto es particula
		[E,x,y,z,dx,dy,dz,w] = np.double(line)
		part = [E,x,y,z,dx,dy,dz]
		return [part,w]
	return None

def Decay_read(line):
	line = line.split()
	E = np.double(line[0])
	w = np.double(line[2])
	E = np.array([E])
	return [E,w]

def SSVtally_read(line):
	line = line.split()
	if len(line) == 4: # Linea de texto es particula
		[x,y,z,w] = np.double(line)
		pos = [x,y,z]
		return [pos,w]
	return None