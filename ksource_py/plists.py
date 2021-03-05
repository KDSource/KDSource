# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class PList:
	def __init__(self, readformat, filename, pt='n', trasl=None, rot=None, switch_x2z=False, set_params=True, **kwargs):
		if isinstance(readformat, str):
			self.readformat = readformat
			exec("self.readfun = "+readformat+"_read")
		else:
			self.readfun = readformat
			self.readformat = readformat.__name__[:-5]
		self.filename = filename
		self.pt = pt
		if trasl is not None:
			trasl = np.array(trasl)
		if rot is not None:
			rot = st.Rotation.from_rotvec(rot)
		self.trasl = trasl
		self.rot = rot
		self.x2z = switch_x2z
		self.start = None
		if "start" in kwargs: self.start = kwargs["start"]
		self.end = None
		if "end" in kwargs: self.end = kwargs["end"]
		if set_params:
			self.set_params()
		else:
			self.I = self.p2 = 1.0
			self.N = 1

	def set_params(self):
		try:
			file = open(self.filename, "r")
		except FileNotFoundError:
			print("Error: {} no encontrado".format(self.filename))
		I = p2 = N = 0
		if self.start is not None: read = False
		else: read = True
		print("Reading file", self.filename, "...")
		for line in file:
			if not read:
				if self.start in line: read = True
			else:
				if self.end is not None and self.end in line: break
				part_w = self.readfun(line)
				if part_w is not None: # Linea de texto es particula
					w = part_w[1]
					I += w
					p2 += w*w
					N += 1
		file.close()
		print("Done")
		self.I = I
		self.p2 = p2
		self.N = N

	def get(self, N=-1, skip=0):
		try:
			file = open(self.filename, "r")
		except FileNotFoundError:
			print("Error: {} no encontrado".format(self.filename))
		parts = []
		ws = []
		cont = 0
		if self.start is not None: read = False
		else: read = True
		for line in file:
			if not read:
				if self.start in line: read = True
			else:
				if self.end is not None and self.end in line: break
				part_w = self.readfun(line)
				if part_w is not None: # Linea de texto es particula
					part,w = part_w
					if cont >= skip: 
						parts.append(part)
						ws.append(w)
					cont += 1
					if N>0 and cont==skip+N:
						break
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

def T4tally_read(line):
	line = line.split()
	if len(line) and len(line[0]) and (line[0][0]=="(" and line[0][-1]==")"): # Linea de texto es particula
		[x,y,z] = np.double(line[0][1:-1].split(sep=","))
		part = [x,y,z]
		w = np.double(line[1])
		return [part,w]
	return None

def Decay_read(line):
	line = line.split()
	E = np.double(line[0])
	w = np.double(line[2])
	E = np.array([E])
	return [E,w]