# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class PList:
	def __init__(self, readfun, filename, pt="n", trasl=None, rot=None, switch_x2z=False, set_params=True):
		self.filename = filename
		self.pt = pt
		self.read = readfun
		if trasl is not None:
			trasl = np.array(trasl)
		if rot is not None:
			rot = st.Rotation.from_rotvec(rot)
		self.trasl = trasl
		self.rot = rot
		self.x2z = switch_x2z
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
			part_w = self.read(line)
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
		file = open(self.filename, "r")
		parts = []
		ws = []
		cont = 0
		for line in file:
			part_w = self.read(line)
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
			parts[:,1:4] += self.trasl
		if self.rot is not None: # Aplico rotacion
			parts[:,1:4] = self.rot.apply(parts[:,1:4]) # Posicion
			parts[:,4:7] = self.rot.apply(parts[:,4:7]) # Direccion
		if self.x2z: # Aplico permutacion (x,y,z) -> (y,z,x)
			E,x,y,z,dx,dy,dz = parts.T
			parts = np.stack((E,y,z,x,dy,dz,dx), axis=1)
		return [parts, ws]

def SSV_read(line):
	line = line.split()
	if len(line) == 8: # Linea de texto es particula
		[E,x,y,z,dx,dy,dz,w] = np.double(line)
		part = [E,x,y,z,dx,dy,dz]
		return [part,w]
	return None

def PTRAC_read(line):
	line = line.split()
	if len(line) == 9: # Linea de texto es particula
		[x,y,z,dx,dy,dz,E,w,t] = np.double(line)
		part = [E,x,y,z,dx,dy,dz]
		return [part,w]
	return None

def Tripoli_read_part(line):
	line = line.split()
	if line[0] == "NEUTRON" or line[0] == "PHOTON": # Linea de texto es particula
		[E,x,y,z,dx,dy,dz,w] = np.double(line[1:])
		part = [E,x,y,z,dx,dy,dz]
		return [part,w]
	return None

def Decay_read(line):
	line = line.split()
	E = np.double(line[0])
	w = np.double(line[2])
	part = [E]
	return [part,w]

def Tripoli_read_part(line):
	line = line.split()
	if line[0] == "NEUTRON" or line[0] == "PHOTON": # Linea de texto es particula
		[E,x,y,z,dx,dy,dz,w] = np.double(line[1:])
		part = [x,y,z]
		return [part,w]
	return None