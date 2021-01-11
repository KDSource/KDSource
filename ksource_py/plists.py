# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class PList:
	def __init__(self, readfun, filename):
		self.filename = filename
		self.read = readfun
		self.set_params()
	def set_params(self):
		file = open(self.filename, "r")
		I = p2 = N = 0
		for line in file:
			part_w = self.read(line)
			if part_w is not None: # Linea de texto es particula
				w = part_w[1]
				I += w
				p2 += w*w
				N += 1
		file.close()
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
				if cont >= skip: 
					parts.append(part_w[0])
					ws.append(part_w[1])
				cont += 1
				if N>0 and cont==skip+N:
					break
		file.close()
		return [np.array(parts), np.array(ws)]

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