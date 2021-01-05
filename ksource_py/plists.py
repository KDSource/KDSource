# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class PList:
	def __init__(self, readfun, filename):
		self.filename = filename
		self.read = readfun
	def get(self, N, skip=0):
		file = open(filename, "r")
		parts = []
		cont = 0
		for line in file:
			if part=self.read(line) is not None: # Linea de texto es particula
				if cont >= skip: 
					parts.append(part)
				cont += 1
			if cont == skip+N:
				break
		return np.array(parts)

def PTRAC_read(line):
	line = line.split()
	if len(line) == 9: # Linea de texto es particula
		[x,y,z,dx,dy,dz,E,w,t] = np.double(line)
		part = [E,x,y,z,dx,dy,dz,w]
		return part
	return None

def Tripoli_read_part(line):
	line = line.split()
	if line[0] == "NEUTRON" or line[0] == "PHOTON": # Linea de texto es particula
		[E,x,y,z,dx,dy,dz,w] = np.double(line[1:])
		part = [E,x,y,z,dx,dy,dz,w]
		return part
	return None

def Decay_read(line):
	line = line.split()
	E = np.double(line[0])
	w = np.double(line[2])
	part = [E,w]
	return part

def Tripoli_read_part(line):
	line = line.split()
	if line[0] == "NEUTRON" or line[0] == "PHOTON": # Linea de texto es particula
		[E,x,y,z,dx,dy,dz,w] = np.double(line[1:])
		part = [x,y,z,w]
		return part
	return None