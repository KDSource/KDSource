# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st

class PList:
	def __init__(self, readformat, filename, pt='n', trasl=None, rot=None, switch_x2z=False, set_params=True):
		if np.isscalar(readformat): readformat = [readformat]
		self.readformats = readformat
		self.readfuns = []
		for rformat in readformat: exec("self.readfuns.append("+rformat+"_read)")
		if np.isscalar(filename): filename = [filename]
		self.filenames = filename
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
		I = p2 = N = 0
		for i in range(len(self.filenames)):
			if self.filenames[i] is not None:
				file = open(self.filenames[i], "r")
				I_ = p2_ = N_ = 0
				print("Reading file", self.filenames[i], "...")
				for line in file:
					part_w = self.readfuns[i](line)
					if part_w is not None: # Linea de texto es particula
						w = part_w[1]
						I_ += w
						p2_ += w*w
						N_ += 1
				file.close()
				print("Done")
				if N_ > N: I=I_; p2=p2_; N=N_
		print("I = {}\np2 = {}\nN = {}".format(I, p2, N))
		self.I = I
		self.p2 = p2
		self.N = N
		self.params_set = True

	def get(self, N=-1, skip=0):
		partss = []
		wss = []
		for i in range(len(self.filenames)):
			parts = []
			ws = []
			if self.filenames[i] is not None:
				file = open(self.filenames[i], "r")
				cont = 0
				if skip > 0:
					for line in file:
						part_w = self.readfuns[i](line)
						if part_w is not None: # Linea de texto es particula
							cont += 1
							if cont == skip: break
				cont = 0
				for line in file:
					part_w = self.readfuns[i](line)
					if part_w is not None: # Linea de texto es particula
						part,w = part_w
						parts.append(part)
						ws.append(w)
						cont += 1
						if cont==N: break
				file.close()
			else:
				for _ in range(self.N):
					part_w = self.readfuns[i]()
					if part_w is not None: # Linea de texto es particula
						part,w = part_w
						parts.append(part)
						ws.append(w)
			if len(ws) == 0:
				raise Exception("No se pudo obtener particulas de archivo {}".format(self.filenames[i]))
			partss.append(np.array(parts))
			wss.append(np.array(ws))
		N = np.max([len(ws) for ws in wss])
		for i in range(len(self.filenames)):
			n = int(np.ceil(N/len(wss[i])))
			partss[i] = np.tile(partss[i], [n,1])[:N]
			wss[i] = np.tile(wss[i], n)[:N]
		parts = np.concatenate(partss, axis=1)
		ws = np.prod(wss, axis=0)
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
		file.write(self.pt+'\n')
		file.write(self.readformats+'\n')
		file.write(self.filenames+'\n')
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
	try:
		line = line.split(sep=',')
		E = np.double(line[0])/1000.
		w = np.double(line[2])
		E = np.array([E])
		return [E,w]
	except:
		return None

def SSVtally_read(line):
	line = line.split()
	if len(line) == 4: # Linea de texto es particula
		[x,y,z,w] = np.double(line)
		pos = [x,y,z]
		return [pos,w]
	return None

def Isotrop_read(line=None):
	dz = np.random.uniform(-1, 1);
	dxy = np.sqrt(1-dz**2);
	phi = np.random.uniform(0, 2.*np.pi);
	dx = dxy*np.cos(phi);
	dy = dxy*np.sin(phi);
	return [[dx,dy,dz], 1]