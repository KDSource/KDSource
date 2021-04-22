# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st
import os, subprocess
import mcpl

from .aux import pt2pdg,pdg2pt

MCPLPATH = "/opt/mcpl/bin/"


def convert2mcpl(filename, readformat):
	# Buscar archivo MCPL con mismo nombre
	already_converted = False
	if os.path.isfile(filename.split('.')[0]+".mcpl"):
		already_converted = True
		filename = filename.split('.')[0]+".mcpl"
	if os.path.isfile(filename.split('.')[0]+".mcpl.gz"):
		already_converted = True
		filename = filename.split('.')[0]+".mcpl.gz"
	if readformat == "mcpl":
		already_converted = True
	if already_converted:
		print("Using existing file {}".format(filename))
		return filename
	# Se debe convertir
	if not readformat in ["ssw", "ptrac", "stock", "ssv"]:
		raise Exception("Formato {} invalido".format(readformat))
	print("Converting {} file to MCPL...".format(readformat))
	filename_mcpl = filename.split(".")[0]+".mcpl"
	result = subprocess.run([MCPLPATH+readformat+"2mcpl", filename, filename_mcpl],
							stdout=subprocess.PIPE,
							stderr=subprocess.PIPE,
							check=True) # if result.returncode != 0: raise exception
	stdout = result.stdout.split()
	filename = stdout[stdout.index(b'Created')+1].decode("utf-8")
	print("Done. Created {}".format(filename))
	return filename

ptmap = {'n':2112, 'p':22, 'e':11, 'e+':-11}

def savessv(pt, parts, ws, outfile, comments=None): # Equivalent to convert2ascii (in mcpl.py)
	fout = open(outfile,'w')
	fout.write("#MCPL-ASCII\n#ASCII-FORMAT: v1\n#NPARTICLES: %i\n"%len(parts))
	if comments is not None:
		fout.write("#NCOMMENTS %d\n"%len(comments))
		for comment in comments:
			fout.write(comment)
	fout.write("#END-HEADER\n")
	fout.write("index     pdgcode               ekin[MeV]                   x[cm]          "
			   +"         y[cm]                   z[cm]                      ux                  "
			   +"    uy                      uz                time[ms]                  weight  "
			   +"                 pol-x                   pol-y                   pol-z  userflags\n")
	pdgcode = ptmap[pt]
	fmtstr="%5i %11i %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g 0x%08x\n"
	for idx,p in enumerate(parts):
		fout.write(fmtstr%(idx,pdgcode,p[0],p[1],p[2],p[3],p[4],p[5],p[6],0,ws[idx],0,0,0,0))

def appendssv(pt, parts, ws, outfile, comments=None): # Equivalent to convert2ascii (in mcpl.py)
	fout = open(outfile,'a')
	pdgcode = ptmap[pt]
	fmtstr="%5i %11i %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g %23.18g 0x%08x\n"
	for idx,p in enumerate(parts):
		fout.write(fmtstr%(idx,pdgcode,p[0],p[1],p[2],p[3],p[4],p[5],p[6],0,ws[idx],0,0,0,0))

class PList:
	def __init__(self, filename, readformat="mcpl", pt='n', trasl=None, rot=None, switch_x2z=False, set_params=True):
		filename = convert2mcpl(filename, readformat) # Convertir formato a MCPL
		self.filename = filename
		if trasl is not None:
			trasl = np.array(trasl).reshape(-1)
			if trasl.shape != (3,):
				raise ValueError("trasl invalido")
		if rot is not None:
			if not isinstance(rot, st.Rotation):
				rot = np.array(rot)
				if rot.shape == (4,):
					rot = st.Rotation.from_quat(rot)
				elif rot.shape == (3,3):
					rot = st.Rotation.from_matrix(rot)
				elif rot.shape == (3,):
					rot = st.Rotation.from_rotvec(rot)
				else:
					raise ValueError("rot invalido")
		self.pt = pt
		self.trasl = trasl
		self.rot = rot
		self.x2z = switch_x2z
		self.N = mcpl.MCPLFile(self.filename).nparticles
		self.params_set = False
		if set_params:
			self.set_params()
		else:
			self.I = self.p2 = 1.0

	def set_params(self):
		pl = mcpl.MCPLFile(self.filename)
		I = p2 = 0
		for pb in pl.particle_blocks:
			I += pb.weight.sum()
			p2 += (pb.weight**2).sum()
		print("I = {}\np2 = {}\nN = {}".format(I, p2, self.N))
		self.I = I
		self.p2 = p2
		self.params_set = True

	def get(self, N=-1, skip=0):
		if N < 0:
			if not self.params_set: self.set_params()
			N = self.N
		pl = mcpl.MCPLFile(self.filename, blocklength=N)
		pl.skip_forward(skip)
		pb = pl.read_block()
		if pb is None:
			return np.zeros(0),np.zeros((0,7))
		parts = np.stack((pb.ekin,pb.x,pb.y,pb.z,pb.ux,pb.uy,pb.uz), axis=1)
		ws = pb.weight
		mask = np.logical_and(ws>0, pb.pdgcode==pt2pdg(self.pt))
		parts = parts[mask]
		ws = ws[mask]
		if self.trasl is not None: # Aplico traslacion
			parts[:,1:4] += self.trasl
		if self.rot is not None: # Aplico rotacion
			parts[:,1:4] = self.rot.apply(parts[:,1:4])
			parts[:,4:7] = self.rot.apply(parts[:,4:7])
		if self.x2z: # Aplico permutacion (x,y,z) -> (y,z,x)
			ekin,x,y,z,ux,uy,uz = parts.T
			parts = np.array([ekin,y,z,x,uy,uz,ux]).T
		return parts,ws

	def save(self, file):
		file.write(self.pt+'\n')
		file.write(os.path.abspath(self.filename)+'\n')
		if self.trasl is not None: np.savetxt(file, self.trasl[np.newaxis,:])
		else: file.write('\n')
		if self.rot is not None: np.savetxt(file, self.rot.as_rotvec()[np.newaxis,:])
		else: file.write('\n')
		file.write("%d\n" % (self.x2z))