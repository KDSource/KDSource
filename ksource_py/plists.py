# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.transform as st
import os, subprocess
import mcpl

MCPLPATH = "/opt/mcpl/bin/"


def convert2mcpl(filename, readformat):
	if not isinstance(readformat, ("ssw", "ptrac", "stock")):
	    raise Exception("Formato {} invalido".format(readformat))
	print("Converting tracks file to MCPL...")
	filename_mcpl = filename.split(".")[0]+".mcpl"
	result = subprocess.run([MCPLPATH+readformat+"2mcpl", filename, filename_mcpl],
		                    stdout=subprocess.PIPE,
		                    stderr=subprocess.PIPE,
		                    check=True) # if result.returncode != 0: raise exception
	filename = result.stdout.split()[-1].decode("utf-8")
	print("Done. Created {}".format(filename))
	return filename

class PList:
	def __init__(self, readformat, filename, pt='n', trasl=None, rot=None, switch_x2z=False, set_params=True):
		if readformat != "mcpl": # Convertir formato a MCPL
			filename = convert2mcpl(filename, readformat)
		self.filename = filename
		if trasl is not None:
			trasl = np.array(trasl)
		if rot is not None:
			rot = st.Rotation.from_rotvec(rot)
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
		I = p2 = N = 0
		for pb in pl.particle_blocks:
			N += len(pb.weight)
			I += pb.weight.sum()
			p2 += pb.weight**2.sum()
		print("I = {}\np2 = {}\nN = {}\nnpart = {}".format(I, p2, N, self.N))
		self.I = I
		self.p2 = p2
		self.params_set = True

	def get(self, N=-1, skip=0):
		if N < 0: N = self.N
		pl = mcpl.MCPLFile(self.filename, blocklength=N)
		pl.skip_forward(skip)
		pb = pl.read_block()
		if pb is not None:
			if self.trasl is not None: # Aplico traslacion
				poss = np.array([pb.x,pb.y,pb.z]).T
				poss += self.trasl
				pb.x,pb.y,pb.z = poss.T
			if self.rot is not None: # Aplico rotacion
				poss = np.array([pb.x,pb.y,pb.z]).T
				poss = self.rot.apply(poss)
				pb.x,pb.y,pb.z = poss.T
				dirs = np.array([pb.ux,pb.uy,pb.uz]).T
				dirs = self.rot.apply(dirs)
				pb.ux,pb.uy,pb.uz = dirs.T
			if self.x2z: # Aplico permutacion (x,y,z) -> (y,z,x)
				x,y,z = pb.x,pb.y,pb.z
				pb.x,pb.y,pb.z = y,z,x
				ux,uy,uz = pb.ux,pb.uy,pb.uz
				pb.ux,pb.uy,pb.uz = uy,uz,ux
		return pb

	def save(self, file):
		file.write(self.pt+'\n')
		file.write("{}\n".format(self.ord))
		for readformat in self.readformats: file.write(readformat+'\n')
		for filename in self.filenames: file.write(os.path.abspath(filename)+'\n')
		if self.trasl is not None: np.savetxt(file, self.trasl[np.newaxis,:])
		else: file.write('\n')
		if self.rot is not None: np.savetxt(file, self.rot.as_rotvec()[np.newaxis,:])
		else: file.write('\n')
		file.write("%d\n" % (self.x2z))