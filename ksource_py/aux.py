# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


# Parametros para KDE

R_gaussian = 1 / (2*np.sqrt(np.pi)) # Roughness of gaussian kernel

def C_gaussian(q): # Silverman costant for gaussian kernel and dimension q
	return (4/(2+q))**(1/(4+q))


# Factores dosimetricos

KSOURCEPY_PATH = "/home/inti/Documents/Maestria/KSource/ksource_py/"

def H10(pt='n', ref='ARN'):
	if ref == 'ARN':
		if pt == 'n':
			E,H10 = np.loadtxt(KSOURCEPY_PATH+"ARN_neutron", unpack=True)
		elif pt == 'p':
			E,H10 = np.loadtxt(KSOURCEPY_PATH+"ARN_photon", unpack=True)
		else:
			raise ValueError("Tipo de particula invalido")
	elif ref == 'ICRP':
		if pt == 'n':
			E,H10 = np.loadtxt(KSOURCEPY_PATH+"ICRP_neutron", unpack=True)
		elif pt == 'p':
			E,H10 = np.loadtxt(KSOURCEPY_PATH+"ICRP_photon", unpack=True)
		else:
			raise ValueError("Tipo de particula invalido (validos: 'n', 'p')")
	else:
		raise ValueError("Referencia invalida (validas: 'ARN', 'ICRP')")
	log_E = np.log(1E-6*E)
	log_H10 = np.log(1E-12*H10)
	lin_interp = interpolate.interp1d(log_E, log_H10)
	log_interp = lambda EE: np.exp(lin_interp(np.log(EE)))
	return log_interp

# Mascaras

class Box:
	def __init__(self, vec1, vec2):
		if vec1 is None:
			vec1 = len(vec2)*[None]
		if vec2 is None:
			vec2 = len(vec1)*[None]
		if len(vec1) != len(vec2):
			raise ValueError("vec1 y vec2 deben tener la misma longitud")
		self.dim = len(vec1)
		self.vec1 = vec1
		self.vec2 = vec2
	def __call__(self, vecs):
		mask = np.ones(len(vecs), dtype=bool)
		for i in range(self.dim):
			if self.vec1[i] is not None:
				mask = np.logical_and(mask, self.vec1[i] < vecs[:,i])
			if self.vec2[i] is not None:
				mask = np.logical_and(mask, vecs[:,i] < self.vec2[i])
		return mask
