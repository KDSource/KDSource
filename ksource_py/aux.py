# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


R_gaussian = 1 / (2*np.sqrt(np.pi)) # Roughness of gaussian kernel

def odd_fact(n):
	if n <= 1:
		return 1
	else:
		return n*odd_fact(n-2)

def C_gaussian(q): # Silverman costant for gaussian kernel and dimension q
	return (4/(2+q))**(1/(4+q))

def H10Factor(pt='n'):
	if pt == 'n':
		E,H10 = np.loadtxt("/home/inti/Documents/Maestria/KSource/ksource_py/ARN_neutron", unpack=True)
	elif pt == 'p':
		E,H10 = np.loadtxt("/home/inti/Documents/Maestria/KSource/ksource_py/ARN_photon", unpack=True)
	else:
		raise ValueError("Tipo de particula invalido")
	log_E = np.log(1E-6*E)
	log_H10 = np.log(1E-12*H10)
	lin_interp = interpolate.interp1d(log_E, log_H10)
	log_interp = lambda EE: np.exp(lin_interp(np.log(EE)))
	return log_interp

class BoxMask:
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

def apply_weight_mask(vecs, ws, weightfun=None, maskfun=None):
	if weightfun is not None:
		ws = ws*weightfun(vecs)
	mask = (ws > 0)
	if maskfun is not None:
		mask = np.logical_and(mask, maskfun(vecs))
	return [vecs[mask,:], ws[mask]]

class TracksStats:
	def __init__(self, vecs, ws, weightfun=None, maskfun=None):
		self.vecs, self.ws = apply_weight_mask(vecs, ws, weightfun, maskfun)
		if len(ws) == 0:
			raise Exception("Lista de tracks vacia")
		self.dim = vecs.shape[1]
	def intensity(self, steps=1, weightfun=None, maskfun=None):
		vecs, ws = apply_weight_mask(self.vecs, self.ws, weightfun, maskfun)
		Is = []
		errs = []
		Ns = list(range(0, len(ws), int(len(ws)/steps)))[1:] + [len(ws)]
		for N in Ns:
			Is.append(np.sum(ws[:N])/N)
			errs.append(np.sqrt(np.sum(ws[:N]**2))/N)
		Is = np.array(Is)
		errs = np.array(errs)
		if(steps > 1):
			plt.errorbar(Ns, Is, errs)
			plt.xlabel("Cantidad de muestras")
			plt.ylabel("I / N (peso medio)")
		return [Ns, Is, errs]
	def mean(self, var, steps=1, weightfun=None, maskfun=None):
		vecs, ws = apply_weight_mask(self.vecs, self.ws, weightfun, maskfun)
		vecs = vecs[:,var]
		mns = []
		errs = []
		Ns = list(range(0, len(ws), int(len(ws)/steps)))[1:] + [len(ws)]
		for N in Ns:
			I = np.sum(ws[:N])
			mns.append(np.sum(ws[:N]*vecs[:N])/I)
			errs.append(np.sqrt((np.sum(ws[:N]*vecs[:N]**2)/I - mns[-1]**2)/I))
		mns = np.array(mns)
		errs = np.array(errs)
		if(steps > 1):
			plt.errorbar(Ns, mns, errs)
			plt.xlabel("Cantidad de muestras")
			plt.ylabel("Media")
		return [Ns, mns, errs]
