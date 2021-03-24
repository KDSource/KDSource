# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class Tracks:
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