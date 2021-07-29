# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def apply_weight_mask(vecs, ws, weightfun=None, maskfun=None):
	if weightfun is not None:
		ws = ws*weightfun(vecs)
	mask = (ws > 0)
	if maskfun is not None:
		mask = np.logical_and(mask, maskfun(vecs))
	return [vecs[mask,:], ws[mask]]

def convergence(vecs, ws, param, fracmin=0.1, steps=10, plot=True):
	Nmax = len(ws)
	Nmin = int(fracmin * Nmax)
	Ns = np.rint(np.linspace(Nmin, Nmax, num=steps+1)).astype(int)
	Ns = Ns[Ns>0]
	params = []
	errs = []
	for N in Ns:
		param_err = param(vecs[:N],ws[:N])
		params.append(param_err[0])
		errs.append(param_err[1])
	params = np.array(params)
	errs = np.array(errs)
	if Ns[0]!=Ns[-1] and plot:
		plt.plot(Ns, params, 'o-')
		plt.fill_between(Ns, params-errs, params+errs, color='blue', alpha=0.3)
		plt.xlabel("Cantidad de muestras")
	return [Ns, params, errs]

def mean_weight(vecs, ws):
	return [ws.mean(), ws.std()/np.sqrt(len(ws))]
def mean(vecs, ws, var):
	mean = np.average(vecs[:,var], weights=ws)
	err = np.sqrt(np.average((vecs[:,var]-mean)**2, weights=ws)/np.sum(ws))
	return [mean, err]
def std(vecs, ws, var):
	mean = np.average(vecs[:,var], weights=ws)
	s2 = np.average((vecs[:,var]-mean)**2, weights=ws)
	std = np.sqrt(s2)
	varerr = np.sqrt(np.average(((vecs[:,var]-mean)**2-s2)**2, weights=ws)/np.sum(ws))
	stderr = varerr * 0.5 / std
	return [std, stderr]

class Stats:
	def __init__(self, vecs, ws, weightfun=None, maskfun=None):
		if len(vecs) != len(ws):
			raise ValueError("Longitud de vecs y ws debe ser igual")
		self.vecs, self.ws = apply_weight_mask(vecs, ws, weightfun, maskfun)
		if len(ws) == 0:
			raise Exception("Lista de tracks vacia")
		self.dim = vecs.shape[1]
		self.N = len(ws)
	def mean_weight(self, fracmin=0.1, steps=10, plot=True):
		Ns,params,errs = convergence(self.vecs, self.ws, mean_weight, fracmin=fracmin, steps=steps, plot=plot)
		if(plot):
			plt.ylabel("Peso medio")
		return [Ns, params, errs]
	def mean(self, var, varname=None, fracmin=0.1, steps=10, plot=True):
		mean_var = lambda vecs,ws: mean(vecs,ws,var)
		Ns,params,errs = convergence(self.vecs, self.ws, mean_var, fracmin=fracmin, steps=steps, plot=plot)
		if(plot):
			if varname is None:
				varname = "v%d"%var
			plt.ylabel("Valor medio de variable {}".format(varname))
		return [Ns, params, errs]
	def std(self, var, varname=None, fracmin=0.1, steps=10, plot=True):
		std_var = lambda vecs,ws: std(vecs,ws,var)
		Ns,params,errs = convergence(self.vecs, self.ws, std_var, fracmin=fracmin, steps=steps, plot=plot)
		if(plot):
			if varname is None:
				varname = "v%d"%var
			plt.ylabel("Desv. estandar de variable {}".format(varname))
		return [Ns, params, errs]
