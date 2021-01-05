# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV, LeaveOneOut

class KSource:
	def __init__(self, metric, plist, N=1e4, bw=0):
		self.metric = metric
		self.plist = plist
		self.kde = KernelDensity(bandwidth=bw)
		
	def optimize_bw(self, parts):
		C_silv = 0.9397
		vecs = self.metric.transform(parts)
		std = self.metric.std(vecs, transform=False)
		bw_silv = C_silv * std * len(parts)**(-1/(4+self.metric.dim)) # Regla del dedo de Silverman
		#
		nsteps = 20 # Cantidad de pasos para bw
		max_fact = 1.5 # Rango de bw entre bw_silv/max_fact y bw_silv*max_fact
		max_log = np.log10(max_fact)
		bw_grid = bw_silv * np.logspace(-max_log,max_log,nsteps).reshape(-1,1) # Grilla de bandwidths
		#
		cv = 10
		grid = GridSearchCV(KernelDensity(kernel='gaussian'),
		                                  {'bandwidth': bw_grid},
		                                  cv=cv,
		                                  verbose=10,
		                                  n_jobs=8)
		grid.fit(vecs)
		bw_best = grid.best_params_['bandwidth']
		#
		self.kde = KernelDensity(bandwidth=bw_best)

	def score(self, parts):
		vecs = self.metric.transform(parts)
		jacs = self.metric.jac(parts)
		scores = self.kde.eval_scores(vecs)
		return jacs * scores