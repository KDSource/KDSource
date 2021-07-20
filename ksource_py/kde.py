# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)

from sklearn.neighbors import NearestNeighbors
from sklearn.model_selection import KFold
from KDEpy import TreeKDE
from joblib import Parallel, delayed


def bw_silv(dim, N_eff):
	return (4/(2+dim))**(1/(4+dim)) * N_eff**(-1/(4+dim))

def bw_knn(data, weights=None, K_eff=100, k=None, batch_size=10000):
	N,dim = data.shape
	batches = int(max(1, np.round(N / batch_size)))
	N_eff = N if weights is not None else np.sum(weights)**2 / np.sum(weights**2)
	if k is None: # Calculo k en base a K_eff
		k_float = batch_size * K_eff / N_eff
		k = int(max(1, round(k_float)))
		f_k = k_float / k # Factor de correccion para k
	else: # k fijado por usuario
		f_k = 1.0
		K_eff = k * N_eff / batch_size
	print("Usando k = {} vecinos por batch (batch_size = {})".format(k, batch_size))
	print("Factor de correccion: f_k = k_float / k = {}".format(f_k))
	print("Vecinos efectivos totales: K_eff = {}".format(K_eff))
	bw_knn = np.array([])
	for batch,batch_data in enumerate(np.array_split(data, batches)):
		print("batch =", batch+1, "/", batches)
		knn = NearestNeighbors(n_neighbors=k+1, n_jobs=-1) # Sumo 1 para descontar propio punto
		knn.fit(batch_data)
		ds,idxs = knn.kneighbors(batch_data)
		bws = ds[:,-1] * (f_k)**(1/dim) # Selecciono k-esima columna y aplico factor de correccion
		bw_knn = np.concatenate((bw_knn, bws))
	return bw_knn

def _kde_cv_score(bw, data, weights=None, cv=10, modelKDE=TreeKDE):
	folds = KFold(n_splits=cv)
	scores = []
	for train_idx, test_idx in folds.split(data):
		data_train, data_test = data[train_idx], data[test_idx]
		if weights is not None:
			weights_train, weights_test = weights[train_idx], weights[test_idx]
		else:
			weights_train = weights_test = None
		if np.ndim(bw) == 1:
			if len(bw) != len(data):
				raise ValueError("Si bw es lista, debe tener la misma longitud que data")
			bw_train = bw[train_idx]
		else:
			bw_train = bw
		kde = modelKDE(bw=bw_train)
		kde.fit(data_train, weights=weights_train)
		if weights_test is not None:
			score = np.mean(weights_test * np.log(kde.evaluate(data_test)))
		else:
			score = np.mean(np.log(kde.evaluate(data_test)))
		scores.append(score)
	return np.mean(scores)

def bw_mlcv(data, weights=None, cv=10, seed=None, grid=None, show=True):
	if seed is None:
		N_eff = len(data) if weights is None else np.sum(weights)**2/np.sum(weights**2)
		seed = bw_silv(np.shape(data)[1], N_eff)
	if grid is None:
		grid = np.logspace(-1, 1, 20)
	bw_grid = np.reshape(grid, (-1,*np.ones(np.ndim(seed),dtype=int))) * seed
	if cv > len(data):
		cv = len(data)

	cv_scores = Parallel(n_jobs=-1, verbose=10)(
		delayed(_kde_cv_score)(bw, data, weights=weights, cv=cv) for bw in bw_grid)
	idx_best = np.argmax(cv_scores)
	
	plt.plot(grid, cv_scores)
	plt.xlabel("Ancho de banda normalizado")
	plt.ylabel("Mean Test Score")
	plt.tight_layout()
	if(show): plt.show()
	if idx_best in (0, len(bw_grid)-1):
		raise Exception("No se encotro maximo en el rango de bw seleccionado. Mueva la grilla e intente nuevamente.")
	bw_mlcv = bw_grid[idx_best]
	return bw_mlcv

def optimize_bw(bw_method, vecs=None, ws=None, weightfun=None, maskfun=None, **kwargs):
	vecs = np.array(vecs)
	if ws is not None:
		ws = np.array(ws)
		mask = (ws > 0)
	else:
		mask = np.ones(len(vecs), dtype=bool)
	if weightfun is not None: # Aplico weightfun
		ws = ws * weightfun(vecs)
	if maskfun is not None: # Aplico maskfun
		mask = np.logical_and(mask, maskfun(vecs))
	vecs = vecs[mask]
	ws = ws[mask]
	#
	if bw_method == 'silv': # Metodo de Silverman
		dim = np.shape(vecs)[1]
		N_eff = len(vecs) if ws is None else np.sum(ws)**2/np.sum(ws**2)
		return bw_silv(dim, N_eff)
	elif bw_method == 'knn': # Metodo K Nearest Neighbours
		return bw_knn(vecs, weights=ws, **kwargs)
	elif bw_method == 'mlcv': # Metodo Maximum Likelihood Cross Validation
		return bw_mlcv(vecs, weights=ws, **kwargs)
	else:
		raise Exception("bw_method invalido. Validos: 'silv', 'mlcv', 'knn'")
