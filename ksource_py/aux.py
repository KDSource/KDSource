# -*- coding: utf-8 -*-

import numpy as np


R_gaussian = 1 / (2*np.sqrt(np.pi)) # Roughness of gaussian kernel

def odd_fact(n):
	if n <= 1:
		return 1
	else:
		return n*odd_fact(n-2)

def C_gaussian(q): # Silverman costant for gaussian kernel and dimension q
	return (4/(2+q))**(1/(4+q))