#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for utility functions
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


# Conversion entre tipo de part (char) y pdgcode
def pt2pdg(pt):
    """
    Convert particle type to pdgcode.

    Valid particle types are:
        - 'n': neutron  (pdgcode = 2112)
        - 'p': photon   (pdgcode = 22)
        - 'e': electron (pdgcode = 11)
    For other particle types, 0 is returned.

    Parameters
    ----------
    pt: str
        Particle type.
    """
    if pt == 'n': return 2112
    if pt == 'p': return 22
    if pt == 'e': return 11
    return 0
def pdg2pt(pdgcode):
    """
    Convert pdgcode to particle type.

    Available pdgcode to be converted are:
        - 2112: neutron  (pt = 'n')
        - 22:   photon   (pt = 'p')
        - 11:   electron (pt = 'e')
    For other pdgcodes, '0' is returned.

    Parameters
    ----------
    pdgcode: int
        Particle PDG code.
    """
    if pdgcode == 2112: return 'n'
    if pdgcode == 22: return 'p'
    if pdgcode == 11: return 'e'
    return '0'


# Dosimetric factors

fact_dos = os.path.dirname(__file__)+"/fact_dos"

def H10(pt='n', ref='ICRP'):
    """
    Object for loading and interpolating dosimetric factors.

    Units are MeV for energy and pSv*cm^2 for dosimetric factors.

    Parameters
    ----------
    pt: str
        Particle type. Available are 'n' (neutron) and 'p' (photon).
    ref: str
        Reference for dosimetric factors. Available are:
        - 'ICRP': International Commission on Radiological Protection.
        - 'ARN': Nuclear Regulatory Authority of Argentina.
    """
    if ref == 'ARN':
        if pt == 'n':
            E,H10 = np.loadtxt(fact_dos+"/ARN_neutron", unpack=True)
        elif pt == 'p':
            E,H10 = np.loadtxt(fact_dos+"/ARN_photon", unpack=True)
        else:
            raise ValueError("Tipo de particula invalido")
    elif ref == 'ICRP':
        if pt == 'n':
            E,H10 = np.loadtxt(fact_dos+"/ICRP_neutron", unpack=True)
        elif pt == 'p':
            E,H10 = np.loadtxt(fact_dos+"/ICRP_photon", unpack=True)
        else:
            raise ValueError("Invalid particle type.")
    else:
        raise ValueError("Invalid reference.")
    with np.errstate(divide = 'ignore'): # Avoid warning
        log_E = np.log(1E-6*E)
        log_H10 = np.log(H10)
    lin_interp = interpolate.interp1d(log_E, log_H10)
    log_interp = lambda EE: np.exp(lin_interp(np.log(EE)))
    return log_interp

# Masks

class Box:
    def __init__(self, vec0, vec1):
        """
        Functor class for box-style masks.

        Parameters
        ----------
        vec0: array-like
            Lower limit for each variable.
        vec0: array-like
            Upper limit for each variable.
        """
        if vec0 is None and vec1 is None:
            self.dim = 0
        else:
            if vec0 is None:
                vec0 = len(vec1)*[None]
            if vec1 is None:
                vec1 = len(vec0)*[None]
            if len(vec0) != len(vec1):
                raise ValueError("vec0 and vec1 must have same len.")
            self.dim = len(vec0)
        self.vec0 = vec0
        self.vec1 = vec1
    def __call__(self, vecs):
        """
        Box-style masking function.

        Parameters
        ----------
        vecs: array-like
            Array of points to evaluate mask.

        Returns
        -------
        mask: array-like
            Array of bools with same len as vecs. Each element has the
            following value:
            - True: if the corresponding vecs element has all its
                    variables above vec0 values and below vec1 values.
            - False: otherwise.
        """
        mask = np.ones(len(vecs), dtype=bool)
        for i in range(self.dim):
            if self.vec0[i] is not None:
                mask = np.logical_and(mask, self.vec0[i] < vecs[:,i])
            if self.vec1[i] is not None:
                mask = np.logical_and(mask, vecs[:,i] < self.vec1[i])
        return mask
