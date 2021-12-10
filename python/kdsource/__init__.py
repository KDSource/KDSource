#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for creating and optimizing KDSource objects

KDSource's are particle sources for Monte Carlo radiation transport
simulations. The full distribution and documentation can be found in the
project GitHub page:

    https://github.com/inti-abbate/KDSource

A KDSource object is based on a particle list in MCPL format (see 
https://mctools.github.io/mcpl/), on which a Kernel Density Estimation
(KDE) is applied, using KDEpy library (see
https://github.com/tommyod/KDEpy).

With kdsource Python library you can create, optimize, and export KDSource
objects. The exported files can be later used as distributional sources
in other Monte Carlo simulations, using the tools available in the full
distribution.

If you use KDSource tools in your work, please add the following
reference:

    O. I. Abbate, Desarrollo de herramientas computacionales para el
    cálculo de blindajes, Master's thesis, Universidad Nacional de Cuyo,
    2021.

"""

__version__ = "0.1.0"
__author__ = "inti-abbate"

from .utils import *
from .kde import *

from .geom import *
from .plist import *
from .kdsource import *

from .tally import *
from .stats import *
from .summary import *