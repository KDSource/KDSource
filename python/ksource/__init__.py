#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for creating and optimizing KSource objects

KSource's are particle sources for Monte Carlo radiation transport
simulations. The full distribution and documentation can be found in the
project GitHub page:

    https://github.com/inti-abbate/KSource

A KSource object is based on a particle list in MCPL format (see 
https://mctools.github.io/mcpl/), on which a Kernel Density Estimation
(KDE) is applied, using KDEpy library (see
https://github.com/tommyod/KDEpy).

With ksource Python library you can create, optimize, and export KSource
objects. The exported files can be later used as distributional sources
in other Monte Carlo simulations, using the tools available in the full
distribution.

If you use KSource tools in your work, please add the following
reference:

    [TODO: Completar referencia]

"""

__version__ = "1.0.0"
__author__ = "inti-abbate"

from .utils import *
from .kde import *

from .geom import *
from .plist import *
from .ksource import *

from .tally import *
from .stats import *
from .summary import *