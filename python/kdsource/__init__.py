#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python module for creating and optimizing KDSource objects.

KDSource's are particle sources for Monte Carlo radiation transport
simulations. The full distribution and documentation can be found in the
project GitHub page:

    https://github.com/KDSource/KDSource

A KDSource object is based on a particle list in MCPL format (see
https://mctools.github.io/mcpl/), on which a Kernel Density Estimation
(KDE) is applied, using KDEpy library (see
https://github.com/tommyod/KDEpy).

With kdsource Python library you can create, optimize, and export
KDSource objects as XML files. These files can be later used as
distributional sources in other Monte Carlo simulations, using the tools
available in the full distribution.

If you use KDSource tools in your work, please add the following
reference:

    Abbate, O. I., Schmidt, N. S., Prieto, Z. M., Robledo, J. I.,
    Dawidowski, J., Márquez, A. A., & Márquez Damián, J. I. KDSource,
    a tool for the generation of Monte Carlo particle sources using
    kernel density estimation [Computer software].
    https://github.com/KDSource/KDSource

"""

__version__ = "0.1.0"
__author__ = "inti-abbate"

from .geom import Geometry, Metric
from .kde import bw_silv
from .kdsource import KDSource, load, run
from .plist import PList, appendssv, convert2mcpl, join2mcpl, savessv
from .surfsource import SurfaceSourceFile, create_source_file
