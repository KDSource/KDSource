<p align="center">
<img src="img/kdsource_logo.png" alt="logo" height="100"/>
</p>

![KDSource](https://img.shields.io/badge/KDSource-0.1.0-brightgreen)
[![Documentation Status](https://readthedocs.org/projects/kdsource/badge/?version=latest)](https://kdsource.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-GPL3-brightgreen.svg)](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
![CMake](https://img.shields.io/badge/Cmake-3.0+-orange)
[![Python 3.10](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/release/python-310/)


This is source version of KDSource, a tool for Monte Carlo particle sources generation
using Kernel Density Estimation. Visit our [Documentation Page](https://kdsource.readthedocs.io/en/latest/) for more details!

KDSource assists Monte Carlo beams and shielding calculations, improving tally results in difficult problems. It allows to model big systems (e.g.: investigation reactor guides hall) thru spatial or temporal coupling of different simulations in different transport codes, implementing as well variance reduction.

It processes particle lists recorded as output of a simulation (e.g.: passing thru a window), to be used as input in another one. It estimates density distribution in energy, position and direction by means of Kernel Density Estimation (KDE) technique, allowing visualizing it as well as using it to produce new particles (artificial, but with the same estimated density). This allows to increase the number of source particles in the second simulation, improving its statistics (variance reduction).

KDSource uses [`MCPL`](https://mctools.github.io/mcpl/) particle lists format.

In TRIPOLI-4 and McStas it is possible "on-the-fly" sampling during simulations, while for other codess it is necessary to record the source particle list as MCPL files before the simulation.


## Contents:

The KDSource package consists in the following tools:

* Python API: Allows creating, optimizing, analyzing, plotting, and saving KDE sources. Optimización consists in automatic selection of bandwidth. Internally, it uses [`KDEpy`](https://kdepy.readthedocs.io/en/latest/) library for KDE.

* C API: Allows loading the sources saved with Python, and generating new synthetic samples. These follow the estimated distribution, and can be saved in a new MCPL file or be introduced directly in a simulation.

* Command line API: Allows easily producing samples, based on sources saved with Python. Also allows to access templates and communication files, as well as MCPL applications.


## Installation:

This is currently being updated with the intention of providing conda and PyPI packages. For the time being, note that you can install via pip if your machine has compilers and libxml2 installed:

```
pip install "git+https://github.com/KDSource/KDSource"
```

Note that on Linux, libxml2 can most likely be installed from your package manager. For instance on Ubuntu:

```bash
   $ sudo apt-get update
   $ sudo apt-get install libxml2-dev
```

On macOS, libxml2 is most likely part of the standard SDKs.

After installation, you should be able to `import kdsource.api as kds` in Python, and have command-line tools like `kdsource-config` and `kdsource-resample` available.


## Usage examples and templates

See the [documentation](https://kdsource.readthedocs.io/en/latest/) page for usage instructions, tutorials, and a detailed documentation of all the functionalities in KDSource.

Usage examples can be found in the [`docs/examples`](docs/examples) subdirectory. At the moment these are:
* [`Verification.ipynb`](docs/examples/Verification.ipynb): Analytic example. KDSource is used to generate a source from a particle list sampled from an known correlated distribution, and the generated particles distributions are compared with the analytical density.

## Issues and contributing

If you are having trouble using the package, please let us know by creating a [Issue on GitHub](https://github.com/KDSource/KDSource/issues) and we'll get back to you.

Contributions are very welcome. To contribute, fork the project, create a branch and submit a Pull Request.


## Reference

Usage of the KDSource package is allowed in the terms detailed in the LICENSE file. However, if you use it for your work, we
would appreciate it if you would use the following reference:

Abbate, O. I., Schmidt, N. S., Prieto, Z. M., Robledo, J. I., Dawidowski, J., Márquez, A. A., & Márquez Damián, J. I. KDSource, a tool for the generation of Monte Carlo particle sources using kernel density estimation [Computer software]. https://github.com/KDSource/KDSource

For how to reference MCPL, refer to: https://mctools.github.io/mcpl/about/
