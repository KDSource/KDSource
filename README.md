<p align="center">
<img src="img/kdsource_logo.png" alt="logo" height="100"/>
</p>

![KDSource](https://img.shields.io/badge/KDSource-0.1.0-brightgreen)
[![Documentation Status](https://readthedocs.org/projects/kdsource/badge/?version=latest)](https://kdsource.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-GPL3-brightgreen.svg)](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
![CMake](https://img.shields.io/badge/Cmake-3.0+-orange)
[![Python 3.6](https://img.shields.io/badge/python-3.6%2B-blue)](https://www.python.org/downloads/release/python-360/)



This is source version of KDSource, a tool for Monte Carlo particle sources generation
using Kernel Density Estimation.

KDSource assists Monte Carlo beams and shielding calculations. It allows to model big systems (e.g.: investigation reactor guides hall) thru spatial or temporal coupling of different simulations in different transport codes, implementing as well variance reduction.

It processes particle lists recorded as output of a simulation (e.g.: passing thru a window), to be used as input in another one. It estimates density distribution in energy, position and direction by means of Kernel Density Estimation (KDE) technique, allowing visualizing it as well as using it to produce new particles (artificial, but with the same estimated density). This allows to increase the number of source particles in the second simulation, improving its statistics (variance reduction).

KDSource uses [`MCPL`](https://mctools.github.io/mcpl/) particle lists format. In its modified version included in this distribution, it allows working with the following Monte Carlo codes:
* MCNP
* PHITS
* McStas
* TRIPOLI-4

In TRIPOLI-4 and McStas it is possible "on-the-fly" sampling during simulations, while for the other formats it is necessary to record the source particle list before the simulation.


## Contents:

The KDSource package consists in the following tools:

* Python API: Allows creating, optimizing, analyzing, plotting, and saving KDE sources. Optimización consists in automatic selection of bandwidth. Internally, it uses [`KDEpy`](https://kdepy.readthedocs.io/en/latest/) library for KDE.

* C API: Allows loading the sources saved with Python, and generating new synthetic samples. These follow the estimated distribution, and can be saved in a new MCPL file or be introduced directly in a simulation.

* Templates and communication files for Monte Carlo codes. Specific files are included for utilization of KDSource tools in McStas y TRIPOLI-4 simulations.

* Command line API: Allows easily producing samples, based on sources saved with Python. Also allows to access templates and communication files, as well as MCPL applications.



## Installation (Linux):

The installation instructions assume this repo has been cloned to a local directory.
	
1. Go to source directory:

   ```bash
   $ cd /path/to/kdsourcesource
   ```

   Where `/path/to/kdsourcesource` is the folder where the source distribution of KDSource was extracted.

2. Install with `cmake`:

   ```bash
   $ mkdir build && cd build
   $ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/kdsourceinstall
   $ make install
   $ cd ..
	```
   Where `path/to/kdsourceinstall` is the folder where you wish to install KDSource internal files.

   It is required to have previously installed `libxml2`.

3. Install Python API with `pip`:

   ```bash
   $ cd python
   $ pip install .
   $ cd ..
   ```

4. KDSource is ready to be used in `/path/to/kdsourceinstall`. For example, you can see the `kdtool` command options with:

   ```bash
   $ /path/to/kdsourceinstall/bin/kdtool --help
   ```

   If you wish to have KDSource tools available in your path, execute:

   ```bash
   $ export PATH=$PATH:/path/to/kdsourceinstall/bin
   ```
   Or add this command to `~/.profile` (and update with `source ~/.profile`).


## Usage examples and templates

Usage examples can be found in the [`docs/examples`](docs/examples) subdirectory. At the moment these are:
* [`Verification.ipynb`](docs/examples/Verification.ipynb): Analytic example. KDSource is used to generate a source from a particle list sampled from an known correlated distribution, and the generated particles distributions are compared with the analytical density.

Moreover, templates for common usage of KDSource in Monte Carlo simulations can be found in the [`templates`](templates) subdirectory, and can be copied to the working directory via the `kdtool templates .` command. They are:
* [`preproc_tracks.ipynb`](templates/preproc_tracks.ipynb): Template for the generation of a KDE source from a particle list registered with any of the [`MCPL`](https://mctools.github.io/mcpl/)-compatible Monte Carlo codes. The generated source can be used as input of any of said codes, generating an unlimited number of particles.
* [`preproc_tally.ipynb`](templates/preproc_tally.ipynb): Template for the generation of a volumetric KDE source from a TRIPOLI-4 reaction tally (usually activation). The generated source can be used as input of any of the [`MCPL`](https://mctools.github.io/mcpl/)-compatible Monte Carlo codes, generating an unlimited number of particles.
* [`postproc.ipynb`](templates/postproc.ipynb): Template for collecting integral results of simulations with McStas and/or TRIPOLI-4.
* [`doseplots.ipynb`](templates/doseplots.ipynb): Template for plotting TRIPOLI-4 volume tallies (usually dose maps).
* McStas templates:
  * [`exe_McStas.sh`](templates/mcstas/exe_McStas.sh): Template for executing McStas with KDSource.
* TRIPOLI-4 templates:
  * [`exe_Tripoli.sh`](templates/tripoli/exe_McStas.sh): Template for executing TRIPOLI-4 with KDSource.
  * [`KDSource.c`](templates/tripoli/KDSource.c): Template for using KDSource as an external source.
  * [`template.t4`](templates/tripoli/template.t4): Template for a TRIPOLI-4 input.


## Issues and contributing

If you are having trouble using the package, please let us know by creating a [Issue on GitHub](https://github.com/inti-abbate/KDSource/issues) and we'll get back to you.

Contributions are very welcome. To contribute, fork the project, create a branch and submit a Pull Request.


## Reference

Usage of the KDSource package is allowed in the terms detailed in the LICENSE file. However, if you use it for your work, we
would appreciate it if you would use the following reference:

Abbate, O. I., Schmidt, N. S., Prieto, Z. M., Robledo, J. I., Dawidowski, J., Márquez, A. A., & Márquez Damián, J. I. KDSource, a tool for the generation of Monte Carlo particle sources using kernel density estimation [Computer software]. https://github.com/KDSource/KDSource