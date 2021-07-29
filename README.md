# KSource

This is source version of KSource package for Monte Carlo radiation calculations.

KSource gives Monte Carlo shielding calculation assistance tools. It allows to model big systems (e.g.: investigation reactor guides hall) thru spatial or temporal coupling of different simulations in different transport codes, implementing as well variance reduction.

It process particle lists recorded as output of a simulation (e.g.: passing thru a window), to be used as input in another one. It estimates density distribution in energy, position and direction by means of Kernel Density Estimation (KDE) technique, allowing visualizing it as well as using it to produce new particles (artificial, but with the same estimated density). This allows to increase the number of source particles in the second simulation, improving its statistics (variance reduction).

KSource uses [`MCPL`](https://mctools.github.io/mcpl/) particle lists format. In its modified version included in this distribution, it allows working with the following Monte Carlo codes:
*	MCNP
*	PHITS
*	McStas
*	TRIPOLI-4

In TRIPOLI y McStas it is possible "on-the-fly" sampling during simulations, while for the other formats it is necessary to record the source particle list before the simulation.


## Contents:

The KSource package consists in the following tools:

*	Python API: Allows creating, optimizing, analyzing, plotting, and saving KDE sources. Optimización consists in automatic selection of bandwidth. Internally, it uses [`KDEpy`](https://kdepy.readthedocs.io/en/latest/) library for KDE.

*	C API: Allows loading the sources saved with Python, and generating new synthetic samples. These follow the estimated distribution, and can be saved in a new MCPL file or be introduced directly in a simulation.

*	Templates and communication files for Monte Carlo codes. Specific files are included for utilization of KSource tools in McStas y TRIPOLI-4 simulations.

*	Command line API: Allows easily producing samples, based on sources saved with Python. Also allows to access templates and communication files, as well as MCPL applications.



## Installation (Linux):
	
1.	Go to source directory:

	```bash
	cd /path/to/ksourcesource
	```

	Where `/path/to/ksourcesource` is the folder where the source distribution of KSource was extracted.

2.	Install with `cmake`:

	```bash
	$ mkdir build && cd build
	$ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/ksourceinstall
	$ make install
	$ cd ..
	```
	Where `path/to/ksourceinstall` is the folder where you wish to install KSource internal files.

	It is required to have previously installed `libxml2`.

3.	Install Python API with `pip`:

	```bash
	$ cd python
	$ pip install .
	$ cd ..
	```

	It is also possible to install Python API directly from PyPI with:

	```bash
	$ pip install ksource
	```

4.	KSource is ready to be used in `/path/to/ksourceinstall`. For example, you can see the `kstool` command options with:

	```bash
	$ /path/to/ksourceinstall/bin/kstool --help
	```

	If you wish to have KSource tools available in your path, execute:

	```bash
	$ export PATH=/path/to/ksourceinstall/bin:$PATH
	```
	Or add this command to `~/.profile` (and update with `source ~/.profile`).
