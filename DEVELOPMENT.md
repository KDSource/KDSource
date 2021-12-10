# Distribution scheme

The main files and directories in this source code distribution are described  as follows:

```
KDSource/
+-- README.md                : Usage instructions.
+-- DEVELOPMENT.md           : Development instructions.
+-- CMakeLists.txt           : Configuration file for installation with cmake.
+-- KDSourceConfig.h.in      : File to generate configurations header for C
|                              library.
+-- docs/                    : Documentation.
|   +-- examples/            : Example of Python and command line APIs usage
|                              with particle sources generated with analytic
|                              distributions.
+-- mcpl/                    : Source code MCPL distribution 1.3.2 with added
|                              hooks for TRIPOLI-4, MCNP's PTRAC and SSV format.
+-- mcstas/                   
|   +-- contrib/             : McStas components for communication with MCPL and
|                              KDSource formats.
+-- python/                  : Python API for creation and optimization of
|   |                          KDSource particle sources, and density estimation.
|   +-- setup.py             : Configuration file for installation of Python API
|   |                          with pip.
|   +-- README.md            : Usage instructions for Python API.
|   +-- MANIFEST.in          : Configuration file for including files to Python
|   |                          distribution.
|   +-- kdsource/            : Source code of kdsource package in Python.
+-- src/                     : Source code of C and command line APIs.
|   +-- kdsource/            : Source code (.h y .c) of C API.
|   +-- kdtool/              : Source code (.c y .sh) of command line API.
+-- templates/               : Template files for KDSource usage in Python and
|                              McStas and TRIPOLI-4 execution.
+-- tests/                   : C API unit testing.
```

# Installation, testing and distribution

### Installation

Install C and command line APIs:
```bash
$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
$ make install
$ cd ..
```
Install Python API:
```bash
$ cd python
$ pip install [-e] .
$ cd ..
```
To have KDSource command line tool in your path execute (or add to `~/.profile` and update with `source ~/.profile`):
```bash
$ export PATH=$PATH:/path/to/install/bin
```

### Testing

Test C and command line APIs:
```bash
$ mkdir build && cd build
$ cmake ..
$ make
$ make test
$ cd ..
```
Test Python API:
```bash
$ cd python
$ pytest -v
$ cd ..
```

### Distribution

At the time, the only implemented distribution is cloning the GitHub repo.