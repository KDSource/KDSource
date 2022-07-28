# Distribution scheme

The main files and directories in this source code distribution are described  as follows:

```
KDSource/
+-- LICENSE                     : KDSource Licence (GNU GPL-3).
+-- README.md                   : Usage instructions.
+-- DEVELOPMENT.md              : Development instructions.
+-- CITATION.cff                : Citation file for the KDSource project.
+-- CMakeLists.txt              : Configuration file for installation with cmake.
+-- KDSourceConfig.h.in         : File to generate configurations header for C.
+-- pyproject.toml              : Configuration file for Python.
|                                 library.
+-- INSTALL.ps1                 : PowerShell installation script for Windows.
+-- docs/                       : Documentation.
|   +-- Makefile                : Configuration file for Sphynx documentation.
|   +-- make.bat                : Configuration file for Sphynx documentation.
|   +-- requirements.txt        : Requirements file for Sphynx documentation.
|   +-- source/                 : Source code for Sphynx documentation.
|   +-- examples/               : Usage examples.
|   +-- KDSource_manual_spanish : KDSource manual (spanish) LaTeX source code.
+-- img/                        : KDSource logo files.
+-- mcpl/                       : Source code MCPL distribution 1.3.2 with added
|                                 hooks for TRIPOLI-4, MCNP's PTRAC and SSV format.
+-- mcstas/
|   +-- contrib/                : McStas components for communication with MCPL and
|                                 KDSource formats.
+-- python/                     : Python API for creation and optimization of
|   |                             KDSource particle sources, and density estimation.
|   +-- README.md               : Usage instructions for Python API.
|   +-- MANIFEST.in             : Configuration file for including files to Python
|   |                             distribution.
|   +-- setup.py                : Configuration file for installation of Python API
|   |                             with pip.
|   +-- pyproject.toml          : Configuration file for Python API installation.
|   +-- tox.ini                 : Configuration file for GitHub Workflow.
|   +-- kdsource/               : Source code of kdsource package in Python.
|   +-- tests/                  : Python API unit testing.
+-- src/                        : Source code of C and command line APIs.
|   +-- kdsource/               : Source code (.h y .c) of C API.
|   +-- kdtool/                 : Source code (.c y .sh) of command line API.
+-- templates/                  : Template files for KDSource usage in Python and
|                                 McStas and TRIPOLI-4 execution.
+-- tests/                      : C API unit testing.
+-- .github/workflows/ci.yml    : Workflow file for GitHub Actions.
```

# Workflow (Linux)

This is the recommended workflow for developing and contributing to KDSource.

## Installation

Clone GitHub repository:
```bash
git clone --recurse-submodules https://github.com/KDSource/KDSource
git remote add origin https://github.com/KDSource/KDSource
```
Install C and command line APIs:
```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/kdsourceinstall
make install
cd ..
```
Install Python API:
```bash
cd python
pip install .
cd ..
```
To have KDSource command line tool in your path execute (or add to `~/.profile` and update with `source ~/.profile`):
```bash
export PATH=$PATH:/path/to/install/bin
```

## Development

Create and switch to feature branch:
```bash
git pull origin master
git branch newfeature
git switch newfeature
```
Here you can modify, delete or create files in the working tree.

## Testing

Test C and command line APIs:
```bash
mkdir build && cd build
cmake ..
make
make test
cd ..
```
Check that make process completes correctly and that all tests are passed.

Test Python API:
```bash
tox -r -c python
```
Check that style, tests and coverage test succeed. If style fails, you can fix the style error manually or use `black` tool to format Python code.

Check documentation:
```bash
sphinx-build -b html docs docs/build
```
The HTML files will be built in `docs/build`, being `index.html` the main page. These are the files that will be built in the [Documentation Page](https://kdsource.readthedocs.io/en/latest/). Check that the presented documentation is correct.

## Contribution

Save your work in feature branch:
```bash
git status
git add -u
git add [new-files]
git git commit
git push origin newfeature
```
Merge chanches to master branch:
```bash
git rebase [-i] master
git switch master
git merge newfeature
git push origin master
```
Delete feature branch:
```bash
git branch -d newfeature
git push origin -d newfeature
```


## Distribution

At the time, the only implemented distribution the GitHub repo.
