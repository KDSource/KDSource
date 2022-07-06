# Installation

Currently, the only implemented installation method is via cloning the GitHub repository and building with CMake and Pip. See bellow for specifical instructions for Linux and Windows.

## Linux and Mac

Requirements: Git 2.14+, GCC 9+, CMake 3+, Pip 22+ (Python 3.8+), LibXml2 2.9.3.

If you are on linux OS you can install `libxml2` with:

```bash
   $ sudo apt-get update
   $ sudo apt-get install libxml2
```

for Ubuntu, or similarly for other Linux distributions, using the corresponding package manager.

In case you are in MacOS `libxml2` can be installed by using [Homebrew](https://brew.sh/):

```bash
   $ brew install libxml2
```
The following mandatory steps are valid for both MacOS and linux systems:

1. First of all, clone this repository with all its submodules to a local repository.

   ```bash
   $ git clone --recurse-submodules https://github.com/KDSource/KDSource
   ```

2. Go to source directory and install with `cmake`:

   ```bash
   $ cd /path/to/kdsourcesource
   $ mkdir build && cd build
   $ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/kdsourceinstall
   $ make install
   $ cd ..
   ```
   Where `/path/to/kdsourcesource` is the folder where the source distribution of KDSource was cloned, and `/path/to/kdsourceinstall` is the folder where you wish to install KDSource internal files.

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


## Windows

Requirements: Git 2.14+, MinGW-GCC 11+, CMake 3+, Pip 22+ (Python 3.8+), LibXml2 2.9.3.

You can install KDSource, including LibXml2, using the `INSTALL.ps1` PowerShell script:
   ```bash
   > .\INSTALL.ps1
   ```
The script requires 7-Zip, with `7z` command in the system path.

If you prefer a step-by-step approach, you can download `libxml2` in the following link:
* 64 bits: http://xmlsoft.org/sources/win32/64bit/
* 32 bits: http://xmlsoft.org/sources/win32/
Download and extract the `libxml2` and `iconv` archives, and add the path to the `bin` subdirectory of each library to the system `PATH` variable.

Important: The architecture (32 vs 64 bits) of the installed `libxml2` and `iconv` must be the same as the MinGW and CMake architecture. Also make sure that other `libxml2` or `iconv` files with different architecture are not in the `PATH`, or at least not ahead of the ones to be used.

The following instructions use the PowerShell command prompt, and therefore assume that the `bin` subdirectory of Git, MinGW and CMake are in the system `PATH`.

1. First of all, clone this repository with all its submodules to a local repository.

   ```bash
   > git clone --recurse-submodules https://github.com/KDSource/KDSource
   ```

2. Go to source directory and install with `cmake`:

   ```bash
   > cd C:\\path\\to\\kdsourcesource
   > mkdir build && cd build
   > cmake .. -DCMAKE_INSTALL_PREFIX=C:\\path\\to\\kdsourceinstall -G "MinGW Makefiles"
   > set C_INCLUDE_PATH=C:\\path\\to\\iconv\\include
   > mingw32-make install
   > cd ..
   ```
   Where `C:\\path\\to\\kdsourcesource` is the folder where the source distribution of KDSource was cloned, and `C:\\path\\to\\kdsourceinstall` is the folder where you wish to install KDSource internal files. `C:\\path\\to\\iconv` is the folder where `iconv` was extracted.

3. Add the `C:\\path\\to\\kdsourceinstall\lib` subdirectory to the system `PATH`.

4. Install Python API with `pip`:

   ```bash
   > cd python
   > pip install .
   > cd ..
   ```

5. KDSource is ready to be used in `C:\\path\\to\\kdsourceinstall`. For example, you can see the `kdtool-resample` command options with:

   ```bash
   > C:\\path\\to\\kdsourceinstall\\bin\\kdtool-resample --help
   ```

   If you wish to have KDSource tools available in your path, add the `bin` subdirectory to the system `PATH`.

   Note: Currently, the `kdtool` application is not available on Windows, but the resampling tool can be invoked directly with `kdtool-resample`, as shown in the example above.

