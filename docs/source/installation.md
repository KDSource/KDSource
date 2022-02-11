# Installation

## Linux

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
   Where `/path/to/kdsourcesource` is the folder where the source distribution of KDSource was cloned, and `path/to/kdsourceinstall` is the folder where you wish to install KDSource internal files.

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

## Windows

Under development