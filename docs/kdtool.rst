kdtool Command Line Tool
========================

Command line tool for particle sampling using a previously generated KDSource XML file, and access to other helpful utilities.
C library for particle generation with KDSource objects.

The general instructions for the ``kdtool`` command are the following:

::

    Usage: kdtool [options]

    KDSource is a Monte Carlo calculations assistance tool. It implements particles
    density estimation and sampling by means of Kernel Density Estimation method.

    Options:
        resample:   Resample particles based on a kdsource XML file.
        templates:  Copy templates for Monte Carlo calculations.
        beamtest:   Test source with simple beam calculation.
        [Any MCPL command]
        -h, --help: Display usage instructions.

Resample
--------

The ``resample`` argument brings a high-level interface for particle sampling, using a prevously created KDSource XML file defining an optimized KDE source. This file is usually generated through the Python API. Internally this command uses the C API.

The instructions for the ``kdtool resample`` command are the following:

::

    Usage: kdtool resample sourcefile [options]

    Resample particles from source defined in XML file sourcefile, and save them in
    a MCPL file.

    Options:
        -o outfile: Name of MCPL file with new samples
                    (default: \"resampled.mcpl\").
        -n N:       Number of new samples (default: 1E5).
        -h, --help: Display usage instructions.

Templates
---------

The ``templates`` argument copies some helpful templates for several usages of the KDSource tool. This includes Jupyter Notebooks files for source fitting and optimization and postprocessing of results, as well as files required for communication and coupling with the Monte Carlo codes McStas and TRIPOLI-4. These files can also be found in the ``templates`` subdirectory of the installed KDSource location.

The instructions for the ``kdtool templates`` command are the following:

::

    Usage: kdtool templates dest [options]

    Copy to dest templates for KDSource usage in Python, or for interacting with
    Monte Carlo codes.

    Options:
        --mcstas:   Copy templates for using McStas.
        --tripoli:  Copy templates for using TRIPOLI-4.
        --all:      Copy all templates.
        -h, --help: Display usage instructions.

Beamtest
--------

The ``beamtest`` argument gives access to a module for comparison between the original tracks source and a KDE source. It performs a small simulation in which particles are sampled from the sources and transported in vacuum, and a tally is computed in a rectangular collimator placed at some distance from the source. The simulation is performed twice, once for each source, and the results for both sources are saved in a text file. The expected result, for a correctly optimized KDE source, is that all tally results errorbands overlap between both sources, and that the KDE source presents a smaller statistical error.

The instructions for the ``kdtool beamtest`` command are the following:

::

    Usage: kdtool beamtest sourcefile [options]

    Executes a simple simulation with source defined in XML file sourcefile, in
    which calculates the number of particles passing thru a rectangular collimator.
    The simulation is repeated using the particle list directly as source, to
    compare the results.

    Results are computed for 4 energy groups, and stored in a results file which
    can be imported from a spreadsheet.

    This tool is designed to be used with flat neutron sources with particles
    propagating towards z direction.

    Options:
        -n N:           Number of source particles (default: 1E6).
        -o results:     Name of file to store results.
        -xwidth value:  Width of the collimator, in cm (default: 7).
        -yheight value: Height of the collimator, in cm (default: 20).
        -z value:       Position of the collimator along z axis, in cm
                        (default: 500).
        -xshift:        Horizontal shift of the center of the collimator,
                        in cm (default: 0)
        -yshift:        Vertical shift of the center of the collimator,
                        in cm (default: 0)
        -h, --help:     Display usage instructions.

MCPL access
-----------

The ``kdtool`` command gives access to any MCPL command from the modified MCPL distribution that is included in KDSource. For example use ``kdtool mcpltool --help`` to see the usage instructions of ``mcpltool``, or use ``kdtool ssv2mcpl input.ssv [output.mcpl]`` to convert a SSV text file to the MCPL format.