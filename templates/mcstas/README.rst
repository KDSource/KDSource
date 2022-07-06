Using this template
-------------------

In order to use this template, both KDSource and McStas must be installed locally.

About the instrument file
*************************

In order to keep things simple, this instrument file only uses a previously built KDSource for a McStas simulation. The instrument file consist of the source component, a position sensitive detector component and an energy detector. 

The KDSource.comp component that will be used needs to link to the correct libraries in order to work. Therefore, a DEPENDENCY line must be added to the instrument file with the following flags:

.. code-block:: console

    DEPENDENCY "-lkdsource -lmcpl -lm -I/path/to/kdsourceinstall/include -Lpath/to/kadsourceinstall/lib"


where `path/to/kdsourceinstall` must be modified to your current KDSource installation folder. An example is given in which KDSource was installed in the `/usr/local/KDsource` folder.

Understanding the exe_McStas.sh
*******************************

We give an executable working on linux (modifications may be needed for Mac and are definetly needed for windows) in which the example instrument `kds_instrument_example.instr` is executed containing a KDSource component. This source uses a KDSource file generated with the kdsource python API named `source.xml`. This source file was based on an mcpl particle list given in `samples.mcpl.gz`. Both these files are needed for the simulation.  

The first lines of the `exe_McStas.sh` file are:

.. code-block:: bash

    ##### Input #####
    INPUT='kds_instrument_example'    # McStas istrument name, without .instr
    OutDir='example' # Directory where execute simulation and store output files
    N=1E6          # Number of particles to 
    KDSourceLibPath=$1
    ### End Input ###


In these lines, the simulation output is saved in the `example` folder and :math:`1\times 10 ^6` particles are simulated. The ``KDsourceLibPath`` bash variable must be given as command line input, or may be modified to KDsource's lib folder in the KDSource installation path. 

The next lines of the bash script are:

.. code-block:: bash

    ## Execute McStas
    rm -rf $OutDir
    export LD_LIBRARY_PATH=$KDSourceLibPath
    cp ../../mcstas/contrib/KDSource.comp .
    mcrun -c $INPUT --dir=$OutDir -n $N
    rm KDSource.comp kds_instrument_example.c kds_instrument_example.out
    echo Simulation ended successfully


in which the example folder is removed if already existing. Then the library path is set for the KDSource.comp component to find libkdsource.so, and then the last version of the component is brought to the current folder. The McStas simulation is compiled and run in the `mcrun` command and finally the all generated files are removed.

For example, if the kdsource installation path is ``/usr/local/KDSource``, then the bash script may be executed in the current folder by command line:

.. code-block:: console

    $./exe_McStas.sh '/usr/local/KDSource/lib'


The simulation output may be explored in the `example` folder with ``mcplot`` or by manually inspecting the files with any software.