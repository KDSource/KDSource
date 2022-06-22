# Usage


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

