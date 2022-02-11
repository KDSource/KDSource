Introduction
============
  
Welcome to KDSource's documentation!
------------------------------------

**KDSource** is a tool for Monte Carlo particle sources generation using Kernel Density Estimation.

KDSource assists Monte Carlo beams and shielding calculations, improving tally results in difficult problems. It allows to model big systems (e.g.: investigation reactor guides hall) through spatial or temporal coupling of different simulations in different transport codes, implementing as well variance reduction.

It processes particle lists recorded as output of a simulation (e.g.: passing thru a window), to be used as input in another one. It estimates density distribution in energy, position and direction by means of Kernel Density Estimation (KDE) technique, allowing visualizing it as well as using it to produce new particles (artificial, but with the same estimated density). This allows to increase the number of source particles in the second simulation, improving its statistics (variance reduction).

.. note::

   Check out the :doc:`usage` section for further information, as well as the :doc:`Tutorial.ipynb` and the :doc:`installation` instructions.

   This project is under active development.

Contents
--------

.. toctree::

   self
   installation
   Tutorial.ipynb
   usage
   documentation
