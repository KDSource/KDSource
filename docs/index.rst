KDSource
========
  
Welcome to KDSource's documentation!
------------------------------------

**KDSource** is a tool for Monte Carlo particle sources generation using Kernel Density Estimation.

KDSource assists Monte Carlo beams and shielding calculations, improving tally results in difficult problems. It allows to model big systems (e.g.: investigation reactor guides hall) through spatial or temporal coupling of different simulations in different transport codes, implementing as well variance reduction.

It processes particle lists recorded as output of a simulation (e.g.: passing thru a window), to be used as input in another one. It estimates density distribution in energy, position and direction by means of Kernel Density Estimation (KDE) technique, allowing visualizing it as well as using it to produce new particles (artificial, but with the same estimated density). This allows to increase the number of source particles in the second simulation, improving its statistics (variance reduction).

How it works
------------

The :doc:`intro` section shows how the KDSource tool works in a conceptual way, as well as the theoretical framework.

Installation
------------

Check out the :doc:`installation` section for installation instructions.

Usage
-----

Check out the `Tutorial <tutorial.ipynb>`_ for learning the basic usage of the KDSource tool. In the `User Guide <user-guide.ipynb>`_ and :doc:`api-ref` you can find more detailed information and documentation for its several components. Finally, the :doc:`examples` section shows several working examples for different applications.

About us
--------

The KDSource team was born in the Bariloche Atomic Center (CAB) in Bariloche, Argentina, mainly from works for the Balseiro Institute. See the :doc:`about-us` section for more information about us.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   :maxdepth: 2

   intro
   installation
   tutorial/tutorial.ipynb
   user-guide/user-guide.ipynb
   examples
   api-ref
   about-us
