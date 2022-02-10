<p align="center">
<img src="../img/kdsource_logo.png" alt="logo" height="100"/>
</p>

[![Python 3.6](https://img.shields.io/badge/python-3.6%2B-blue)](https://www.python.org/downloads/release/python-360/)
[![License](https://img.shields.io/badge/License-GPL3-blue.svg)](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))

This is the Python API of [KDSource package](https://github.com/inti-abbate/KDSource) for Monte Carlo radiation calculations.

KDSource models distribucional particle sources thru Kernel Density Estimation (KDE), based on particle lists with [`MCPL`](https://mctools.github.io/mcpl/) format. It allows to couple different simulations implementing variance reduction, helping Monte Carlo calculations in big systems.

This API allows creating, optimizing, analyzing, plotting, and saving `KDSource` sources. Internally, it uses the [`KDEpy`](https://kdepy.readthedocs.io/en/latest/) library for (adaptive) KDE. The optimization consists in automatic selection of KDE bandwidth. For this purpose, the following methods are included:
*	Silverman's Rule: Simple, fast and not very precise. It is based only on the number of particles and geometry dimension.
*	K Nearest Neighbors: For each sample, computes its bandwidth as the distance to the k-th neighbor.
*	Maximum Likelihood Cross-Validation: Builds a grid of bandwidths, and computes the cross-validation likelihood score for each, which is an indicator of the quality of the estimation. Then chooses the bandwidth that maximizes it. It requires a seed bandwidth, which can be provided by any of the previous methods.

## Usage example

```python
# Create and fit KDSource

plist = kds.PList("mcplfile.mcpl")      # PList: Wrapper for MCPL file
geom = kds.Geometry([kds.Lethargy(),     # Geometry: define metrics for variables
					kds.SurfXY(),
					kds.Polar()])
s = kds.KDSource(plist, geom, bw='silv') # Create KDSource
s.fit(N=1E4)                           # Fit with N particles, optimizing bandwidth

# Create plots

# Energy plot
EE = np.logspace(-3,0,30)
fig,[scores,errs] = s.plot_E(EE)

# Custom variable 1D plot
tt = np.linspace(0,180,30)
fig,[scores,errs] = s.plot_integr('theta', tt)

# Custom variable 2D plot
xx = np.linspace(-10,10,30)
yy = np.linspace(-10,10,30)
fig,[scores,errs] = s.plot2D_integr(['x','y'], [xx,yy])

# Save KDSource to XML file
s.save("xmlfile.xml")

```
