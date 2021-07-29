# KSource

Esta es la API en Python del [paquete KSource](https://github.com/inti-abbate/KSource) para cálculos de radiación por Monte Carlo.

KSource modela fuentes distribucionales de partículas mediante el método Kernel Density Estimation (KDE), basadas en listas de partículas en formato [`MCPL`](https://mctools.github.io/mcpl/). Permite acoplar distintas simulaciones aportando reducción de varianza, ayudando el cálculo Monte Carlo en sistemas grandes.

Esta API permite crear, optimizar, analizar, graficar, y guardar fuentes `KSource`. Internamente, utiliza la librería [`KDEpy`](https://kdepy.readthedocs.io/en/latest/) para el KDE (adaptativo). La optimización consiste en la selección automática del ancho de banda del método KDE. Para ello se incluyen los siguientes métodos:
*	Regla de Silverman: Método simple, rápido y poco preciso. Se basa únicamente en la cantidad de partículas y la dimensión de la geometría.
*	K Vecinos Más Cercanos: Para cada muestra, computa su ancho de banda como la distancia al k-ésimo vecino.
*	Validación Cruzada de Máxima Probabilidad: Crea una grilla de anchos de banda, y sobre cada uno de ellos computa el "cross validation likelihood score", el cual es un indicador de la calidad de la estimación. Luego elige el ancho de banda que lo maximiza. Requiere un ancho de banda semilla, que puede provenir de cualquiera de lo métodos anteriores.

## Ejemplo de uso

```python
# Crear y ajustar KSource

plist = ks.PList("mcplfile.mcpl")      # PList: Wrapper para archivo MCPL
geom = ks.Geometry([ks.Lethargy(),     # Geometry: define metricas para variables
					ks.SurfXY(),
					ks.Polar()])
s = ks.KSource(plist, geom, bw='silv') # Crear KSource
s.fit(N=1E4)                           # Ajustar con N particulas, optimizando ancho de banda

# Crear graficos

# Grafico en energia
EE = np.logspace(-3,0,30)
fig,[scores,errs] = s.plot_E(EE)

# Grafico 1D de variable arbitraria
tt = np.linspace(0,180,30)
fig,[scores,errs] = s.plot_integr('theta', tt)

# Grafico 2D de variables arbitrarias
xx = np.linspace(-10,10,30)
yy = np.linspace(-10,10,30)
fig,[scores,errs] = s.plot2D_integr(['x','y'], [xx,yy])

# Guardar KSource en archivo XML
s.save("xmlfile.xml")

```
