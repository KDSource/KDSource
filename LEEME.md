# KSource

Esta es la versión de fuente del paquete KSource para cálculos de radiación por Monte Carlo.

KSource brinda herramientas de asistencia al cálculo de blindajes por el método Monte Carlo. Permite modelar grandes sistemas (por ej.: hall de guías de reactor de investigación) mediante el acople, espacial o temporal, de diferentes simulaciones provenientes de diferentes códigos de transporte, implementando además reducción de varianza.

Procesa listas de partículas registradas como output de una simulación (por ej.: al atravesar una ventana), a ser utilizadas como input en otra. Estima la distribución de densidad en energía, posición y dirección por la técnica de Kernel Density Estimation (KDE), permitiendo tanto visualizarla como utilizarla para producir nuevas partículas (artificiales, pero con la misma densidad estimada). Esto permite aumentar la cantidad de partículas de fuente en la segunda simulación, mejorando su estadística (reducción de varianza).

KSource utiliza el formato de listas de partículas [`MCPL`](https://mctools.github.io/mcpl/). Dicho formato, en su versión modificada que se incluye en esta distribución, permite trabajar con los siguientes códigos Monte Carlo:
*	MCNP
*	PHITS
*	McStas
*	TRIPOLI-4

En TRIPOLI y McStas es posible el muestreo "on-the-fly" durante la corrida, mientras que para los demás formatos es necesario grabar la lista de partículas de fuente antes de la simulación.


## Contenidos:

El paquete KSource se compone de las siguientes herramientas:

*	API en Python: Permite crear, optimizar, analizar, graficar, y guardar fuentes KDE. La optimización consiste en la selección automática del ancho de banda. Internamente, usa la biblioteca [`KDEpy`](https://kdepy.readthedocs.io/en/latest) para el KDE.

*	API en C: Permite cargar las fuentes guardadas con Python, y generar nuevas muestras artificiales. Estas respetan la distribución estimada, y pueden guardarse en un nuevo archivo MCPL o bien introducirse directamente en una simulación.

*	Plantillas y archivos de comunicación con códigos Monte Carlo. Se incluyen archivos específicos para la utilización de las herramientas KSource en simulaciones con McStas y TRIPOLI-4.

*	API en línea de comando: Permite producir muestras de forma simple, en base a fuentes guardadas con Python. Permite además acceder a las plantillas y archivos de comunicación, así como a las aplicaciones de MCPL.



## Instalacion (Linux):
	
1.	Ir a carpeta de fuente:

	```bash
	cd /path/to/ksourcesource
	```

	Donde `/path/to/ksourcesource` es la carpeta donde se extrajo la versión de código fuente de KSource.

2.	Instalar con `cmake`:

	```bash
	$ mkdir build && cd build
	$ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/ksourceinstall
	$ make install
	$ cd ..
	```
	Donde `path/to/ksourceinstall` es la carpeta donde se desea instalar los archivos internos de KSource.

	Es requerido haber previamente instalado `libxml2`.

3.	Instalar API en Python con `pip`:

	```bash
	$ cd python
	$ pip install .
	$ cd ..
	```

	También es posible instalar el API en Python directamente desde PyPI con:

	```bash
	$ pip install ksource
	```

4.	KSource está listo para utilizar en `/path/to/ksourceinstall`. Por ejemplo, puedes ver las opciones del comando `kstool` con:

	```bash
	$ /path/to/ksourceinstall/bin/kstool --help
	```

	Si quieres tener disponibles las herramientas de KSource en tu path, ejecuta:

	```bash
	$ export PATH=/path/to/ksourceinstall/bin:$PATH
	```
	O agrega este comando a `~/.profile` (y actualiza con `source ~/.profile`).
