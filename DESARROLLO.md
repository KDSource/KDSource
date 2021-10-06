# Esquema de la distribución

A continuación se describen los principales archivos y directorios en esta distribución de código fuente:

```
KSource/
+-- LEEME.md                 : Instrucciones de uso.
+-- DESARROLLO.md            : Instrucciones de desarrollo.
+-- CMakeLists.txt           : Archivo de configuración para instalación con
|                              cmake.
+-- KSourceConfig.h.in       : Archivo para generar cabecera de configuraciones
|                              para biblioteca en C.
+-- docs/                    : Documentación.
+-- examples/                : Ejemplo de uso de APIs en Python y línea de
|                              comando con fuentes de partículas generadas con
|                              distribuciones analíticas.
+-- mcpl/                    : Distribución de código fuente de MCPL 1.3.2 con
|                              hooks agregados para TRIPOLI-4, PTRAC de MCNP y
|                              formato SSV.
+-- mcstas/                   
|   +-- contrib/             : Componentes de McStas para comunicación con
|                              formatos MCPL y KSource.
+-- python/                  : API en Python para creación y optimización de
|   |                          fuentes KSource, y estimación de densidad.
|   +-- setup.py             : Archivo de configuración para instalación de API
|   |                          en Python con pip.
|   +-- LEEME.md             : Instrucciones de uso de API en Python.
|   +-- MANIFEST.in          : Archivo de configuración para incluir archivos
|   |                          a la distribución en Python.
|   +-- ksource/             : Código fuente del paquete ksource en Python.
+-- src/                     : Código fuente de APIs en C y línea de comando.
|   +-- ksource/             : Código fuente (.h y .c) de API en C.
|   +-- kstool/              : Código fuente (.c y .sh) de API en línea de
|                              comando.
+-- templates/               : Archivos plantilla para uso de KSource en Python
|                              y ejecución de McStas y TRIPOLI-4.
+-- tests/                   : Unit testing de API en C.
```

# Instalación, testeo y distribución

### Instalación

Instalar APIs en C y línea de comando:
```bash
$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
$ make install
$ cd ..
```
Instalar API en Python:
```bash
$ cd python
$ pip install [-e] .
$ cd ..
```
Para tener la herramienta de línea de comando de KSource en el path ejecutar (o agregar a `~/.profile` y actualizar con `source ~/.profile`):
```bash
$ export PATH=$PATH:/path/to/install/bin
```

### Testeo

Testear APIs en C y línea de comando:
```bash
$ mkdir build && cd build
$ cmake ..
$ make
$ make test
$ cd ..
```
Testear API en Python:
```bash
$ cd python
$ pytest -v
$ cd ..
```

### Distribución

???