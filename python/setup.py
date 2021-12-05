#!/usr/bin/env python3

"""KSource installer script using setuptools.
"""

import setuptools
import os, re

HERE = os.path.abspath(os.path.dirname(__file__))

# Get long description from README
with open(os.path.join(HERE, "README.md"), "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Get version
with open(os.path.join(HERE, "ksource/__init__.py"), encoding="utf-8") as file:
    VERSION = re.search(r"__version__ = \"(.*?)\"", file.read()).group(1)


setuptools.setup(
    name="ksource",
    version=VERSION,
    author="Osiris Inti Abbate",
    author_email="intiabbate@gmail.com",
    description="Python API for KSource package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/inti-abbate/KSource",
    project_urls={},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License ::",
        "Operating System :: OS Independent",
    ],
    packages=["ksource"],
    install_requires=[
        "numpy>=1.20.3,<2.0",
        "scipy>=1.6.3,<2.0",
        "matplotlib>=3.4.2,<4.0",
        "KDEpy>=1.1.0,<2.0",
        "joblib>=1.0.1,<2.0",
        "scikit-learn>=0.24.2,<1.0",
        "mcpl>=1.3.2,<2.0",
        "Pillow>=8.2.0,<9.0"
    ],
    python_requires=">=3.6",
)
