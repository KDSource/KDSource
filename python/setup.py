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
        "numpy>=1.14.2",
        "scipy>=1.0.1",
        "matplotlib>=2.2.0",
        "KDEpy>=1.1.0",
        "joblib>=1.0.1",
        "scikit-learn>=0.24.2",
        "mcpl>=1.3.2",
        "Pillow>=8.2.0"
    ],
    python_requires=">=3.6",
)
