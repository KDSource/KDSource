# Configuration file for the Sphinx documentation builder.
import os
import sys

sys.path.insert(0, os.path.abspath("../../python"))
sys.path.insert(0, os.path.abspath("../tutorial"))

# -- Project information

project = "KDSource"
copyright = "2021, Abbate"
author = "Abbate"

release = "1.0"
version = "0.1.0"

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "myst_parser",
    "nbsphinx",
    "numpydoc"
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

autodoc_member_order = 'bysource'

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
epub_show_urls = "footnote"

html_logo = "img/kdsource_logo.svg"
html_theme_options = {
    'logo_only': True,
    'display_version': False,
}

# -- Options for numpydoc
numpydoc_show_class_members = False
