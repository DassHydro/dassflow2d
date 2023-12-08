# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'DassFlow-2D'
copyright = 'Under Cecil licence'
author = 'INSA Toulouse, INRAE-Aix-en-Provence'

# The full version, including alpha/beta/rc tags
release = '1.0'
version = release

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.todo',
              'sphinx.ext.autosectionlabel', 
              'sphinx.ext.viewcode', 
              'sphinx.ext.autodoc', 
              'sphinx.ext.duration', 
              'sphinx.ext.autosummary', 
              'numpydoc', 
              'sphinx_panels', 
              'IPython.sphinxext.ipython_directive', 
              'IPython.sphinxext.ipython_console_highlighting', 
              'sphinx_autosummary_accessors', 
              'sphinxcontrib.bibtex', 
              'sphinx.ext.intersphinx',
              "sphinx.ext.napoleon",
              "pyvista.ext.plot_directive",
              'jupyter_sphinx'
              ]


bibtex_bibfiles = ['references.bib']


pygments_style = 'sphinx'

numpydoc_show_class_members=True
autosummary_generate = True  # Turn on sphinx.ext.autosummary

default_role = "autolink"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"
html_theme_options = {
  "collapse_navigation": False,
  "use_edit_page_button": True,
}
html_context = {
    "display_gitlab":True,
    "github_url": "https://github.com/DassHydro/dassflow2d", # or your self-hosted GitLab
    "github_host":"https://github.com/DassHydro/dassflow2d",
    "github_user": "pag13",
    "github_repo": "dassflow2d",
    "github_version": "documentation",
    "doc_path": "/doc/SPHINX_DOCUMENTATION/source/",
}

html_use_modindex = True

panels_add_bootstrap_css = False


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [
    "css/dassflow.css",
]

