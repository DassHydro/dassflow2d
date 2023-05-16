# DassFlow2D-WRAP





\tableofcontents



# General information

This git contains the development version of DasssFlow2D-Wrapped

The documentation can be found at the following adress: **./doc/DOC_SPHINX/build/index.html**



## Repository organisation

- code  : contains source code and contains the bin directory where the simulation happen
- doc   : contains documentation
- case  : contains reference cases
- Tools : contains scripts for pre and post processing of dassflow 2d wrap

## How to compile and open documentation

Necessary installs:

````
pip install -U Sphinx
pip install numpydoc
pip install pydata_sphinx_theme
pip install sphinx-panels
pip install IPython
pip install sphinxcontrib-bibtex
python3 -m pip install sphinx-autosummary-accessors
````

Open a terminal at the following adress: **./doc/DOC_SPHINX/**

**WARNING** Dassflow2d must have been installed (``make install`` in ./code/ directory )

````
make clean html
````
