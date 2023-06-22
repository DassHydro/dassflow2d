# DassFlow2D-WRAP
This git repository contains the development version of DassFlow2D-Wrapped.

## Table of Contents
1. [ General Information. ](#geninfo)
2. [ Repository organisation. ](#reporg)
3. [ Requirements. ](#req)
4. [ Install project. ](#ins)
5. [ Compile and open documentation. ](#doc)
6. [ First steps and more. ](#next)

<a name="geninfo"></a>
## 1. General information
<strong>WARNING: </strong>For the moment, dassflow2d and other versions work on <strong><ins>Linux only</ins></strong>.

You can clone this project to your own machine using the following command:
````
git clone http://github.com/DassHydro-dev/dassflow2d
````

<strong>Full</strong> and <strong>up-to-date</strong> documentation can be found in [index.html](doc/SPHINX_DOCUMENTATION/build/html/index.html). To open it, enter the following command in the terminal (<em>in your repository directory</em>):
````
open ./doc/SPHINX_DOCUMENTATION/build/html/index.html
````

<a name="reporg"></a>
## 2. Repository organisation
<ul>
  <li>Tools: contains scripts for pre and post processing of dassflow 2d wrap.</li>
  <li>cases: contains reference cases.</li>
  <li>code: contains source code & the bin directory where the simulations happen.</li>
  <li>doc: contains SPHINX documentation (<em>See <a href="#doc">Compile and open documentation</a> for more information</em>).</li>
</ul>

<a name="req"></a>
## 3. Requirements
### 3.1 For dassflow2d installation
<em>**Note:** some of the modules below might already be installed on your Linux machine.</em>
#### 3.1.1 For the Fortran code
- <strong>python 3.8</strong> or above (<em>check your python version with <code>python3 --version</code> in the terminal</em>)
- an up-to-date Java Development Kit (JDK) (<em>check your java version with <code>java --version</code> in the terminal or install it with <code>sudo apt install default-jdk</code></em>)
- an MPI library : mpich
````
sudo apt install -y mpich
````
- Follow the <a href="https://tapenade.gitlabpages.inria.fr/tapenade/distrib/README.html">tutorial</a> to download and install Tapenade and add the following lines to your ````~/.bashrc```` to add tapenade to your PATH:
````
alias tapenade="tapenade_dir/bin/tapenade"
TAPENADE_HOME=tapenade_dir/bin
export PATH=$PATH:$TAPENADE_HOME
export PATH=$PATH:$"tapenade_dir"
````
<em>**Note:** tapenade_dir is the absolute path to the directory containing the tapenade files you just downloaded.</em>

#### 3.1.2 For the Python wrapped code
- pip3
- f90wrap
````
  pip install f90wrap
````
Add f90wrap to your PATH in ````~/.bashrc````. For example:
````
F90WRAP_HOME=~/.local/bin
export PATH=$PATH:$F90WRAP_HOME
````
<em>**Note:** don't forget to enter the command <code>source ~/.bashrc</code> to reload your .bashrc after modifying it.</em>

### 3.2 For SPHINX documentation compilation
In your terminal, execute the following commands:
````
pip install -U Sphinx
pip install numpydoc
pip install pydata_sphinx_theme
pip install sphinx-panels
pip install IPython
pip install sphinxcontrib-bibtex
pip install jupyter_sphinx
python3 -m pip install sphinx-autosummary-accessors
````

<a name="ins"></a>
## 4. Install project
<ul>
  <li> Make sure all the <a href="#req">requirements</a> are met.</li>
  <li> Execute the following commands in the terminal (<em>in your repositoy directory</em>):</li>
</ul>

````
cd ./code
make install
````
<em>**Note:** project installation has default parameters that you can change in [Makefile.inc](code/Makefile.inc) before installation</em>.

<a name="doc"></a>
## 5. Compile and open documentation
<ol>
  <li> Make sure all the <a href="#req">requirements</a> are met.</li>
  <li> Make sure the project has correctly been installed (<em>see <a href="#ins">Install Project</a></em>).</li>
  <li> Execute the following commands in the terminal (<em>in your repository directory</em>):</li>
</ol>

````
cd ./doc/SPHINX_DOCUMENTATION/
make clean html
````
<a name="next"></a>
## 6. First steps and more
<em>Please refer to the <a href="#doc">SPHINX documentation</a> for more information, simple test cases and further details.</em>

<em>\*\*TBA : website, for more on dassflow2d\*\*</em>
