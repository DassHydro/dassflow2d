.. _Installation:

==============================
Installation
==============================

`dassflow2d` can be used under Linux and should work with most distributions. The installation instructions are detailed for Ubuntu.
Developpement was carried out on Ubuntu 20.04.4 LTS, with Python 3.8 and 3.9.

Depending on your distribution, you will need to use the correct package manager and insert the appropriate packages. Note that Ubuntu 22.04 LTS relies on Python 3.10, with which the code as not yet been tested.

-----------------------------------
Download `dassflow2d` source code
-----------------------------------

This can be done either via the git repository or via SVN.

+++++++++++++++++++++++++++++++++++++++
From git repository
+++++++++++++++++++++++++++++++++++++++

To download `dassflow2d`:

1. Open a terminal in the directory where you want to install DassFlow
2. Write the command:

.. code-block:: bash

	git clone git@gitlab-ssh.irstea.fr:lilian.villenave/dassflow2d-wrap.git

+++++++++++++++++++++++++++++++++++++++
From SVN repository
+++++++++++++++++++++++++++++++++++++++

.. note::

		SVN repository containing DassFlow is hosted on Sourcesup, which is	a free French academic service.
		You need a Sourcesup account to access the code this way.


-----------------------------------
Requirements for the Fortran code
-----------------------------------
To take advantage of full capabilities of DassFlow2D, some librairies must be installed.

A fortran compiler have to be installed: either the `GNU FORTRAN compiler <http://gcc.gnu.org/fortran/>`_  or the   `INTEL FORTRAN compiler <http://software.intel.com/en-us/intel-compilers>`_.

The automatic code differentiation tool  `TAPENADE <http://www-sop.inria.fr/tropics/>`_.

`TAPENADE` tool depends on JAVA, thus you must have an appropriate JAVA version, and specify some information about its localisation to TAPENADE

.. dropdown:: Detailed install - Fortran requirements
   :animate: fade-in-slide-down
   :container: + shadow
   :title: font-weight-bolder

	**Fortran Compiler**
	You sould have a Fortran compiler installed by default on your machine, if not you must install it from official website: `<http://gcc.gnu.org/fortran/>`_

	**JAVA**

	1. First, update the apt package index with:

	.. code-block:: bash

	 	sudo apt update

	2. Once the package index is updated install the default Java OpenJDK package with:

 	.. code-block:: bash

    		sudo apt install default-jdk

	3. Verify the installation, by running the following command which will print the Java version

 	.. code-block:: bash

	   	java -version

   	4. ADD JAVA_HOME to your environement

   	- Open your .bashrc file:

 	.. code-block:: bash

    	  	gedit ~/.bashrc

	- Add within your file the following lines :

 	.. code-block:: bash

    		JAVA_HOME="/your path to the directory of java executable"
    		# for example:
    	  JAVA_HOME="/usr/bin"
    	  export PATH=$PATH:$JAVA_HOME

	- Source your modification: in command line write:

 	.. code-block:: bash

    		source ~/.bashrc


	- Check your definition is correct

 	.. code-block:: bash

	 echo $JAVA_HOME

  	**TAPENADE**

	- First download and install TAPENADE following the tutorial on [TAPENADE WEBSITE](http://www-sop.inria.fr/tropics/)

	- Open your .bashrc file:

	.. code-block:: bash

    		gedit ~/.bashrc

	- Add within your file the following lines :

	.. code-block:: bash

		 alias tapenade= "/tapenade_dir/bin/tapenade"
    		 TAPENADE_HOME=tapenade_dir/bin`
    		 export PATH=$PATH:$TAPENADE_HOME
    		 export PATH=$PATH:$"tapenade_dir"

	.. tip::

			``tapenade_dir``: is the path to the tapenade directory you downloaded, for exemple:
			`/home/livillenave/Documents/software/tapenade_3.16`

    **MPI librairy**

  	you must install mpich library

    	.. code-block:: bash

    		sudo apt-get install -y mpich



-----------------------------------
Requirements for the wrapped code
-----------------------------------

To take advantage of wrapped version of `dassflow2d`,  a Python interpreter is needed, as well as some specific Python librairies that are used in the code:


.. dropdown:: Detailed install - Wrapped requirements
   :animate: fade-in-slide-down
   :container: + shadow
   :title: font-weight-bolder


    **Python and pip3**

    Python version must be 3.8 or higher

    .. code-block:: bash

            sudo apt install python3.8 python3-pip

    **f90wrap**

    .. code-block:: bash

			pip install f90wrap


Multiple Python librairies are required. They will be installed automatically if needed during the compilation process.
