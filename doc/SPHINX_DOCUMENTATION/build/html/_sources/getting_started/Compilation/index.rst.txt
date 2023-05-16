.. _Compilation:

===================
compilation
===================

In this section, we make our first compilation of the code.

	- Open a terminal in *dassflow2/code/*
	- run the following commands (that will import the example test case):

    .. code-block:: bash

                rm -r ./bin_A/*
                cp -r ../cases/tuto_case/1_lake-at-rest/bin_A/* ./bin_A

Once you provided the needed files (only the file 'm_user_data.f90' is necessary for this step), you can execute the following commands to compile the code and generate **dassflow2d** package:

    .. code-block:: bash

                make install


.. danger::

	when compiling, even if **not used**, the makefile must see an empty file name m_user_data.f90 in the **bin_dir** (identical to the one you can see in the imported bin directory)

.....................
Compilation options
.....................


You can modify your file '/code/makefile.inc' to change compilation options. The various availaible options
are defined directly within the file.


To summurize,with the default configuration, you dont never have to change anything except if you want to
perform parrallel run, in this case, you will set up:

- **USE_MPI=1**
- **NB_PROC=2** if you want to use 2 processors


.. note:

	Parrallel mode is functional in direct and inference run up to 4 cores


--------------------------------------------------------------------------------
Check your installation/compilation of DassFlow
--------------------------------------------------------------------------------

Open a terminal and open a python terminal:

    .. code-block:: bash

		python

In your python terminal, import dassflow2d:

    .. code-block:: python

    		import dassflow2d

    .. code-block:: python

            import dassflow2d
            import inspect
            inspect.getfile(dassflow2d) # this give you the path of the module imported
