.. _3_make_your_first_parallel:

======================================================================
Upgrade to parallel
======================================================================

In this tutorial, we show how to use **parallel mode** :

------------------------------
Set up your environment
------------------------------

You need to have the test case  `/dassflow2d-wrap/cases/tuto_cases/3_inference-python` imported


-----------------------------------
Set up the makefile and compile
-----------------------------------
We aim to set up the correct compilation options, and make our first compilation of the code. To do so, follow the following steps:

    - Open your makefile.inc, which is located at **/dassflow2d-wrap/code/Makefile.inc**
    - set the parameters  ``USE_MPI=1`` and ``NB_PROC=2``


To compile the model,   open a terminal in the following directory: **/dassflow2d-wrap/code/**, and run the following command:

.. code-block:: bash

    make install


-----------------------------------
Launch your first run using make
-----------------------------------

You can either perform direct simulation or minimization:

.. code-block:: bash

  make rundirect
  make runmin



----------------------------------------------
In python script and using python optimizer
----------------------------------------------



To perform the script ``dassflow2d-wrap/cases/tuto_case/3_inference-python/1_tuto_inference-python.py``, in parallel mode,
you must execute the following command, in a terminal opened at ``dassflow2d-wrap/cases/tuto_case/3_inference-python`` :



.. code-block:: bash

  mpirun -np 2 python3 1_tuto_inference-python.py
