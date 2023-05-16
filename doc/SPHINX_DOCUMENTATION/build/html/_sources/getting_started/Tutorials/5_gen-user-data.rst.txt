.. _5_gen-user-data:

==============================================================================
Generate  Basic channel
==============================================================================

you can generate basic meshing and
necessary input files (defining other parameters files for friction, boundary conditions, observed data etc...) usgin the source code at `/dassflow2d-wrap/Tools/1_pre-treatment/0_gen_basic_channel` .

++++++++++++++++++++++++++++++
Generate input files
++++++++++++++++++++++++++++++

Go to the directory `/dassflow2d-wrap/Tools/1_pre-treatment/0_gen_basic_channel` and open a terminal.


run the following command in your terminal (still in directory `/dassflow2d-wrap/Tools/2_gen_basic_channel`).

.. code-block:: bash

  python3 -m numpy.f2py -c -m gen_channel_case _gen_channel_case.f90


This generates a wrapped package called *gen_channel_case* of subroutines written in *_gen_channel_case.f90*. Then, perform the command:

.. code-block:: bash

  python3 run.py



To generate your own parameterized mesh, most parameter you can define are define in :ref:`User Guide<user_guide>`.


Then, copy the files generated in  ``/dassflow2d-wrap/Tools/2_gen_basic_channel/files`` to your bin directory: ``/dassflow2d-wrap/code/bin_A``.

.. note::

  you might have to modify the last line of input.txt, either to add or remove a coma


then to perform simulation, open a terminal in  ``/dassflow2d-wrap/code/``, and execute the command:


.. code-block:: bash

  make rundirect
