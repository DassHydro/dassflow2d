:notoc:

.. module:: dassflow

**************************
DassFlow-2d documentation
**************************

**Date**: |today| **Version**: |version|


:Package: `dassflow2d`  is an open source, python library interfacing the Fortran DassFlow2D model. All source files are written in Fortran 2003 and the software is wrapped into Python using f90wrap :cite:p:`Kermode_2020`.

:license: DassFlow2d is under :ref:`CeCILL license<license>`

:Numercial guide: Numerical guide can be downloaded here:  :download:`./numerical_guide/docs/doc-dassflow2d-num.lyx <./numerical_guide/docs/doc-dassflow2d-num.lyx>`



**DassFlow2D**  is a free computational sotfware framework initially dedicated to river hydraulic simulation and especially designed for variational data assimilation (4D-VAR). The forward as the adjoint discrete models can be runned in parallel calling the MPI library.


.. panels::
    :card: + intro-card text-center
    :column:  col-lg-6 col-md-6 col-sm-12 col-xs-12 d-flex p-3

    ---
    :img-top: _static/index_api.svg


    Introduction
    ^^^^^^^^^^^^^^^^^^^^
    Get a general presentation of Dassflow 2D

    +++

    .. link-button:: Introduction
              :type: ref
              :text: To the introduction of Dassflow 2d model
              :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/index_getting_started.svg

    Getting started
    ^^^^^^^^^^^^^^^^^^^^

    New to `dassflow2d`? Check out the getting started guides.

    +++

    .. link-button:: getting_started
              :type: ref
              :text: To the getting started guides
              :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/index_user_guide.svg

    
    User Guide
    ^^^^^^^^^^^^^^^^^^^^

    The user guide provides in depth-information of the *dassflow2d* library


    +++

    .. link-button:: user_guide
              :type: ref
              :text: To the user guide
              :classes: btn-block btn-secondary stretched-link

    ---
    :img-top: _static/index_api.svg

    API reference
    ^^^^^^^^^^^^^^^^^^^^

    The reference guide contains a detailed description of
    the *dassflow2d* API.

    +++

    .. link-button:: api
              :type: ref
              :text: To the reference guide
              :classes: btn-block btn-secondary stretched-link



      ---

      :img-top:  _static/index_contribute.svg

      Additional documentation
      ^^^^^^^^^^^^^^^^^^^^^^^^^

      Numerical guide and more


      +++

      .. link-button:: numerical_guide
                      :type: ref
                      :text: To the Numerical guide
                      :classes: btn-block btn-secondary stretched-link


.. toctree::
   :maxdepth: 4
   :hidden:
   :titlesonly:

   Introduction/index
   getting_started/index
   user_guide/index
   reference/index
   numerical_guide/index
   license/index


.. rubric:: References
.. bibliography::
