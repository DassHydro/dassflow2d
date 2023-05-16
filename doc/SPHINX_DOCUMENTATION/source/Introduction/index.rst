.. _Introduction:


============================
Presentation
============================


DassFlow (Data Assimilation for Free Surface Flows) is an open computational software designed to model shallow free surface flows and perform Variational Data Assimilation. The hydrodynamic flow model can be coupled with a hydrological model, therefore providing a complete hydraulic-hydrological chain applicable at basin scale.

The computational kernel is coded in Fortran 2003, with MPI, and wrapped in Python. The following physical models are available: the 2D Shallow Water (SW) equations in variables :math:`(h,q_{x},q_{y})`.
**where** :math:`x, y  \text{ denotes x and y direction in global coordinate system on mesh} \forall(x,y)\in\mathcal{D}_{\Omega}`,  **at time** :math:`t`, physical time.

Note that a 1D version of SW system in variables :math:`(S,Q)` (Saint-Venant's equations) is available but in another DassFlow version (please consult DassFlow webpage for details: https://www.math.univ-toulouse.fr/DassFlow/index.html).

In the latest version (DassFlowV3), the DassFlow platform includes a Herschel-Bulkley rheology version of the hydraulic equations (for non-Newtonian fluids) and the GR4H rainfall-runoff model implemented in a semi-distributed manner.

.. hint::

    **Question**

    - Is Herschel-Bulkley rheology tested in this version ?
.. hint::

    **Note**

    - GR4 is being recoded/improved : not commited yet to the wrapped version


The numerical guide presenting the various numerical schemes and boundary condition treatments available in DassFlow can be downloaded here: :download:`../numerical_guide/docs/doc-Dassflow2d-num.lyx <../numerical_guide/docs/doc-dassflow2d-num.lyx>`

.. DassFlow (Data Assimilation for Free Surface Flows) denotes a set of a few computational codes aiming at modeling free surface geophysical flows (water, ice, lavas etc) with data assimilation capabilities. The flow models are mainly shallow ones (long-wave assumption), non turbulent.
.. The presented code is in its wrapped version. The kernel code is written in fortran and and the fortran code is wrapped and made accessible, executable, and can communicate with a Python interface. This enables   easy use of any other Python libraries.

The present wrapped version is dedicated to river flows simulation and especially designed for variational data assimilation (4D-VAR)

============================
Why DassFlow2D ?
============================

----------------------------------------------------------
A model dedicated to river flows simulation
----------------------------------------------------------

DassFlow2D solves the Shallow Water equations using various finite-volume method on a computational mesh either of triangular or rectangular elements.

The numerical schemes are a mix of classical ones and original ones. They are Finite Volume schemes (1st order/2nd order) for shallow flow models and Finite Element schemes (2nd order).

Boundary conditions required for real-world flows are available. See detailed numerical models statements in :download:`../numerical_guide/docs/doc-dassflow2d-num.lyx <../numerical_guide/docs/doc-dassflow2d-num.lyx>`.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DassFlow2D can take into account the following phenomena:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Bed friction (Manning-Strickler law and others implemented)
- Boundary conditions required for real-world flows are available
- Mixed flow regimes (fluvial and torrential) and transitions
- Wet/dry front propagations
- Spatially distributed rainfall and infiltration treatment (Green-Ampt, and SCS-curve number are implemented)
- 1D-like/2D river networks with a single hydraulic solver
- Coupling with hydrological modeling (spatially distributed or not) directly implemented in the forward/inverse algorithms (INRAE models available GR4H and SMASH)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Additionnal capabilities of DassFlow2D:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

- Data assimilation (4D-var)
- Cartesian coordinate system for meshing (Mesh still can be irregular with triangles or quadrangles)
- Python tools for data manipulation and visualisation ( see  :ref:`Main packages used` section) .

.. hint::

    **to be implemented**

    - Direct inferfacing with Telemac's .vtk output format
    - Telemac interfacing (.slf and .cli files correspondance with Dassflow mesh and bc file format.)
    - SMASH interfacing

------------------------------------------------------------------------------
A model dedicated to inverse problems: Variational Data Assimilation (VDA)
------------------------------------------------------------------------------
An important feature of DassFlow is its capability to address inverse problems through VDA, an adjoint-based method).


.. note::


  Since the Python wrapped versions, combining VDA method (physics-based method) with Deep Learning (data driven estimations) might be much easier.

  Hybrid methods to be developped.


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DassFlow2D can ingest the following data (uneven samplings):
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Input data required for a forward run:

- DEM

.. up to very high resolution (**todo** Lidar, link IGN, also lower dem MERIT)  :math:`b(x,y)\forall(x,y)\in\mathcal{D}_{\Omega}`

- Boundary conditions (discharge, water levels, rating curve)

**Additional input data** for a forward run

- Soil occupation for defining patches (friction and infiltration parameters)
- Precipitation, evapotranspiration and hydrological parameters (hydrological module)

Data types that can be assimilated are any number of combination of the following:

- Water levels
- Discharges (hydraulic and hydrological)
- 2D flow velocities

.. hint::

    **to be implemented**

    - WIP Depth-dependent porosity
    - WIP Hydraulic structures (i.e. weirs, bridges, ...)
    - WIP Water extents assimilation
    - WIP Optical video fluxes assimilation
    - WIP Soil occupation maps reading


.. hint::

  Remark on multisensor + obs operator and error matrices ?

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DassFlow2D enables to simultaneously infer large control vectors:
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Parameters types that can be inferred from observations are any number of combination of the following:

- Distributed two-parameter friction power law
- Distributed bathymetry
- Semi-distributed hydrological parameters
- Boundary conditions

.. hint ::

  The inference algorithm can be performed in large dimension and regularisations can be used.
  --> implémentation matrices covariance lilian

++++++++++++++++++++++++++++++++++++++++++++++++
Inverse algorithm principle
++++++++++++++++++++++++++++++++++++++++++++++++

VDA aims at minimizing the discrepancy between the model output and some flow observations, measured via a cost function :math:`J`. This is achieved via the optimization of the flow model parameters with a quasi-Newton algorithm. It uses :math:`\nabla J`, the gradient of the cost with regard to parameters, computed by soving the adjoint model (obtained by algorithmic diffrenciation or not). The complete optimization process provides parameters identification, reduced uncertainties, calibrated model(s). VDA formulations are based on priors such as first guess values and covariances operators. Classical and original covariances operators are available. Different regularization terms can be easily introduced in the cost function.

+++++++++++++++++++++++++++++++++++++++++++++++++
Gradient computation and optimization algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++

For VDA, the gradients are computed using the adjoint method. The adjoint codes are generated either using the automatic differentiation tool Tapenade (INRIA) or by coding the adjoint equations.
Optimization routines (1st order methods) are either the ones from NumPy library or the L-BFGS quasi-Newton algorithm M1QN3 (from INRIA).
Local sensitivities are available by plotting the gradients (spatially distributed gradient values).


Assimilated data can be combinations of in-situ and remote-sensed observations, either time-series or spatially distributed.

A twin experiment mode (generation next inversions of synthetic data).  :ref:`is available in getting started section<2_make_your_first_4Dvar>`

.. Numerous benchmarks are available both for the forward models and the inverse ones. --> LATER


The Python wrapped versions can be easily interfaced or coupled with Python Neural Networks (deep learning) codes.

.. hint::

  **to be implemented**

  - WIP Assimilation of Lagrangian data e.g. extracted from video images.
  - WIP Superposition of local 2D models over 1D model, with simultaneous data assimilation process.

================================================
What you will get from this documentation ?
================================================

- Examples of Python scripts are presented throughout

**Installation**

- Python packages and software requirements: :ref:`install section in getting started<Installation>`.

**How to run basic direct and inverse runs from Python interface or using commandline**:

- How to run basic direct and inverse runs from Python
- Direct simulation: :ref:`1_make_your_first_run`
- Inverse run with variational assimilation using Fortran minimizer :ref:`2_make_your_first_4Dvar`
- Inverse run with variational assimilation using Python minimizer :ref:`3_inference python`
- Parallel computation :ref:`3_make_your_first_parallel`

**More in depth presentation of Python/DassFlow2D capabilities**:

- The methods relative to DassFlow package :ref:`api`
- Configuration and input file information  :ref:`input_data_guide`
- Output files information  :ref:`output_data_guide`
- How to pass information from Python to Fortran or from Fortran to python : presented via getters and setter thoughout :ref:`user_guide`


.. What to fill ? question de lilian ?
**Information on code architecture and philosophy:**
- Fortran architecture
- Wrapping architecture

==========================================
Main packages used
==========================================

The main external Python librairies used are:

- ``scipy`` for minimization algorithms (lbfgs). See their documentation here: https://docs.scipy.org/doc/scipy/
- ``pyvista`` for mesh plotting, and vtk interfacing. See their documentation here: https://docs.pyvista.org/
- ``vtk`` for vtk interfacing (is mainly used by ``pyvista``).
- ``numpy`` for classical array manipulation. See their documentation here:  https://numpy.org/doc/stable/
- ``matplotlib`` for classical plots. See heir documentation here:  https://matplotlib.org/stable/
- ``sphinx`` for documentation. See heir documentation here:  https://pydata-sphinx-theme.readthedocs.io/en/stable/index.html#
- ``rst``  for documentation. See heir documentation here:  https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html
  - https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst

Additionnaly librairies from base Python are used, such as:

- os : for console "communication". See documentation here:  https://docs.python.org/fr/3/library/os.html

==========================================
Acknowledgments
==========================================

The recent DassFlow results rely on a collaboration between IMT-INSA Toulouse (J. Monnier et al., math. modeling-comput. sc.), INRAE (P.-A. Garambois et al., hydrology modeling), CNES-CS group (K. Larnier et al., Senior Engineer, programming-assessing-datasets), with great PhD students and engineers: L. Pujol, L. Villenave, T. Malou, J. Verley and P. Brisset (IMT-INRAe-INSA-Univ. Strasbourg-CNES-CLS group). See publications list.
Moreover the HiVDI algorithm (K. Larnier, J. Monnier) aiming at estimating river discharges from SWOT data only (forthcoming space mission NASA-CNES et al.) is illustrated.
