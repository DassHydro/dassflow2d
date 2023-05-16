.. _Coding Philosophy:

This short section introduces notions of code requirements pertaining to the Fortran core code and its interfacing with Python.

===================
Coding Philosophy
===================
------------------------------------------------------------------------------------
Variational assimilation through adjoint-based method and automatic differentiation
------------------------------------------------------------------------------------

In this code, Variational Data Assimilation is used as an adapted approach to deal with non-linear dynamic models and with heterogeneous observations in space and time, but also to deal with multivariate data assimilation problems and large control vectors.

It relies on the gradient of a cost function (e.g. misfit to observations). This gradient is obtained by applying the automatic differentiation tool Tapenade to Fortran source codes found in code/src. The resulting files are called the adjoint code and appear in code/tap. The file names correspond to several routines found in code/src files. In the source code, the **NOADJ** tag indicates that the line or lines are intended to be skipped during this process.

For the total code to be differentiable, the mathematical formula, or combination thereof, linking the TAP_INVARS and TAP_OUTVARS structures list in the Makefile (not .inc), needs to be differentiable.

As a rule of thumb, any part of the code that does not affect this link between variables and parameters may be flagged as **NOADJ**. For example, many computational checks or result writing routines are currently flagged as **NOADJ**.

Errors in the differentiation process may or may not cause it to abort. Instead, errors from within /tap files may be raised during a run (direct or inverse). In this case, it is preferable to address the issue from the corresponding /src file or routine than to modify the /tap files.

Any parameter that is not linked to the **cost** variable and written/read in the **control** will not be inferable, although it may not forbid differentiation and inverse modeling.


------------------------------------------------------------------------------------
Passing structures from the computational core to the corresponding Python classes
------------------------------------------------------------------------------------

Python can be used to pre-process and post-process data for DassFlow. The Fortran core will, at the minimum, handle the gradient computation, which benefits from the compiled code's efficiency, and may also handle gradient minimization (see `2_make_your_first_4Dvar`) as well as data handling. Python scripts can replace th Fortran code in the latter two of these tasks (see `3_inference python`).

Currently, the Fortran code handles most of the reading of input data. If Python is used to read data and passes it to Fortran, it will overwrite potentially existing data read by Fortran earlier in the process. Thus, additional input data may be read using Python without modifying the Fortran code.

You may need to pass the following (non-exhaustive) list of structures:

- Initial hydraulic states
dof0%h ,
dof0%u ,
dof0%v

- Current hydraulic states
dof%h ,
dof%u ,
dof%v

- Model parameters
manning ,
manning_beta ,
bathy_cell ,
infil%GA%Ks ,
infil%GA%PsiF ,
infil%GA%DeltaTheta ,
infil%SCS%lambda ,
infil%SCS%CN

- Boundary condition definition data
bc%hyd%t ,
bc%hyd%q ,
bc%rat%h ,
bc%rat%q ,
bc%rain%t ,
bc%rain%q

- Inverse modeling structures
control,
control_back,
cost

- MPI parameters
np (number of CPU cores)

Note that the order in which you pass data through the interface is important.
