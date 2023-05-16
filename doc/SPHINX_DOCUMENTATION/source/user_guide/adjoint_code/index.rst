.. _adjoint_code:

============
Adjoint code
============



In order to calculate the gradient of the cost function :math:`\nabla J` in a relative minimum of time computation,
and so to perform sensibility analysis **subsec:Sensibility-analysis** or data assimilation **subsec:4D-Var-data-assimilation**,
the ``dassflow2d`` software has been originally designed, in sense of source code structuration,
to generate automatically the discrete model adjoint using the differentiation tool **Tapenade**.

Some final tricks calling a Perl program ./src/adjoint/finish_to_gen_adjoint.pl have also been implemented to override some Tapenade manual operations
(dynamic array sizes, adjoint variables management and the adjoint of the MPI standard
communication operations) in order to really achieve to a complete automatic generation.

User can for example define a new cost function and so generate the associated discrete model adjoint.

+++++++++++++++++++++++++++++++++
Adjoint code generation
+++++++++++++++++++++++++++++++++

The automatic generation is called typing the command 'make install' at ``dassflow2d`` root directory. All discrete adjoint model source files are placed in the **code/src/tap/** directory.


• it calls make tap_files instance that generate differiated code

• In order to compile/link these source files, the ADJOINT variable must be set to '1' in the Makefile.inc.

• Then, typing the commane 'make' generates the executable \exe{} in the \bin{} directory.



+++++++++++++++++++++++++++++++++
Sensibility analysis
+++++++++++++++++++++++++++++++++


The generation of the discrete model adjoint provides a way to perform sensibility analysis and estimate the influence of the different parameters.
If :math:`\mathbf{k}` denotes the control vector containing all the parameters, then the absolute model sensibility to a given parameter :math:`k_{i}`,

:math:`{k_{i}}={\displaystyle \frac{\partial J}{\partial k_{i}}}`

are written in output result files providing a measure of the model change response due to a change of a given parameter.

This mode is called typing 'make rungrad' at the  ``dassflow2d`` root directory.

In context of the Shallow Water model, the control vector :math:`\mathbf{k}` and the corresponding output result files are,

• the hydrograph time series and its gradient in the 'hydographxxx_grad' files. The x characters denotes the label of the hydograph.

• the bed elevation scalar field :math:`z_{b}` and its gradient in the 'bathy_grad.yyy' file. The y characters correspond to the file format.

• the Manning-Strickler roughness coefficient n and its gradient in the 'manning_grad.yyy' file. The y characters correspond to the file format.

All these output result files are placed in the \bingrad{} directory at the end of the calcultation of the cost function gradient \nabla J.


+++++++++++++++++++++++++++++++++
4D-Var data assimilation
+++++++++++++++++++++++++++++++++


In this mode, a local minimum of the cost function is seeked using the discrete model adjoint to calculate the cost function gradient \nabla J and a local descent algorithm, i.e. m1qn3.
One can identify the model input control variables “matching at best” with observation data.

The identification process principle is sketched in :download:`figure downloadable here <files/daorganisation.pdf>`

Given a first guess :math:`\mathbf{k}_{0}` (user is free to configure his first guess as he would have configured a direct simulation),
we week the iterates :math:`\mathbf{k}_{i}` which make decrease the cost function using a descent algorithm. To do so, at each iterate,

1. The cost function :math:`J(\mathbf{k}_{i})` and its gradient :math:`\nabla J(\mathbf{k}_{i})` are computed calling the discrete forward model (from 0 to T) and its adjoint (from T to 0, reverse in time).

2. Given :math:`\mathbf{k}_{i}` , :math:`J(\mathbf{k}_{i})` and :math:`\nabla J(\mathbf{k}_{i})`, the m1qn3 library is invoked in order to find a new iterate such that :math:`J(\mathbf{k}_{i+1})<J(\mathbf{k}_{i})`.

3. The convergence criteria is tested, i.e., the ``eps_min`` parameter in the ``input.txt`` file.




This mode is called typing 'make runmin' at the  ``dassflow2d``  root directory. The iterate informations are written on screen and in the 'min_cost.txt' file in the **bin_A/min/** directory.

The control variables bulding the vectors math:`\mathbf{k}_{i}` are actived setting to '1' the following parameters in the ``input.txt`` file,

• <c_manning> for the Manning-Strickler roughness coefficient and written in the 'manning.xxx' output result files.

• <c_bathy> for the bed elevation and written in the bathy.xxx' output result files.

• <c_ic> for the initial condition and written in the 'ic.xxx' output result files.

• <c_hydrograph> for the hydrograph(s) and written in the 'hydograph_yyy.xxx' output result files.

• <c_ratcurve> for the rating curve(s) and written in the 'ratcurve_yyy.xxx' output result files.

• **TO BE COMPLETED**

The identification process could be very long in time and the ``restart_min`` parameter in the ``input.txt`` file can be setted to the maximum number of iterations to perform.
A file ``restart_min.bin`` is generated at each iteration and is read if it exists to restart the identification process.
