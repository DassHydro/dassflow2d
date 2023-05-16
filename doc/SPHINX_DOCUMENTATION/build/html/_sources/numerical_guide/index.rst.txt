.. _numerical_guide:

=================================
Numerical Guide
=================================
**Renommer en Numerical Guide, place l'onglet avant le dev guide dans le bandeau**

=================================
DassFLow capabilities at a glance
=================================

.. dropdown:: DassFlow capabilities at a glance
    :animate: fade-in-slide-down
    :container: + shadow
    :title: font-weight-bolder



    **Direct/forward model**

        The first forward model is based on the bidimensional Shallow Water equations written here with Manning-Strickler friction parameterization.

        On a computational domain :math:`\Omega\in\mathbb{\mathbb{\mathbb{R}}}^{2}` and for a time interval :math:`\left[0,T\right]`, the equations numerically resolved are:

        .. math::

            \left\{ \begin{array}{lclcc}
                \partial_{t}h+div(\mathbf{q}) & = & 0 & \textrm{in} & \text{Ω}\times]0,T]\\
                \\
                \partial_{t}\mathbf{q}+div\left(\dfrac{\mathbf{q}\otimes\mathbf{q}}{h}+g\dfrac{h^{2}}{2}\right) & = & -gh\mathbf{\nabla}z_{b}-g\dfrac{n^{2}\left\Vert \mathbf{q}\right\Vert }{h^{7/3}}\mathbf{q} & \textrm{in} & \text{Ω}\times]0,T]
                \end{array}\right.
            :label: SWE

        with provided initial and boundary conditions. 	Different types of boundary conditions can be prescribed at user-defined subsets of the computational domain to consider walls, inflows or outflows (see section **todo reference**).



        The state variables are the water depth h and the local discharge :math:`\mathbf{q}=h\mathbf{u}`, where :math:`\mathbf{u}=(u,v)^{T}` is the depth-averaged velocity vector. :math:`g` is the magnitude of the gravity, :math:`z_{b}` the bed elevation and :math:`n` the Manning-Strickler roughness coefficient.

        The model is numerically resolved by a Finite Volume method considering either structured or unstructured grids of discretization of the computational domain. (see doc **todo ref doc**).

        Several numerical schemes have been implemented and give the possibility to use a globally first or second order numerical solver with the well-balanced property.

        The global stability and robustness of the schemes insure a good treatment of dynamic wet/dry fronts without a water depth cut-off.

    **Adjoint model**

        The adjoint code can be used to perform a parametric sensitivity analysis or a parameter identification process (4D-Var data assimilation) by calculating the gradient with respect to parameters of a cost function :math:`J` measuring the discrepency between simulated and observed flow quantities. Only the first order scheme is available for the adjoint model.

        The adjoint code is automatically generated using the differentiation tool Tapenade [http://www-sop.inria.fr/tropics/] and some final tricks.

        *Gradient values : sensitivity analysis* - The computation of the cost gradient :math:`\nabla J` enables to perform spatially-distributed sensitivity analysis of the flow model to its input parameters; this is the 'local sensitivity mode'.

        *Variational data assimilation (4D-Var)* - Complete Variational Data Assimilation algorithm is available. This enables solving high dimensional optimization problems such as the identification of some input parameters values (eg. the bathymetry, some inflows or physical parametrizations) of the shalow water model from flow observations. The VDA process is based on local minimization algorithms adapted to minimize functions depending on large numbers of variables, such as the quasi-Newton technique L-BFGS method (**todo ref**). The optimization routine M1QN3 (**todo ref**) is directly implemented in the fortran source code but DassFLow wrapping also enables to use python optimization package instead (e.g. scipy).

    .. The original \DassFlow{} contained m1qn3 (written in Fortran), while the method lbfgs is introduced by using the wrapping and by applying the optimizer defined in the Python package scipy. --> MYTHO AT THE MOMENT


=================================
Full documentation
=================================

Here, you can find detailed documentation about the mathematical models and numerical solvers used for their resolution as well as about inverse algorithms. Furthermore, DassFlow website can be found here: <https://www.math.univ-toulouse.fr/DassFlow/>

.. _Math_num_doc:

----------------------------------------
Mathematical and numerical documentation
----------------------------------------

- :download:`doc-dassflow2d-num.lyx <./docs/doc-dassflow2d-num.lyx>`

--------------------------
Previous documentation
--------------------------

- :download:`doc_dassflow2d 2013 ("current") <./docs/doc-dassflow2d-current.pdf>`

- :download:`Tolosa-sw.pdf <./docs/Tolosa-sw.pdf>`


--------------------------
shenyuan documentation
--------------------------
- :download:`doc-shenyuan_DF2D-wrappe.lyx <./docs/doc-shenyuan_DF2D-wrappe.lyx>`
- :download:`doc-shenyuan_doc_developer.lyx <./docs/doc-shenyuan_doc_developer.lyx>`
- :download:`doc-shenyuan_doc_user.lyx <./docs/doc-shenyuan_doc_user.lyx>`
