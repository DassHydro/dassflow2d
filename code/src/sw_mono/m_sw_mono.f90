!======================================================================================================================!
!
!                    DassFlow Version 3.0
!
!======================================================================================================================!
!
!  Copyright University of Toulouse-INSA, Univ. of Strasbourg, INRAE & CNRS (France)
!
!  This file is part of the DassFlow software (Data Assimilation for Free Surface Flows).
!  DassFlow is a computational software whose purpose is to simulate geophysical free surface flows,
!  designed for variational sensitivities and data assimilation (4D-var). Inverse capabilities are
!  based on the adjoint code generation by a source-to-source algorithmic differentiation (Tapenade software used).
!
!  DassFlow software includes few mostly independent "modules" with common architectures and structures.
!  Please consult the DassFlow webpage for more details: http://www-gmm.insa-toulouse.fr/~monnier/DassFlow/.
!
!  Many people have contributed to the DassFlow development from the initial version to the latest ones.
!  Current contributions:
!               L. Pujol (PhD Unistra)
!               L. Villenave (PhD student)
!               P.-A. Garambois (INRAE Aix-en-Provence)
!               J. Monnier (INSA & Mathematics Institute of Toulouse IMT).
!               K. Larnier (CS group - IMT-INSA).
!  Former scientific or programming contributions of:
!               F. Couderc (CNRS & Mathematics Institute of Toulouse IMT).
!               J.-P. Vila (INSA & Mathematics Institute of Toulouse IMT).
!               R. Madec   (Mathematics Institute of Toulouse IMT).
!  plus less recent other developers (M. Honnorat and J. Marin).
!
!  Contact : see the DassFlow webpage
!
!  This software is governed by the CeCILL license under French law and abiding by the rules of distribution
!  of free software. You can use, modify and/or redistribute the software under the terms of the CeCILL license
!  as circulated by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
!
!  As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the
!  license, users are provided only with a limited warranty and the software's author, the holder of the economic
!  rights, and the successive licensors have only limited liability.
!
!  In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or
!  developing or reproducing the software by the user in light of its specific status of free software, that may
!  mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
!  experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
!  software's suitability as regards their requirements in conditions enabling the security of their systems and/or
!  data to be ensured and, more generally, to use and operate it in the same conditions as regards security.
!
!  The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you
!  accept its terms.
!
!======================================================================================================================!
!> \file m_sw_mono.f90
!! \brief This file includes the module m_model
!! \details The file includes only module m_model (see doc m_model module ).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module m_model
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> Module m_model
!!
!! \details  module specific to the use of DassFlow-2D to  solve  Shallow Water Equations
MODULE m_model

   USE m_common
   USE m_linear_algebra
   USE m_mesh
   USE m_mpi
   USE m_time_screen                                                                                              !NOADJ
!  USE m_user_data                      ! added because of friction (manning_user)

   implicit none

   integer(ip)  ::  sw_nb = 3    ! Usefull to create arrays dimensioned with the number of unknows

   !===================================================================================================================!
   !  Discrete Model Unknows Structure
   !===================================================================================================================!
    !> Discrete Model Unknowns Structure
    !!
    !! \details h,u,v are the unknowns of the model
    !! 2Dvectors (x,y directions ) of gradient of h,u,v are defined here also
   TYPE unk
   !> unknowns of the problem

      real(rp), dimension(:), allocatable  ::  h               !> water heigth (m)
      real(rp), dimension(:), allocatable  ::  u               !> speed along x (m.s-1)
      real(rp), dimension(:), allocatable  ::  v               !> speed along y (m.s-1)

      real(rp), dimension(:), allocatable  ::  infil               !> infiltrated water heigth (m)

!      real(rp), dimension(:), allocatable  ::  entropy               !> entropy for low froude scheme
      real(rp)     			   ::  t_display

      type(vec2d), dimension(:), allocatable  ::  grad_h        !> gradient of h (spatial at time t ?)
      type(vec2d), dimension(:), allocatable  ::  grad_u        !> gradient of u (spatial at time t ?)
      type(vec2d), dimension(:), allocatable  ::  grad_v        !> gradient of v (spatial at time t ?)
      type(vec2d), dimension(:), allocatable  ::  grad_z        !> gradient of bathymetry (spatial at time t ?)
!      type(vec2d), dimension(:), allocatable  ::  grad_entropy       !> gradient of ENTROPY (spatial at time t ?)

   END TYPE unk

! entropy for low froude scheme that should be added to dof


   !===================================================================================================================!
   !  Discrete Variables specific to Model
   !===================================================================================================================!

   real(rp), dimension(:), allocatable  ::  bathy_node         !> Bathymetry at mesh nodes
   real(rp), dimension(:), allocatable  ::  bathy_cell         !> Bathymetry at mesh cells gravity center
   real(rp), dimension(:), allocatable  ::  manning            !> Manning coefficient for mesh cells
   real(rp), dimension(:), allocatable  ::  manning_beta           !> Manning's beta (for power of h) coefficient for mesh cells

   integer(ip)  ::  nland                                      !> Total number of land associated to Manning

   integer(ip), dimension(:), allocatable  ::  land            !> Cells land number associated to Manning

   type(vec2d), dimension(:), allocatable  ::  grad_z        !> Cell Gradient of bathy_cell
   type(vec2d), dimension(:), allocatable  ::  grad_z2       !> Cell Gradient of bathy_cell^2
   type(vec2d), dimension(:), allocatable  ::  z_eq          !> Equivalent Bathymetry

   real(rp)  ::  mass_cut                                    !> ??????? mystery
   integer(ip)  ::  manning_data_glob                        !> ??????? mystery


   real(rp), dimension(:), allocatable :: slope_y !Added for Andromede, to expand
   real(rp), dimension(:), allocatable :: slope_x !Added for Andromede, to expand

   TYPE xsshp

    real(rp) :: xleft
    real(rp) :: xcenter
    real(rp) :: xright
    real(rp) :: s
    real(rp) :: hmax
    real(rp) :: topz

   END TYPE xsshp

   type(xsshp), dimension(:), allocatable :: XSshape

!    type(bathy_params) :: bathy_params

   !===================================================================================================================!
   !  Infiltration parameters Structure
   !===================================================================================================================!

   TYPE greenampt

    real(rp) :: PsiF              ! matric pressure
    real(rp) :: Ks                ! saturated hydraulic conductivity
    real(rp) :: DeltaTheta        ! variation of moisture content (saturated-initial)

   END TYPE

   TYPE scs_cn

    real(rp) :: lambdacn             ! initial abstraction ratio
    real(rp) :: CN                 ! Curve Number

   END TYPE

   TYPE infiltration_data

    integer(ip) :: nland

    integer(ip), dimension(:), allocatable  ::  land      ! Cells land number associated to infilration
    !real(rp)   , dimension(:), allocatable  ::  infil_qty     ! Infiltrated quantity at each cell NOW IN DOF%INFIL

    type(greenampt), dimension(:), allocatable :: GA
    type(scs_cn)   , dimension(:), allocatable :: SCS

    real(rp), dimension(:,:), allocatable :: coord

    real(rp), dimension(:), allocatable :: h_infil_max

   END TYPE

   type( infiltration_data ), target :: infil



   !===================================================================================================================!
   !  Friction parameters Structure
   !===================================================================================================================!

   TYPE friction_data
   !> derived type friction_data
    integer(ip)  ::  nland                                      !> Total number of land associated to Manning
    real(rp), dimension(:), allocatable  ::  manning            !> Manning coefficient each land (size is nland)
    real(rp), dimension(:), allocatable  ::  manning_beta       !> Manning's beta for each land (size is nland)

    integer(ip), dimension(:), allocatable  ::  land            !> nland value associated to cell k (land is ordered same as mesh)
   END TYPE friction_data


   ! bathy param structure

      TYPE param_model
   !> bathy_cell
      real(rp), dimension(:), allocatable  ::  bathy_cell               !> b

     end TYPE param_model

   !===================================================================================================================!
   !  Boundary Condition Structures
   !===================================================================================================================!

   integer(ip)  ::  feedback_inflow
   real(rp)  ::  coef_feedback


   TYPE gr4

      integer(ip)  ::  cell_id
      real(rp)  ::  surf
      real(rp), dimension(:), allocatable  ::  t , P , E , Q, P0, E0
      real(rp), dimension(4)  ::  params
      real(rp), dimension(13)  ::  state

   END TYPE

   !> Hydrograph definition
   TYPE hydrograph

      integer(ip)  ::  group
      real(rp), dimension(:), allocatable  ::  t,q

   END TYPE

	TYPE hpresc

	integer(ip)  ::  group
	real(rp), dimension(:), allocatable  ::  t,h

	END TYPE

	TYPE zspresc

	integer(ip)  ::  group
	real(rp), dimension(:), allocatable  ::  t,z

	END TYPE

   !> rating curve definition
   TYPE ratcurve

      integer(ip)  ::  group
      real(rp), dimension(:), allocatable  ::  h , q
      real(rp)  ::  z_rat_ref , zout , c1 , c2 , pow(2)

   END TYPE

   TYPE rain

      real(rp) :: x_min !Rain tile corners
      real(rp) :: x_max
      real(rp) :: y_min
      real(rp) :: y_max

      integer(ip) :: tile_index

      real(rp), dimension(:), allocatable :: t , q

      real(rp) :: qin       ! Current rain
      real(rp) :: cumul     ! Cumulated rain for SCS infiltration coefficient

   END TYPE

   !> Boundaries conditions structure
   TYPE bcs

      integer(ip)  ::  nb , nb_in , nb_out , nb_rn, nb_rn_t, nb_gr4in

      character(len=lchar), dimension(:,:), allocatable  ::  typ

      integer(ip), dimension(:), allocatable  ::   grpf

      real(rp), dimension(:), allocatable  ::  inflow
      real(rp), dimension(:), allocatable  ::  outflow

      type(hydrograph), dimension(:), allocatable::  hyd
      type(gr4), dimension(:), allocatable  ::  gr4
      type(ratcurve), dimension(:), allocatable::  rat
      type(hpresc), dimension(:), allocatable::  hpresc
      type(zspresc), dimension(:), allocatable::  zspresc

      type(rain), dimension(:), allocatable :: rain
      integer(ip), dimension(:), allocatable  ::  rain_land

      real(rp), dimension(:), allocatable  ::  sum_mass_flux

   END TYPE bcs
   !> definition of OBJECT boundary condition used along the code
   type(bcs), target  ::  bc

   !===================================================================================================================!
   !  Recording Structures
   !===================================================================================================================!
   !>  Mesurement Structures (to compare with simulation)
   TYPE station_obs

      type( point_in_mesh ), dimension(:), allocatable  ::  pt

      real(rp)  :: weight    ! weight of observations

      real(rp)  :: length    ! Length of river

      real(rp)  :: dt_offset ! Time of first observation

      real(rp)  :: dt        ! Frequency of observation ( satellite time repetitiveness)

      real(rp),dimension(:), allocatable :: dt_obs   ! Array with observation time

      integer(ip) :: ind_t   ! Index observation time

      integer(ip) :: nb_dt   ! Number of observation time

      real(rp), dimension(:), allocatable  ::  t , h , u , v , q, w

   END TYPE station_obs

   !>  observed section (tranche en travers pour calculer le débit au travers de cette tranche de rivière ?)
   TYPE section_obs

      type( point_in_mesh ), dimension(:), allocatable  ::  pt

      real(rp)  ::  dt , dx

      real(rp), dimension(:), allocatable  ::  t , h , u , v , q

      type( vec2d )  ::  normal

   END TYPE section_obs

   TYPE station_obsQ

      type( point_in_mesh ), dimension(:), allocatable  ::  pt

      real(rp)  :: weight    ! weight of observations

      real(rp)  :: length    ! Length of river

      real(rp)  :: dt_offset ! Time of first observation

      real(rp)  :: dt        ! Frequency of observation ( satellite time repetitiveness)

      real(rp),dimension(:), allocatable :: dt_obs   ! Array with observation time

      integer(ip) :: ind_t   ! Index observation time

      integer(ip) :: ind_bc   ! related BC index, as defined in bc.txt

      integer(ip) :: nb_dt   ! Number of observation time

      real(rp), dimension(:), allocatable  ::  t , h , u , v , Q, w

   END TYPE station_obsQ


   !===================================================================================================================!
   !  Recording Variables
   !===================================================================================================================!

   type(station_obs), dimension(:), allocatable  ::  station    !> definition stations obs
   type(section_obs), dimension(:), allocatable  ::  section    !> definition section (obs?)
   type(station_obsQ), dimension(:), allocatable ::  stationQ

   !===================================================================================================================!
   !  Input data and physical descriptors type
   !===================================================================================================================!


   TYPE soil_data

    real(rp)  ::  clay !>  Soil percentage of clay at each cell
    real(rp)  ::  silt !>  Soil percentage of silt at each cell
    real(rp)  ::  sand !>  Soil percentage of sand at each cell

    real(rp)  ::  ThetaS
    real(rp)  ::  ThetaR

   END TYPE soil_data

   TYPE ptf_data

    real(rp), dimension(9) :: Kappa !>  Coefficients for pedotransfer function

   END TYPE ptf_data

   TYPE surface_data

    real(rp)  ::  imperm  !>  Surface impermeabilisation percentage at each cell
!     integer(ip) :: imperm_group !>  Group of impermeabilization percentage derived from topographic maps pre-processed in Python

    real(rp)  ::  Dmax !>  Support rugosity at each cell
!     integer(ip) :: Dmax_group !>  Group of support rugosity

!     integer(ip)  ::  soil_occ_type  !>  Index of the occupation type

   END TYPE surface_data


   TYPE structure_data

!     integer(ip) ::  s_type !>  Index of the structure type

    character(len=lchar), dimension(:,:), allocatable :: name  !>  Name of the structure

    real(rp)  :: C1  !>  Coefficients for structure laws
    real(rp)  :: C2
    real(rp)  :: C3

    real(rp)  ::  true_x  !>  x-axis coordinate of the full structure
    real(rp)  ::  true_y  !>  y-axis coordinate of the full structure

   END TYPE structure_data


   TYPE input_data !> This structure should contain model parameters that are not meant to be inferred by VDA (physical descriptors)

    type(soil_data), dimension(:), allocatable  ::  soil      !>  Subsurface soil data
    integer(ip)  ::  soil_nland            !> Number of distinct lands associated to soil
    integer(ip), dimension(:), allocatable  ::  soil_land            !> Cells land number associated to soil

    type(ptf_data), dimension(:), allocatable  ::  ptf
    integer(ip)  ::  ptf_nland
    integer(ip), dimension(:), allocatable  ::  ptf_land

    type(surface_data), dimension(:), allocatable  ::  surf   !>  Surface data
    integer(ip)  ::  surf_nland            !> Number of distinct lands associated to surfaces
    integer(ip), dimension(:), allocatable  ::  surf_land            !> Cells land number associated to soil

    type(structure_data), dimension(:), allocatable  ::  structures !>  Hydraulic structures data
    integer(ip)  ::  struct_nland            !> Number of distinct lands associated to structures
    integer(ip), dimension(:), allocatable  ::  struct_land            !> Cells land number associated to soil

   END TYPE input_data

   type(input_data), target ::  phys_desc
   type(ptf_data), dimension(:), allocatable  ::  PTF

   !===================================================================================================================!
   !  Input variables specific to model (in addition to m_common)
   !===================================================================================================================!

   real(rp)     ::  g                                 !> Gravity constant
   real(rp)     ::  heps                              !> Cut-off of water depth to stabilize numerical scheme (if h<heps --> h=0)
   integer(ip)  ::  friction                          !> Activation of a Friction Law in Model

   ! Variables in control vector ( X )
   ! derivated by adjoint model ( grad(X) )
   integer(ip)  ::  c_shape_s
   integer(ip)  ::  c_xcenter
   integer(ip)  ::  c_hmax
   integer(ip)  ::  c_manning                         !> activate inference of manning alpha parameter
   integer(ip)  ::  c_manning_beta                    !> activate inference of manning_beta parameter
   integer(ip)  ::  c_bathy                           !> activate inference of bathymetry
   integer(ip)  ::  c_slope_y
   integer(ip)  ::  c_slope_x
   integer(ip)  ::  c_ic                              !> activate inference of ????
   integer(ip)  ::  c_hydrograph                      !> activate inference of hydrograph
   integer(ip)  ::  c_ratcurve                        !> activate inference of rating curve
   integer(ip)  ::  c_gr4params                       !> activate inference of all 4 gr4 parameters
   integer(ip)  ::  c_rain                            !> activate inference of rain
   integer(ip)  ::  c_infil_max
   integer(ip)  ::  c_Ks                             ! GA Infiltration parameter
   integer(ip)  ::  c_PsiF                           ! GA Infiltration parameter
   integer(ip)  ::  c_DeltaTheta                     ! GA Infiltration parameter
   integer(ip)  ::  c_lambda                         ! SCS-CN Infiltration parameter
   integer(ip)  ::  c_CN                             ! SCS-CN Infiltration parameter
   integer(ip)  ::  c_ptf                            ! Pedotransfer coefficients
    ! Each variable eps in perturbation control vector
    ! used to test validity of the adjoint model
   real(rp)     ::  eps_manning                       !
!   real(rp)     ::  eps_manning_beta                     ! TO ADD
   real(rp)     ::  eps_bathy                         !
   real(rp)     ::  eps_ic                            !
   real(rp)     ::  eps_hydrograph                    !
   real(rp)     ::  eps_ratcurve                      !
   real(rp)     ::  eps_gr4params                     !
   real(rp)     ::  eps_rain                         !
   real(rp)     ::  eps_Ks                           !
   real(rp)     ::  eps_PsiF                         !
   real(rp)     ::  eps_DeltaTheta                   !
   real(rp)     ::  eps_lambdacn                       !
   real(rp)     ::  eps_CN                           !
   real(rp)     ::  eps_ptf                           !

   real(rp)     ::  regul_manning                     !
!   real(rp)     ::  regul_manning_beta                     ! TO ADD
   real(rp)     ::  regul_bathy                       !
   integer(ip)     ::  regul_bathy_grad                       !
   integer(ip)     ::  regul_bathy_shape                       !
   real(rp)     ::  regul_ic                          !
   real(rp)     ::  regul_hydrograph                  !
   real(rp)     ::  regul_ratcurve                    !
   real(rp)     ::  regul_gr4params                   !

   integer(ip)  ::  fix_time_step_serie               !	affected in m_adjoint.f90 useful for only for the adjoint


   !===================================================================================================================!
   !  Input variables namelist (m_common + model specific upper ones)
   !===================================================================================================================!

   namelist/list_input/ &
      mesh_type, &
      mesh_name, &
      lx, &
      ly, &
      nx, &
      ny, &
      bc_N, &
      bc_S, &
      bc_W, &
      bc_E, &

      bc_rain, &
      bc_infil, &

      ts, &
      dt, &
      dtw, &
      dtp, &
      dta, &
      cfl, &
      adapt_dt, &

      do_warmup, &

      w_vtk, &
      w_tecplot, &
      w_gnuplot, &
      w_exact, &
      w_norm, &
      w_obs, &

      use_obs,&
      use_Zobs,&
      use_UVobs,&
      use_HUVobs,&
      use_Qobs,&
      use_Qobs_gr4,&
      use_NSE,&

      use_xsshp,&
      xsshp_along_x,&
      xsshp_along_y,&

      use_ptf,&

      spatial_scheme, &
      temp_scheme, &

      max_nt_for_direct , &
      max_nt_for_adjoint, &

      restart_min, &
      eps_min, &

      g, &
      heps, &
      friction, &
      feedback_inflow, &
      coef_feedback, &
      c_shape_s, &
      c_xcenter, &
      c_hmax, &
      c_manning, &
      c_manning_beta, &
      c_bathy, &
      c_slope_y,&
      c_slope_x,&
      c_ic, &
      c_hydrograph, &
      c_ratcurve, &
      c_gr4params, &
      c_rain, &
      c_infil_max, &
      c_Ks, &
      c_PsiF, &
      c_DeltaTheta, &
      c_lambda, &
      c_CN, &
      c_ptf,&

      eps_min, &
      eps_manning, &
      eps_bathy, &
      eps_ic, &
      eps_hydrograph, &
      eps_ratcurve, &
      eps_gr4params, &
      eps_rain, &
      eps_Ks, &
      eps_PsiF, &
      eps_DeltaTheta, &
      eps_lambdacn, &
      eps_CN, &
      eps_ptf, &

      regul_manning, &
      regul_bathy, &
      regul_bathy_grad, &
      regul_bathy_shape, &
      regul_ic, &
      regul_hydrograph, &
      regul_ratcurve, &
      regul_gr4params

     !> Input parameters (for wrap (BEWARE REPETITION WITH M_COMMOM.F90  + not used yet ?) )
 	TYPE Input_Param
 		character(len=lchar)  ::  mesh_type                        !> calling cartesian mesh or reader type ('dassflow' or 'basic'), 'gmsh' to finish
 		character(len=lchar)  ::  mesh_name                        !> mesh file name if not basic

 		character(len=lchar)  ::  bc_N                             !> (if mesh_type='basic') Type of boundary condition at North mesh boundary
 		character(len=lchar)  ::  bc_S                             !> (if mesh_type='basic') Type of boundary condition at South mesh boundary
 		character(len=lchar)  ::  bc_W                             !> (if mesh_type='basic') Type of boundary condition at West  mesh boundary
 		character(len=lchar)  ::  bc_E                             !> (if mesh_type='basic') Type of boundary condition at East  mesh boundary

 		integer(ip)  ::  bc_rain
 		integer(ip)  ::  bc_infil

 		real(rp)     ::  lx                                        !> (if mesh_type='basic') Lenght of computational domain x horizontal direction
 		real(rp)     ::  ly                                        !> (if mesh_type='basic') Lenght of computational domain y vertical   direction

 		integer(ip)  ::  nx                                        !> (if mesh_type='basic') Number of nodes in x horizontal direction
 		integer(ip)  ::  ny                                        !> (if mesh_type='basic') Number of nodes in y vertical   direction

 		real(rp)     ::  ts                                        !> User defined simulation time  (total time of simulation)

 		integer(ip)  ::  adapt_dt                                  !> Adaptative time step
 		real(rp)     ::  dt                                        !> (if adapt_dt=0) Time step value
 		real(rp)     ::  cfl                                       !> (if adapt_dt=1,2) CFL value

 		logical  ::  do_warmup                                      !> Toggle Warmup run for GR4 module

 		real(rp)     ::  dtw                                       !> Time step to Output Result Files  (write in external file)
		real(rp)     ::  dtp                                       !> Time step to Output Post Variables
 		real(rp)     ::  dta                                       !> Time Step to Generate BC (for Data Assimilation)
                                                                     !! (-- for linear interpolation of not defined input ? --- )

 		integer(ip)  ::  w_tecplot                                  !> Tecplot Output File (to check + update manning beta value)
 		integer(ip)  ::  w_vtk                                      !> VTK Output File     (to check + update manning beta value))
 		integer(ip)  ::  w_gnuplot                                  !> Gnuplot Output File (to check + update manning beta value))
 		integer(ip)  ::  w_bin                                      !> Binary Output File  (to check + update manning beta value))

 		integer(ip)  ::  w_exact                                    !> Exact Solution Output File
 		integer(ip)  ::  w_norm                                     !> Error Norms Calculation

 		integer(ip)  ::  w_obs                                     !< Gen Observation Output File
        integer(ip)  ::  use_obs                                   !< Use Observations in cost function definition
        integer(ip)  ::  use_Zobs                                  !< Write Water Surface elevation observations  and use them in cost function definition
        integer(ip)  ::  use_UVobs                                 !< Write Flow velocity observations  and use them in cost function definition
        integer(ip)  ::  use_HUVobs                                !< UNUSED Write At-a-cell-flow observations  and use them in cost function definition
        integer(ip)  ::  use_Qobs                                  !< Use Boundary flow observations in cost function definition
        integer(ip)  ::  use_Qobs_gr4                              !< Use Hydrological flow in cost function definition
        integer(ip)  ::  use_NSE                                   !< Use Nash-Sutcliffe Efficiency instead of RMSE for cost function definition using Flow observations

        integer(ip)  ::  use_xsshp                                  !< Use channel shape parameter "geometry_params.txt" file to parameterize cross-section and slope
        integer(ip)  ::  xsshp_along_x                              !< Toogle whether channel is defined along x-axis
        integer(ip)  ::  xsshp_along_y                              !< Toogle whether channel is defined along y-axis

        integer(ip)  ::  use_ptf                                    !< Toogle whether a pedotransfer function is used to calculate infil parameters from phys_desc parameters

		character(len=lchar)  ::  spatial_scheme                    !> Name of Spatial  Discretization Scheme ('first_b1' only at the moment)
 		character(len=lchar)  ::  temp_scheme                       !> Name of Temporal Discretization Scheme ('euler' or 'imex' at the moment )

 		character(len=lchar), dimension(:), allocatable  ::  args   !> Arguments passed on the command line

 		integer(ip)  ::  max_nt_for_direct                          !> Maximum iterations to perform the direct model
 		integer(ip)  ::  max_nt_for_adjoint                         !> Maximum iterations to perform the direct model in view
 																   !! to bound the memory of the adjoint model

 		real(rp)     ::  g                                 !> Gravity constant
 		real(rp)     ::  heps                              !> Cut-off of water depth to stabilize numerical scheme
 		integer(ip)  ::  friction                          !> Activation of a Friction Law in Model

    ! Variables in control vector ( X )
    ! derivated by adjoint model ( grad(X) )
        integer(ip)  ::  c_shape_s
        integer(ip)  ::  c_xcenter
        integer(ip)  ::  c_hmax
 		integer(ip)  ::  c_manning                         !> activate inference of manning alpha parameter (if c_xxx = 1)
 		integer(ip)  ::  c_manning_beta                    !> activate inference of manning beta parameter (if c_xxx = 1)
 		integer(ip)  ::  c_bathy                           !> activate inference of bathymetry (if c_xxx = 1)
 		integer(ip)  ::  c_ic                              !> activate inference of ???(if c_xxx = 1)
 		integer(ip)  ::  c_hydrograph                      !> activate inference of hydrograph r (if c_xxx = 1)
 		integer(ip)  ::  c_ratcurve                        !> activate inference of rating curve (if c_xxx = 1)
 		integer(ip)  ::  c_rain                            !> activate inference of rain r (if c_xxx = 1)
 		integer(ip)  ::  c_gr4params                       !> activate inference of all 4 gr4 parameters
 		integer(ip)  ::  c_infil_max
 		integer(ip)  ::  c_Ks                             ! GA Infiltration parameter
 		integer(ip)  ::  c_PsiF                           ! GA Infiltration parameter
 		integer(ip)  ::  c_DeltaTheta                     ! GA Infiltration parameter
 		integer(ip)  ::  c_lambda                         ! SCS-CN Infiltration parameter
 		integer(ip)  ::  c_CN                             ! SCS-CN Infiltration parameter
 		integer(ip)  ::  c_ptf                            ! Pedotransfer coefficients

     ! Each variable eps in perturbation control vector
     ! used to test validity of the adjoint model
 		real(rp)     ::  eps_min
 		real(rp)     ::  eps_manning                       !
 		real(rp)     ::  eps_bathy                         !
		real(rp)     ::  eps_ic                            !
 		real(rp)     ::  eps_hydrograph                    !
 		real(rp)     ::  eps_ratcurve                      !
 		real(rp)     ::  eps_rain                         !
 		real(rp)     ::  eps_Ks                           !
 		real(rp)     ::  eps_PsiF                         !
 		real(rp)     ::  eps_DeltaTheta                   !
 		real(rp)     ::  eps_lambdacn                       !
 		real(rp)     ::  eps_CN                           !
 		real(rp)     ::  eps_ptf                           !
     ! regularisation coeff ??? (weight given to each control variable?)
 		real(rp)     ::  regul_manning                     !
 		real(rp)     ::  regul_bathy                       !
 		integer(ip)     ::  regul_bathy_grad                       !
 		integer(ip)     ::  regul_bathy_shape                       !
		real(rp)     ::  regul_ic                          !
 		real(rp)     ::  regul_hydrograph                  !
 		real(rp)     ::  regul_ratcurve                    !
 	END TYPE Input_Param

CONTAINS

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Default values for Input variables namelist
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!
    !> Define default values   for Input variables namelist
   SUBROUTINE Default_values

      cfl       =  0.8_rp
      adapt_dt  =  1_ip

      dtw       =  2000.
      dtp       =  60.
      dt	=  5.

      w_tecplot =  0_ip
      w_vtk     =  0_ip
      w_gnuplot =  0_ip

      w_exact   =  0_ip
      w_norm    =  0_ip
      w_obs     =  0_ip

      use_obs   =  0_ip
      use_Zobs  =  0_ip
      use_UVobs =  0_ip
      use_HUVobs =  0_ip
      use_Qobs   =  0_ip
      use_Qobs_gr4   =  0_ip
      use_NSE   =  0_ip

      use_xsshp = 0_ip
      xsshp_along_x = 0_ip
      xsshp_along_y = 0_ip

      use_ptf = 0_ip

      do_warmup = .True.

      spatial_scheme  =  'first_b1'
      temp_scheme     =  'euler'

      max_nt_for_direct   =  100000000_ip
      max_nt_for_adjoint  =  2500_ip

      g         =  9.81_rp
      heps      =  0.00000001_rp
      friction  =  1_ip

      c_shape_s = 0_ip
      c_xcenter = 0_ip
      c_hmax = 0_ip
      c_manning     =  0_ip
      c_manning_beta   =  0_ip
      c_bathy       =  0_ip
      c_slope_y     =  0_ip
      c_slope_x     =  0_ip
      c_ic          =  0_ip
      c_hydrograph  =  0_ip
      c_ratcurve    =  0_ip
      c_gr4params   =  0_ip
      c_rain        =  0_ip
      c_infil_max   =  0_ip
      c_Ks          =  0_ip
      c_PsiF        =  0_ip
      c_DeltaTheta  =  0_ip
      c_lambda      =  0_ip
      c_CN          =  0_ip
      c_ptf         =  0_ip

      eps_manning     =  0.2_rp
      eps_bathy       =  0.01_rp
      eps_ic          =  0.1_rp
      eps_hydrograph  =  0.2_rp
      eps_ratcurve    =  0.1_rp
      eps_gr4params   =  0.1_rp
      eps_rain        =  0.2_rp
      eps_Ks            =  0.01_rp
      eps_PsiF          =  0.01_rp
      eps_DeltaTheta    =  0.01_rp
      eps_lambdacn        =  0.01_rp
      eps_CN            =  0.01_rp
      eps_ptf            =  0.01_rp

      regul_manning     =  0._rp
      regul_bathy       =  0._rp
      regul_bathy_grad       =  0_ip
      regul_bathy_shape      =  0_ip
      regul_ic          =  0._rp
      regul_hydrograph  =  0._rp
      regul_ratcurve    =  0._rp
      regul_gr4params   =  0._rp

      inquire( iolength = length_real ) tc

      tc0  =  0._rp
      nt0  =  0_ip

      feedback_inflow   =	1_ip
      coef_feedback     =	0.1_rp
      verbose  =  0_ip

      restart_min  =  0_ip
      eps_min      =  1.d-4

      fix_time_step_serie  =  0_ip

      is_file_open(:)  =  ''
      file_open_counter  =  0

   END SUBROUTINE Default_values



!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation of Model unk
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   subroutine alloc_dof(dof,mesh)
        !> Allocation of Model unk
        !> Notes
        !> -----
        !> **alloc_dof** :
        !>
        !> - Allocate dof.
        !>   ``common_array_allocation`` [mod_common_data.f90]
        !> - Run smash core.
        !>   ``smash_core`` [mod_smash_interface.f90]
        !> - Write results in ``.txt`` format.
        !>   ``write_smash_results`` [WRITE_RESULTS.f90]
        !>
        ! ============================= ===================================
        ! Parameters                    Description
        ! ============================= ===================================
        ! ``setup``                     model_setup Derived Type
        ! ``domain``                    mesh Derived Type
        ! ``watershed``                 catchments Derived Type
        ! ``inputdata``                 input_data Derived Type
        ! ``lois``                      loi_ouvrage Derived Type
        ! ``param``                     spatialparam Derived Type
        ! ``model_states``              spatialstates Derived Type
        ! ``model_routing_states``      spatiotemporalstates Derived Type
        ! ``outputs``                   smash_outputs Derived Type
        ! ============================= ===================================
      implicit none
      type(msh), intent(in)  ::  mesh
      type(unk), intent(out)  ::  dof
      allocate(dof%h(mesh%nc + mesh%ncb))
      allocate(dof%u(mesh%nc + mesh%ncb))
      allocate(dof%v(mesh%nc + mesh%ncb))

      allocate(dof%infil(mesh%nc + mesh%ncb))

 !     allocate(dof%entropy(mesh%nc )) ! entropy for low froude scheme

      allocate(dof%grad_h(mesh%nc + mesh%ncb))
      allocate(dof%grad_u(mesh%nc + mesh%ncb))
      allocate(dof%grad_v(mesh%nc + mesh%ncb))
      allocate(dof%grad_z(mesh%nc + mesh%ncb))
!      allocate(dof%grad_entropy(mesh%nc )) ! entropy for low froude scheme

      dof%h(:)  =  0._rp
      dof%u(:)  =  0._rp
      dof%v(:)  =  0._rp

      dof%infil(:)  =  0._rp
 !     dof%entropy(:) = 0._rp ! entropy for low froude scheme

      dof%grad_h(:)%x  =  0._rp
      dof%grad_h(:)%y  =  0._rp

      dof%grad_u(:)%x  =  0._rp
      dof%grad_u(:)%y  =  0._rp

      dof%grad_v(:)%x  =  0._rp
      dof%grad_v(:)%y  =  0._rp

      dof%grad_z(:)%x  =  0._rp
      dof%grad_z(:)%y  =  0._rp

 !     dof%grad_entropy(:)%x  =  0._rp
 !     dof%grad_entropy(:)%y  =  0._rp
   end subroutine




!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Deallocation of Model unk
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   SUBROUTINE dealloc_dof(dof)
      implicit none
      type(unk), intent(inout)  ::  dof
      if (allocated(dof%h)) deallocate(dof%h)
      if (allocated(dof%u)) deallocate(dof%u)
      if (allocated(dof%v)) deallocate(dof%v)

      if (allocated(dof%infil)) deallocate(dof%infil)
!      if (allocated(dof%entropy)) deallocate(dof%entropy)  ! entropy for low froude scheme
      if (allocated(dof%grad_h)) deallocate(dof%grad_h)
      if (allocated(dof%grad_u)) deallocate(dof%grad_u)
      if (allocated(dof%grad_v)) deallocate(dof%grad_v)
      if (allocated(dof%grad_z)) deallocate(dof%grad_z)
 !     if (allocated(dof%grad_entropy)) deallocate(dof%grad_entropy)     ! entropy for low froude scheme
   END SUBROUTINE


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Deallocation of Model unk
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE dealloc_model

      implicit none

      if ( allocated( bathy_node ) ) 		deallocate( bathy_node )
      if ( allocated( bathy_cell ) ) 		deallocate( bathy_cell )
      if ( allocated( XSshape ) ) 		    deallocate( XSshape )
      if ( allocated( slope_y ) ) 		    deallocate( slope_y )
      if ( allocated( slope_x ) ) 		    deallocate( slope_x )
      if ( allocated( land       ) ) 		deallocate( land       )
      if ( allocated( manning    ) ) 		deallocate( manning    )
      if ( allocated( manning_beta ) ) 		deallocate( manning_beta )
      if ( allocated( infil%land ) ) 		deallocate( infil%land )
      if ( allocated( infil%coord ) ) 		deallocate( infil%coord )
      if ( allocated( infil%GA   ) ) 		deallocate( infil%GA   )
      if ( allocated( infil%SCS  ) ) 		deallocate( infil%SCS  )
      if ( allocated( infil%h_infil_max ) ) 		deallocate( infil%h_infil_max )
      if ( allocated( phys_desc%soil_land ) ) 		 deallocate( phys_desc%soil_land )
      if ( allocated( phys_desc%soil ) ) 		     deallocate( phys_desc%soil )
      if ( allocated( phys_desc%ptf_land ) ) 		 deallocate( phys_desc%ptf_land )
      if ( allocated( phys_desc%ptf ) ) 		     deallocate( phys_desc%ptf )
      if ( allocated( PTF ) ) 		     deallocate( PTF )
!       if ( allocated( phys_desc%surf_land ) ) 		 deallocate( phys_desc%surf_land )
!       if ( allocated( phys_desc%surf ) ) 		     deallocate( phys_desc%surf )
!       if ( allocated( phys_desc%struct_land ) ) 	 deallocate( phys_desc%struct_land )
!       if ( allocated( phys_desc%structures ) ) 		 deallocate( phys_desc%structures )

      !------------------------------------------------!
      ! millascenious forgoten variables to deallocate
      !------------------------------------------------!

      ! BC
      if ( allocated( bc%typ ) ) 			deallocate(bc%typ)
      if ( allocated( bc%grpf ) ) 			deallocate(bc%grpf)
      if ( allocated( bc%inflow ) ) 		deallocate(bc%inflow)
      if ( allocated( bc%outflow ) ) 		deallocate(bc%outflow)
      if ( allocated( bc%hyd ) ) 			deallocate(bc%hyd)
      if ( allocated( bc%hpresc ) ) 		deallocate(bc%hpresc)
      if ( allocated( bc%zspresc ) ) 		deallocate(bc%zspresc)
      if ( allocated( bc%rain ) ) 			deallocate(bc%rain)
      if ( allocated( bc%rat ) ) 			deallocate(bc%rat)
      if ( allocated( bc%sum_mass_flux ) ) 	deallocate(bc%sum_mass_flux)
      if ( allocated( bc%rain_land ) ) 	deallocate(bc%rain_land)
      if ( allocated( bc%rain ) ) 	deallocate(bc%rain)


      if ( allocated( swap_index ) ) deallocate(swap_index)
      if ( allocated( inv_swap_index ) ) deallocate(inv_swap_index)
      if ( allocated( part ) ) deallocate(part)
      if ( allocated( part_size ) ) deallocate(part_size)
      if ( allocated( part_neighbs ) ) deallocate(part_neighbs)

      ! initialization.f90 things
		! zspresc
      if ( allocated( grad_z ) ) deallocate(grad_z)
      if ( allocated( grad_z2 ) ) deallocate(grad_z2)
      if ( allocated( z_eq ) ) deallocate(z_eq)

		! minimization
      !if ( allocated( innovation ) ) deallocate(innovation)
      !if ( allocated( innovW ) ) deallocate(innovW)


		if ( allocated( station ) ) deallocate(station)
		if ( allocated( stationQ ) ) deallocate(stationQ)
		if ( allocated( section ) ) deallocate(section)
		! measure stations
!~       if ( allocated( station%pt ) ) deallocate(station%pt)
!~       if ( allocated( station%dt_obs ) ) deallocate(station%dt_obs)
!~       if ( allocated( station%t ) ) deallocate(station%t)
!~       if ( allocated( station%h ) ) deallocate(station%h)
!~       if ( allocated( station%u ) ) deallocate(station%u)
!~       if ( allocated( station%v ) ) deallocate(station%v)
!~       if ( allocated( station%q ) ) deallocate(station%q)
!~       if ( allocated( station%w ) ) deallocate(station%w)

!~ 		! measure sections
!~       if ( allocated( section%pt ) ) deallocate(section%pt)
!~       if ( allocated( section%t ) ) deallocate(section%t)
!~       if ( allocated( section%h ) ) deallocate(section%h)
!~       if ( allocated( section%u ) ) deallocate(section%u)
!~       if ( allocated( section%v ) ) deallocate(section%v)
!~       if ( allocated( section%q ) ) deallocate(section%q)
!~       if ( allocated( section%w ) ) deallocate(section%w)

   END SUBROUTINE dealloc_model



!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calling MPI and communicating dof to fill ghost cells
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


   SUBROUTINE com_dof( dof , mesh )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in   )  ::  mesh
      type( unk ), intent(inout)  ::  dof

      !================================================================================================================!
      !
      !================================================================================================================!

      #if defined USE_MPI || USE_MPI_ADJ

         call Time_Init_Part(80)                                                                                  !NOADJ

         call com_var_r( dof%h(:) , mesh )
         call com_var_r( dof%u(:) , mesh )
         call com_var_r( dof%v(:) , mesh )

         if (bc_infil .ne. 0) call com_var_r( dof%infil(:) , mesh )

         call Time_End_Part(80)                                                                                   !NOADJ

      #endif

   END SUBROUTINE com_dof


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Allocation or Reallocate Stations Type
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   SUBROUTINE alloc_or_realloc_station( station_inout , new )

      implicit none

      type( station_obs ), dimension(:), allocatable, intent(inout)  ::  station_inout

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old , iobs , pt

      type( station_obs ), dimension(:), allocatable  ::  station_tmp

      intrinsic move_alloc

      if ( .not. allocated( station_inout ) ) then

         allocate( station_inout( new ) )

         return

      end if

      old = size( station_inout )

      if     ( new == old ) then

         return

      else if ( new > old ) then

         allocate( station_tmp( new ) )

         do iobs = 1,old

            station_tmp( iobs )%dt      =  station_inout( iobs )%dt
           station_tmp( iobs )%weight  =  station_inout( iobs )%weight

            allocate( station_tmp( iobs )%pt( size(station_inout( iobs )%pt ) ) )

            do pt = 1,size( station_inout( iobs )%pt )

               station_tmp( iobs )%pt( pt )%cell   =  station_inout( iobs )%pt( pt )%cell
               station_tmp( iobs )%pt( pt )%coord  =  station_inout( iobs )%pt( pt )%coord

            end do

         end do

         call move_alloc( station_tmp , station_inout )

      else

         call Stopping_Program_Sub( 'Wrong Station Dimension for Allocation' )

      end if

   END SUBROUTINE alloc_or_realloc_station

   SUBROUTINE alloc_or_realloc_stationQ( station_inout , new )

      implicit none

      type( station_obsQ ), dimension(:), allocatable, intent(inout)  ::  station_inout

      integer(ip), intent(in)  ::  new

      integer(ip)  ::  old , iobs , pt

      type( station_obsQ ), dimension(:), allocatable  ::  station_tmp

      intrinsic move_alloc

      if ( .not. allocated( station_inout ) ) then

         allocate( station_inout( new ) )

         return

      end if

      old = size( station_inout )

      if     ( new == old ) then

         return

      else if ( new > old ) then

         allocate( station_tmp( new ) )

         do iobs = 1,old

            station_tmp( iobs )%dt      =  station_inout( iobs )%dt
           station_tmp( iobs )%weight  =  station_inout( iobs )%weight

            allocate( station_tmp( iobs )%pt( size(station_inout( iobs )%pt ) ) )

            do pt = 1,size( station_inout( iobs )%pt )

               station_tmp( iobs )%pt( pt )%cell   =  station_inout( iobs )%pt( pt )%cell
               station_tmp( iobs )%pt( pt )%coord  =  station_inout( iobs )%pt( pt )%coord

            end do

         end do

         call move_alloc( station_tmp , station_inout )

      else

         call Stopping_Program_Sub( 'Wrong Station Dimension for Allocation' )

      end if

   END SUBROUTINE alloc_or_realloc_stationQ

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Routines for wrapping
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

																													!<NOADJ
	!==================================================================================================================!
	! initialise and finalize  dof
	!==================================================================================================================!
     subroutine unk_initialise(dof,mesh)
      implicit none
      type(msh), intent(in)  ::  mesh
      type(unk), intent(out)  ::  dof
      allocate(dof%h(mesh%nc + mesh%ncb))
      allocate(dof%u( mesh%nc + mesh%ncb))
      allocate(dof%v(mesh%nc + mesh%ncb))

      allocate(dof%infil(mesh%nc + mesh%ncb))
!      allocate(dof%entropy(mesh%nc )) ! entropy for low froude scheme

      allocate(dof%grad_h(mesh%nc + mesh%ncb))
      allocate(dof%grad_u(mesh%nc + mesh%ncb))
      allocate(dof%grad_v(mesh%nc + mesh%ncb))
      allocate(dof%grad_z(mesh%nc + mesh%ncb))
!      allocate(dof%grad_entropy(mesh%nc )) ! entropy for low froude scheme

      dof%h(:)  =  0._rp
      dof%u(:)  =  0._rp
      dof%v(:)  =  0._rp

      dof%infil(:)  =  0._rp
 !     dof%entropy(:)  =  0._rp            ! entropy for low froude scheme

      dof%grad_h(:)%x  =  0._rp
      dof%grad_h(:)%y  =  0._rp

      dof%grad_u(:)%x  =  0._rp
      dof%grad_u(:)%y  =  0._rp

      dof%grad_v(:)%x  =  0._rp
      dof%grad_v(:)%y  =  0._rp

      dof%grad_z(:)%x  =  0._rp
      dof%grad_z(:)%y  =  0._rp

!      dof%grad_entropy(:)%x  =  0._rp  ! entropy for low froude scheme
!      dof%grad_entropy(:)%y  =  0._rp    ! entropy for low froude scheme
   end subroutine

	! destroy a dof dof
	SUBROUTINE unk_finalise(dof)
      implicit none
      type(unk), intent(inout)  ::  dof
      call mpi_wait_all
!       write(*,*) proc, "(allocated(dof%h)) deallocate(dof%h)", allocated(dof%h), size(dof%h)
      if (allocated(dof%h)) deallocate(dof%h)
!       write(*,*) proc, "(allocated(dof%h)) deallocate(dof%u)", allocated(dof%u), size(dof%u)
      if (allocated(dof%u)) deallocate(dof%u)
!       write(*,*) proc, "(allocated(dof%h)) deallocate(dof%v)", allocated(dof%v), size(dof%v)
      if (allocated(dof%v)) deallocate(dof%v)
      call mpi_wait_all
! write(*,*) proc, "(allocated(dof%h)) deallocate(dof%infil)", allocated(dof%infil), size(dof%infil)
      if (allocated(dof%infil)) deallocate(dof%infil)
!       write(*,*) proc, "(allocated(dof%h)) deallocate(dof%grad_h)", allocated(dof%grad_h), size(dof%grad_h)
!      if (allocated(dof%entropy)) deallocate(dof%entropy)   ! entropy for low froude scheme
      if (allocated(dof%grad_h)) deallocate(dof%grad_h)
      call mpi_wait_all
!       write(*,*) proc, "(allocated(dof%h)) deallocate(dof%grad_u)", allocated(dof%grad_u), size(dof%grad_u)
      if (allocated(dof%grad_u)) deallocate(dof%grad_u)
!       write(*,*) proc, "(allocated(dof%h)) deallocate(dof%grad_v)", allocated(dof%grad_v), size(dof%grad_v)
      if (allocated(dof%grad_v)) deallocate(dof%grad_v)
!       write(*,*) proc, "(allocated(dof%h)) deallocate(dof%grad_z)", allocated(dof%grad_z), size(dof%grad_z)
      if (allocated(dof%grad_z)) deallocate(dof%grad_z)
 !     if (allocated(dof%grad_entropy)) deallocate(dof%grad_entropy)   ! entropy for low froude scheme
   END SUBROUTINE




!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Get cell spatial index from tile coordinates, for land attribution
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

   !>  Get cell spatial index from square tile coordinates
   !!
   !! \details From a file name and a extension ('tecplot', 'gnuplot', ...) complete file name.
   !! \param[in]  mesh Mesh
   !! \param[in]  xmin, xmax, ymin, ymax Coordinates of the tile
   !! \param[in]  spatial_index Cell index
   !! \return The spatial index of the cell.

   SUBROUTINE spatial_index_fromxy(mesh, xmin, xmax, ymin, ymax, spatial_index)

    type(msh), intent(  in) :: mesh
    real(rp),  intent(  in) :: xmin, xmax, ymin, ymax
    integer(ip), intent( out) :: spatial_index

    spatial_index = 0_ip

    if ((mesh%cell(i)%grav%x > xmin) .and. (mesh%cell(i)%grav%x < xmax) .and. &
        (mesh%cell(i)%grav%y > ymin) .and. (mesh%cell(i)%grav%y < ymax)) then

        spatial_index = k

    endif

   END SUBROUTINE spatial_index_fromxy


																													!>NOADJ

END MODULE m_model
