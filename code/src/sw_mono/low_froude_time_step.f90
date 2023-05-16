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
!\brief  Perform Low froude Time Step dedicated to Shallow-Water Equations
!\return dof updated after this new timestep
SUBROUTINE low_froude_time_step( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_time_screen                                                                                              !NOADJ
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in   )  ::  mesh

   type( unk ), intent(inout)  ::  dof

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!
   real(rp), dimension( sw_nb , mesh%nc )  ::  tflux     !> Finite Volume total flux for each cell

   integer(ip)  ::  iL , iR                              !> Left and Right cells indexes to edge


   ! ======================================== LEVEL 1 -common-  Variables ======================================== !
   integer(ip)  ::  surfL , surfR                              !> Left and Right SURFACES of cells  (for low_froude scheme)
   integer(ip)  ::  periL , periR                              !> Left and Right PERIMETERS of cells  (for low_froude scheme)
   integer(ip)  :: edge_length                                 !> edge length

   ! DOF
   real(rp)  ::  hL , uL(2) , vL(2)             !> Left  State in edge cell normal coordinates
   real(rp)  ::  hR , uR(2) , vR(2)             !> Right State in edge cell normal coordinates, zR correspond to bathymetry
   ! The first element of each variable correspond to its value in global coordinate system
   !     second elemment correpond to its value in local coordinate systeme (relative to the edge, for cell k to cell ke)
   ! The second elements are used to define necessary input values for calc_boundary_state subroutine


   ! ZR ZL
   real(rp)  :: zR                          !> Bathymetry Left cell
   real(rp)  :: zL                          !> Bathymetry Right cell

   ! weight rusanob
   real(rp)  :: weight_rusanov             !> weight given to rusanov scheme, for one edge treatment (compare to entropy dissipative part)
   real(rp), dimension(:), allocatable  :: all_weight_rusanov

   ! normal(x,y) is the normal vector of the edge ie, from cell k to cell ke
   real(rp)  :: normal_x    !
   real(rp)  :: normal_y    !

   ! g (acceleration pensenteur) is defined in a module
   real(rp)  :: phiL       ! pressure +gravity potential on cell k
   real(rp)  :: phiR       ! pressure +gravity potential on cell ke

  ! ======================================== LEVEL 2 -specific-  Variables ======================================== !
   real(rp)  ::  hL_rusanov             !> reconstructed water height from andusse et al (Left)
   real(rp)  ::  hR_rusanov             !> reconstructed water height from andusse et al (Right)

   real(rp)  :: phiL_rusanov       !> reconstructed pressure potential on cell k
   real(rp)  :: phiR_rusanov       !> reconstructed pressure potential on cell ke

   real(rp)  :: lambda_rusanov       !> wave speed for rusanov scheme

   real(rp)  :: alpha_lf       !> weigth in entropy dissipative scheme
   real(rp)  :: gamma_lf       !> weigth in entropy dissipative scheme

  ! ======================================== LEVEL 3 -FLUXES-  Variables ======================================== !
 ! >>> For rusanov scheme
   real(rp)  ::  F_rusanov      !> F part of rusanov flux
   real(rp)  ::  G_rusanov(2)   !> G part of rusanov flux, first element correpond to x-axis part, second to y axis part
   real(rp)  ::  phi_rusanov    !> Phi part of rusanov flux
 ! >>> For entropy disspative scheme
   real(rp)  ::  F_lf      !> F part of entropy dissipative       flux
   real(rp)  ::  G_lf(2)   !> G part of entropy dissipative       flux, first element correpond to x-axis part, second to y axis part
   real(rp)  ::  phi_lf    !> Phi part of entropy dissipative       flux

   ! >>> Temporary variable to store combination of both flux defined above
   real(rp)  ::  ftot_edge(3) !> Temporary variable to store combination of both flux defined above

   integer(rp)       :: tmp_index_cells_1(4)
   character(rp)  :: tmp_index_cells_2(4)


  ! ======================================== LEVEL 0 -temporary to store intermediate values -  Variables ======================================== !
   real(rp)  :: c1       !> temporary wave speed for rusanov scheme
   real(rp)  :: c2       !> temporary wave speed for rusanov scheme
   real(rp)  :: Fmax_lf, Fmin_lf !> temporary variable for G flux calculation in LF scheme

   real(rp)  :: h_cell    !>       tmp h value on one cell, used to interpolate weight rusanov coefficient
   real(rp)  :: tmp      !> temporary variable to make tests
    ! ======================================== UNUSED ?  ======================================== !
   real(rp)  ::  h , u , v                               !> Temporal primitive variables
   real(rp)  ::  vel                                     !> Velocity norm
   real(rp)  ::  sfl                                     !> Manning

   real(rp), dimension( sw_nb )  ::  nflux               !> Finite Volume normal edge flux
   real(rp), dimension( sw_nb )  ::  lflux               !> Finite Volume edge flux in (x,y) coordinates
   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!


   !-------------------------------------------------------------------------------------------------------------------!
   !  Initialise variables (that can be initialised)
   !-------------------------------------------------------------------------------------------------------------------!

   tflux(:,:)  =  0._rp


   !force weight of scheme repartition: a value is estimated on each cell (also ghost cells)
   allocate(all_weight_rusanov(mesh%nc + mesh%ncb) )
   all_weight_rusanov(:) = 0

   do i = 1,mesh%nc
    if (mesh%cell(i)%boundary) then ! In practice, \delta_{k} is set to 1 for cells with an edge at the physical boundary ..
        all_weight_rusanov(i) = 1
    else if ( dof%h(i) < heps ) then ! \delta_{k} is set to 1 for cells for which there is a neighbouring one with a dry state
        all_weight_rusanov(i) = 1
    else
        all_weight_rusanov(i) = 1 ! to interpolate :  then going from cell to cell with common edges from a given number of (typically from 3 to 5), the value goes to 0 linearly.
    end if
   end do

 ! force alpha and gamma parameters
   alpha_lf = demi
   gamma_lf = demi

!   write(*,*) "starting loop on edges "
!   write(*,*) "mesh%ne", mesh%ne
   open(20,file="flux_edge.txt",status='replace',form='formatted')
   write(20,*) "ie, iL, iR, mesh%edge(ie)%boundary, mesh%edge(ie)%subdomain, Ftot_edge(:), hL * (( 1 - weight_rusanov) * phi_lf + weight_rusanov * phi_rusanov)"
   open(30,file="phi.txt",status='replace',form='formatted')
   write(30,*) "phiL, phiR, phiL_rusanov, phiR_rusanov"

   open(40,file="flux_rusanov.txt",status='replace',form='formatted')
   write(40,*) "F, G(1), G(2), phi_rusanov"
   open(50,file="flux_lf.txt",status='replace',form='formatted')
   write(50,*) "F, G(1), G(2), phi_lf"

   !-------------------------------------------------------------------------------------------------------------------!
   !  Calculate fluxes on each edge
   !-------------------------------------------------------------------------------------------------------------------!

   do ie = 1,mesh%ne
    
   !================================================================================================================!
   !  Calculate Left and Right States
   !================================================================================================================!

      ! identify cell indices
      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)


   !-------------------------------------------------------------------------------------------------------------------!
   !  get water height values
   !-------------------------------------------------------------------------------------------------------------------!
      ! Store h values for the two cells neibouring the edge number ie.
      hL  =  dof%h( iL )
      hR  =  dof%h( iR )

   !-------------------------------------------------------------------------------------------------------------------!
   !  get geometrical properties
   !-------------------------------------------------------------------------------------------------------------------!

       ! Store surface and permieter of both cells
      surfL = mesh%cell(iL)%surf
      periL = mesh%cell(iL)%peri

      if(iR<mesh%nc) then
         ! if this is not a ghost cell, extract directly information
         surfR = mesh%cell(iR)%surf
         periR = mesh%cell(iR)%peri
      else
         ! the ghost cell geometry is the same as its associated boundary cell (WARNING supposition lilian)
         surfR = mesh%cell(iL)%surf
         periR = mesh%cell(iL)%peri
      end if

      ! get normal vector values for this edge
      normal_x = mesh%edge(ie)%normal%x
      normal_y = mesh%edge(ie)%normal%y
      ! get edge length
      edge_length = mesh%edge(ie)%length


   !-------------------------------------------------------------------------------------------------------------------!
   !  If there is water (for one of the two edges), estimate fluxes on the edge.
   !  Else, the flux on this edge is null  -> go to "Euler Time Step" section directly.
   !-------------------------------------------------------------------------------------------------------------------!

if ( hL > heps .or. hR > heps ) then
      ! store bathymetry values
      zL  =  bathy_cell( iL )
      zR  =  bathy_cell( iR )

   !-------------------------------------------------------------------------------------------------------------------!
   !  get flow speeds and their reprojection on the edge (necessary for boundary treatment)
   !-------------------------------------------------------------------------------------------------------------------!
   ! >>> Left states
      uL(1)  =  dof%u( iL )
      vL(1)  =  dof%v( iL )

      uL(2)  =  normal_x * uL(1) + normal_y * vL(1)
      vL(2)  =  normal_x * vL(1) - normal_y * uL(1)

      if ( mesh%edge(ie)%boundary )  then
      ! ---------- For ghost cells ---------- !
      	call calc_boundary_state( mesh , hL , zL , uL(2) , vL(2) , &
                                 hR , zR , uR(2) , vR(2) )
      ! reproject in global coordinate system
      	uR(1) =  normal_x * uR(2)   -  normal_y * vR(2)
      	vR(1) =  normal_y * uR(2)   +  normal_x * vR(2)

      else
      ! ---------- For classic cells ---------- !
      	uR(1)  =  dof%u( iR )
      	vR(1)  =  dof%v( iR )

      ! project in local (edge) coordinate system --> NOT USED in the code
      	uR(2)  =  mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
      	vR(2)  =  mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)

      end if ! end if  ( mesh%edge(ie)%boundary )


      !-------------------------------------------------------------------------------------------------------------------!
      ! Estimates of PHI (pressure + gravity potential)
      !-------------------------------------------------------------------------------------------------------------------!
      phiR = g * (hR + zR)
      phiL = g * (hL + zL)

   !================================================================================================================!
   !  Calculate Specific variables (necessary either to rusanov or entropy dissipatrive fluxes)
   !================================================================================================================!

   !-------------------------------------------------------------------------------------------------------------------!
   !   Estimates  well balanced variables - for rusanov like scheme
   !-------------------------------------------------------------------------------------------------------------------!
      ! ---------- Water height  ---------- !
      hL_rusanov  =  max( 0._rp , hL + zL - max( zL , zR ) )
      hR_rusanov  =  max( 0._rp , hR + zR - max( zL , zR ) )

      ! ---------- Potential (phi_rusanov) ---------- !
      ! left state
      if(hR + zR < zL  .AND. hL < hR ) then
          phiL_rusanov = phiR
      else
          phiL_rusanov = phiL
      end if

      ! Right state
      if(hL + zL < zR  .AND. hR < hL ) then
          phiR_rusanov = phiL
      else
          phiR_rusanov = phiR
      end if

        ! write(30,*) phiL, phiR, phiL_rusanov, phiR_rusanov

   !-------------------------------------------------------------------------------------------------------------------!
   !    Wave speed calculation for rusanov like scheme
   !-------------------------------------------------------------------------------------------------------------------!
      c1 =  sqrt(uL(1)**2 + vL(1)**2) + sqrt(g * hL_rusanov)
      c2 =  sqrt(uR(1)**2 + vR(1)**2) + sqrt(g * hR_rusanov)

      lambda_rusanov = max( c1, c2)
   !-------------------------------------------------------------------------------------------------------------------!
   !    Weight given to each scheme calculated for edge ie
   !-------------------------------------------------------------------------------------------------------------------!
      weight_rusanov = max(all_weight_rusanov(iL), all_weight_rusanov(iR) )

   !=============================================================================================================!
   !  Calculate FLUXES
   !=============================================================================================================!

   !-------------------------------------------------------------------------------------------------------------------!
   !  RUSANOV LIKE SCHEME
   !-------------------------------------------------------------------------------------------------------------------!


      F_rusanov = demi * ( &
                            hL_rusanov * (uL(1) * normal_x + vL(1) * normal_y )   &
                          + hR_rusanov * (uR(1) * normal_x + vR(1) * normal_y ) ) &
                 - lambda_rusanov * (hR_rusanov - hL_rusanov)

      G_rusanov(1) = demi * (&
                        hL_rusanov * (uL(1) * uL(1) * normal_x + uL(1) * vL(1) * normal_y )   &
      			      + hR_rusanov * (uR(1) * uR(1) * normal_x + uR(1) * vR(1) * normal_y ) ) &
                   - lambda_rusanov * (hR_rusanov * uR(1) - hL_rusanov * uL(1))

      G_rusanov(2) = demi * (&
                        hL_rusanov * (uL(1) * vL(1) * normal_x + vL(1) * vL(1) * normal_y ) &
      			      + hR_rusanov * (uR(1) * vR(1) * normal_x + vR(1) * vR(1) * normal_y )  ) &
                   - lambda_rusanov * (hR_rusanov * vR(1) - hL_rusanov * vL(1))

      phi_rusanov = demi * (phiR_rusanov + phiL_rusanov )

   !-------------------------------------------------------------------------------------------------------------------!
   !  ENTROPY DISSPATIVE SCHEME
   !-------------------------------------------------------------------------------------------------------------------!

  ! if the edge is at a boundary, the entropy disspative flux doesnt need to be calculated (because 1-weight_rusanov = 0)
    if (mesh%edge(ie)%boundary) then
      F_lf = 0._rp
      G_lf(1) = 0._rp
      G_lf(2) = 0._rp
      phi_lf = 0._rp
    else

      F_lf = demi * ( &
                    hL * (uL(1) * normal_x + vL(1) * normal_y ) &
      		     +  hR * (uR(1) * normal_x + vR(1) * normal_y )  ) &
            - 0.125_rp * gamma_lf * (periL/surfL + periR/surfR ) * (phiR - phiL)


        ! ---------- intermediate calculus -------------
      Fmax_lf = max(0._rp, F_lf)
      Fmin_lf = min(0._rp, F_lf)
        ! ---------- back to fluxes ------------

      G_lf(1) = uL(1) * Fmax_lf + uR(1) * Fmin_lf
      G_lf(2) = vL(1) * Fmax_lf + vR(1) * Fmin_lf

      phi_lf =   demi * (phiR + phiL )  &
                   - 0.25_rp * alpha_lf * g * (periL/surfL + periR/surfR ) * &
                   ( (uL(1) * normal_x + vL(1) * normal_y) - (uR(1) * normal_x + vR(1) * normal_y) )

    end if

!write(40,*) F_rusanov,  G_rusanov(1), G_rusanov(2), phi_rusanov
!write(50,*) F_lf, G_lf(1), G_lf(2), phi_lf



 ! ---- Combine fluxes
      Ftot_edge(1) = ( 1 - weight_rusanov) * F_lf          + weight_rusanov * F_rusanov
      
      Ftot_edge(2) = ( 1 - weight_rusanov) * G_lf(1)       + weight_rusanov * G_rusanov(1)   &
      +  normal_x * (( 1 - weight_rusanov) * phi_lf       + weight_rusanov * phi_rusanov) !hL *
      
      Ftot_edge(3) = ( 1 - weight_rusanov) * G_lf(2)       + weight_rusanov * G_rusanov(2) &
      +  normal_y * (( 1 - weight_rusanov) * phi_lf  + weight_rusanov * phi_rusanov) ! hL*


         !=============================================================================================================!
         !  Boundary post treatment :
         !    - Feedback control of bathy_cell in ghost cells to properly control the Qin imposed
         !    - Calculation of nflux sum for each inflow
         !=============================================================================================================!
! SOMETHING TO DO ?
!          if ( mesh%edge(ie)%boundary ) call boundary_post( Ftot_edge(1) , iR , mesh )

         !=============================================================================================================!
         !  Flux rotation and summation (as antisymmetric part to save time computation)
         !=============================================================================================================!

! write(20,*) ie, iL, iR, mesh%edge(ie)%boundary, mesh%edge(ie)%subdomain, Ftot_edge(:)

      Ftot_edge(1:3)  =  Ftot_edge(1:3)  *  edge_length


      tflux( 1 , iL )  =  tflux( 1 , iL )  +  Ftot_edge(1)
      tflux( 2 , iL )  =  tflux( 2 , iL )  +  Ftot_edge(2)
      tflux( 3 , iL )  =  tflux( 3 , iL )  +  Ftot_edge(3)

        if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then

      tflux( 1 , iR )  =  tflux( 1 , iR )  -  Ftot_edge(1)
      tflux( 2 , iR )  =  tflux( 2 , iR )  -  Ftot_edge(2)
      tflux( 3 , iR )  =  tflux( 3 , iR )  -  Ftot_edge(3)

        end if     	! end if  ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain )
      end if        ! end if  ( hL(1) > heps .or. hR(1) > heps )
   end do


!    write(*,*) "tflux( 1 , : )", tflux( 1 , : )
!    write(*,*) "tflux( 2 , : )", tflux( 2 , : )
!    write(*,*) "tflux( 3 , : )", tflux( 3 , : )

   open(10,file="flux.txt",status='replace',form='formatted')
   do i = 1, mesh%nc

   if(mesh%cell(i)%boundary) then
      tmp_index_cells_1 = mesh%cell(i)%cell
    do j = 1,4
        if(tmp_index_cells_1(j) > mesh%nc) then
                tmp_index_cells_2(j) =   mesh%cellb(tmp_index_cells_1(j) - mesh%nc)%typlim
        else
                tmp_index_cells_2(j) = "no_bc"
        end if
    end do
!           write(10,*) tmp_index_cells_2(:), tflux( 1 , i ), tflux( 2 , i ), tflux( 3 , i )
    else
!           write(10,*) "no_bc ","no_bc ","no_bc ","no_bc ", tflux( 1 , i ), tflux( 2 , i ), tflux( 3 , i )
    end if
   end do
   close(10)
   close(20)
   close(30)
   close(40)
   close(50)
!   read(*,*)
   !===================================================================================================================!
   !  Euler Time Step
   !===================================================================================================================!

   do i = 1,mesh%nc

      h  =  dof%h(i)
      u  =  dof%u(i)
      v  =  dof%v(i)

      dof%h(i)  =  max( 0._rp , h  -  dt * tflux(1,i) * mesh%cell(i)%invsurf )
      
      !================================================================================================================!
      !   Positivity cut-off
      !================================================================================================================!

      if ( dof%h(i) <= heps ) then

         dof%u(i)  =  0._rp
         dof%v(i)  =  0._rp

      else

         dof%u(i)  =  (  h * u  -  dt * ( tflux(2,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
         dof%v(i)  =  (  h * v  -  dt * ( tflux(3,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)

         !=============================================================================================================!
         !   Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
         !=============================================================================================================!

         if      ( friction == 1 ) then

            vel  =  sqrt( dof%u(i)**2 + dof%v(i)**2 )

            sfl  =  dof%h(i)**d2p3 + sqrt( dof%h(i)**d4p3 + 4._rp * dt * g * &
                   (manning(land(i))*dof%h(i)**manning_beta(land(i)))**2 * vel )

            sfl  =  2._rp * dof%h(i)**d2p3 / sfl

         else if ( friction == 2 ) then

            sfl  =  one - dt * manning( land(i) )

         else

            sfl  =  1._rp

         end if

         dof%u(i)  =  dof%u(i) * sfl
         dof%v(i)  =  dof%v(i) * sfl

      end if

   end do

   !===================================================================================================================!
   !  Calling MPI and filling ghost cells
   !===================================================================================================================!

   call com_dof( dof , mesh )

   call com_var_r( bathy_cell , mesh )                   ! Required MPI Communication due to inverse variable dependency

END SUBROUTINE low_froude_time_step







!~ !> \brief  Perform Low froude Time Step dedicated to Shallow-Water Equations
!~ !! \return dof updated after this new timestep
!~ SUBROUTINE low_froude_time_step( dof , mesh )

!~    USE m_common
!~    USE m_mesh
!~    USE m_mpi
!~    USE m_time_screen                                                                                              !NOADJ
!~    USE m_model

!~    implicit none

!~    !===================================================================================================================!
!~    !  Interface Variables
!~    !===================================================================================================================!

!~    type( msh ), intent(in   )  ::  mesh

!~    type( unk ), intent(inout)  ::  dof

!~    !===================================================================================================================!
!~    !  Local Variables
!~    !===================================================================================================================!
!~    real(rp), dimension( sw_nb , mesh%nc )  ::  tflux     !> Finite Volume total flux for each cell
   
!~    integer(ip)  ::  iL , iR                              !> Left and Right cells indexes to edge
   
   
!~    ! ======================================== LEVEL 1 -common-  Variables ======================================== !
!~    integer(ip)  ::  surfL , surfR                              !> Left and Right SURFACES of cells  (for low_froude scheme)
!~    integer(ip)  ::  periL , periR                              !> Left and Right PERIMETERS of cells  (for low_froude scheme)
!~    integer(ip)  :: edge_length                                 !> edge length  

!~    ! DOF
!~    real(rp)  ::  hL , uL(2) , vL(2)             !> Left  State in edge cell normal coordinates
!~    real(rp)  ::  hR , uR(2) , vR(2)             !> Right State in edge cell normal coordinates, zR correspond to bathymetry
!~    ! The first element of each variable correspond to its value in global coordinate system 
!~    !     second elemment correpond to its value in local coordinate systeme (relative to the edge, for cell k to cell ke)
!~    ! The second elements are used to define necessary input values for calc_boundary_state subroutine
   
!~    ! ZR ZL
!~    real(rp)  :: zR                          !> Bathymetry Left cell
!~    real(rp)  :: zL                          !> Bathymetry Right cell
   
!~    ! weight rusanob
!~    real(rp)  :: weight_rusanov             !> weight given to rusanov scheme, for one edge treatment (compare to entropy dissipative part)
!~    real(rp), dimension(:), allocatable  :: all_weight_rusanov    
   
!~    ! normal(x,y) is the normal vector of the edge ie, from cell k to cell ke
!~    real(rp)  :: normal_x    ! 
!~    real(rp)  :: normal_y    !
   
!~     ! g (acceleration pensenteur) is defined in a module

!~    real(rp)  :: phiL       ! pressure potential on cell k
!~    real(rp)  :: phiR       ! pressure potential on cell ke
   
!~   ! ======================================== LEVEL 2 -specific-  Variables ======================================== !

!~ ! >>> For rusanov scheme

!~    real(rp)  ::  hL_rusanov            !> reconstructed water height from andusse et al (Left)
!~    real(rp)  ::  hR_rusanov             !> reconstructed water height from andusse et al (Right)
   
!~    real(rp)  :: phiL_rusanov       !> reconstructed pressure potential on cell k
!~    real(rp)  :: phiR_rusanov       !> reconstructed pressure potential on cell ke   
   
!~    real(rp)  :: lambda_rusanov       !> wave speed for rusanov scheme
  
!~ ! >>> For entropy disspative scheme
!~    real(rp)  :: alpha_lf       !> weigth in entropy dissipative scheme    
!~    real(rp)  :: gamma_lf       !> weigth in entropy dissipative scheme     
   
!~   ! ======================================== LEVEL 3 -FLUXES-  Variables ======================================== !
!~   ! >>> RUSANOV SCHEME
!~    real(rp)  ::  F_rusanov      !> F part of rusanov flux
!~    real(rp)  ::  G_rusanov(2)   !> G part of rusanov flux, first element correpond to x-axis part, second to y axis part
!~    real(rp)  ::  phi_rusanov    !> Phi part of rusanov flux  
   
!~   ! >>> ENTROPY DISSIPATIVE SCHEME
!~    real(rp)  ::  F_lf      !> F part of entropy dissipative       flux
!~    real(rp)  ::  G_lf(2)   !> G part of entropy dissipative       flux, first element correpond to x-axis part, second to y axis part
!~    real(rp)  ::  phi_lf    !> Phi part of entropy dissipative       flux  
   
!~    ! >>> Temporary variable to store combination of both flux defined above
!~    real(rp)  ::  ftot_edge(3) !> Temporary variable to store combination of both flux defined above
!~ !
!~ !
!~ !


!~   ! ======================================== LEVEL 0 -temporary to store intermediate values -  Variables ======================================== !
!~       real(rp)  :: c1       !> temporary wave speed for rusanov scheme
!~       real(rp)  :: c2       !> temporary wave speed for rusanov scheme
!~       real(rp)  :: Fmax_lf, Fmin_lf !> temporary variable for G flux calculation in LF scheme
      
!~       real(rp)  :: h_cell    !>       tmp h value on one cell, used to interpolate weight rusanov coefficient
!~       real(rp)  :: tmp          !> temporary variable to make tests
   
!~     ! ======================================== UNUSED ?  ======================================== !
!~     real(rp)  ::  h , u , v                               !> Temporal primitive variables

!~    real(rp)  ::  vel                                     !> Velocity norm
   

!~    real(rp)  ::  sfl                                     !> Manning
 
!~    real(rp), dimension( sw_nb )  ::  nflux               !> Finite Volume normal edge flux  
!~    real(rp), dimension( sw_nb )  ::  lflux               !> Finite Volume edge flux in (x,y) coordinates
!~    !===================================================================================================================!
!~    !  Begin Subroutine
!~    !===================================================================================================================!
   
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  Initialise variables (that can be initialised)
!~    !-------------------------------------------------------------------------------------------------------------------!  
!~ write(*,*) "in low froude time step"

!~ tflux(:,:)  =  0._rp

!~ allocate(all_weight_rusanov(mesh%nc + mesh%ncb) )
!~ all_weight_rusanov(:) = 0._rp

!~ do i = 1,mesh%nc
!~     if (mesh%cell(i)%boundary) then ! In practice, \delta_{k} is set to 1 for cells with an edge at the physical boundary ..
!~         all_weight_rusanov(i) = 1
!~     else if ( dof%h(i) < heps ) then ! \delta_{k} is set to 1 for cells for which there is a neighbouring one with a dry state
!~         all_weight_rusanov(i) = 1
!~     else
!~     ! TODO: to interpolate :  then going from cell to cell with common edges from a given number of (typically from 3 to 5), the value goes to 0 linearly.
!~         all_weight_rusanov(i) = 0.5
!~     end if
!~ end do

!~  ! force alpha and gamma parameters
!~  alpha_lf = demi
!~  gamma_lf = demi
 
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  Calculate fluxes on each edge
!~    !-------------------------------------------------------------------------------------------------------------------!
!~   write(*,*) "start loop"
!~    do ie = 1,mesh%ne

!~       !================================================================================================================!
!~       !  Calculate Left and Right States
!~       !================================================================================================================!
      
!~       ! identify cell indices
!~       iL  =  mesh%edge(ie)%cell(1)
!~       iR  =  mesh%edge(ie)%cell(2)
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  get water height values
!~    !-------------------------------------------------------------------------------------------------------------------!
!~       ! Store h values for the two cells neibouring the edge number ie.
!~       hL  =  dof%h( iL )
!~       hR  =  dof%h( iR )
!~       ! store bathymetry values
!~       zL  =  bathy_cell( iL )
!~       zR  =  bathy_cell( iR )

!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  get geometrical properties
!~    !-------------------------------------------------------------------------------------------------------------------!
!~        !  STORE AREA AND PERIMETER OF CELLS FOR LOW_FROUDE SCHEME
!~       surfL = mesh%cell(iL)%surf
!~       surfR = mesh%cell(iR)%surf
!~       periL = mesh%cell(iL)%peri
!~       periR = mesh%cell(iR)%peri  
      
!~       ! get normal vector values for this edge
!~       normal_x = mesh%edge(ie)%normal%x 
!~       normal_y = mesh%edge(ie)%normal%y 
!~       ! get edge length
!~       edge_length = mesh%edge(ie)%length
      
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  If there is water (for one of the two edges), estimate fluxes on the edge.
!~    !  Else, flux on this edge is null  -> go to "Euler Time Step" section directly.
!~    !-------------------------------------------------------------------------------------------------------------------!
!~     if ( hL > heps .or. hR > heps ) then

!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  get flow speeds and their reprojection on the edge (necessary for boundary treatment)
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    ! >>> Left states
!~       uL(1)  =  dof%u( iL )
!~       vL(1)  =  dof%v( iL )

!~       uL(2)  =  normal_x * uL(1) + normal_y * vL(1)
!~       vL(2)  =  normal_x * vL(1) - normal_y * uL(1)


!~     ! >>> Right states

!~     if ( mesh%edge(ie)%boundary ) then
!~     ! ---------- For ghost cells ---------- !

!~         call calc_boundary_state( mesh , hL , zL , uL(2) , vL(2) , &
!~                                  hR , zR , uR(2) , vR(2) )

!~         ! reproject in global coordinate system
!~         uR(1) =  normal_x * uR(2)   -  normal_y * vR(2)
!~         vR(1) =  normal_y * uR(2)   +  normal_x * vR(2)


!~     else
!~     ! ---------- For classic cells ---------- !

!~         uR(1)  =  dof%u( iR )
!~         vR(1)  =  dof%v( iR )

!~         ! project in local (edge) coordinate system --> NOT USED in the code
!~         uR(2)  =  mesh%edge(ie)%normal%x * uR(1) + mesh%edge(ie)%normal%y * vR(1)
!~         vR(2)  =  mesh%edge(ie)%normal%x * vR(1) - mesh%edge(ie)%normal%y * uR(1)

!~     end if ! end if  ( mesh%edge(ie)%boundary )

!~    !-------------------------------------------------------------------------------------------------------------------!
!~    ! Estimates of PHI (pressure + gravity potential)
!~    !-------------------------------------------------------------------------------------------------------------------!
!~     ! theory
!~     phiR = g * (hR + zR)
!~     phiL = g * (hL + zL)

!~     ! en pratique pour être bien balancé !TODO TO DELETE and find why
!~     phiR =  (hR + zR)
!~     phiL =  (hL + zL)



!~    !================================================================================================================!
!~    !  Calculate Specific variables (necessary either to rusanov or entropy dissipatrive fluxes)
!~    !================================================================================================================!
!~ write(*,*) "estimate variables"

!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !   Estimates  well balanced variables - for rusanov like scheme
!~    !-------------------------------------------------------------------------------------------------------------------!
!~     ! ---------- Water height  ---------- !
!~     hL_rusanov  =  max( 0._rp , hL + zL - max( zL , zR ) )
!~     hR_rusanov  =  max( 0._rp , hR + zR - max( zL , zR ) )
          
!~     ! ---------- Potential (phi_rusanov) ---------- !
!~ !    if( mesh%edgeb( mesh%edge(ie)%lim )%typlim .EQ. 'wall' ) then
!~ !        ! if this ege is a boundary and is define as 'wall' there is no potential
!~ !        phiL_rusanov = 0
!~ !        phiR_rusanov = 0
!~ !    else
!~         ! else official treatment
!~         ! --- left state ---!
!~         if(hR + zR < zL  .AND. hL < hR ) then
!~             phiL_rusanov = phiR
!~         else
!~             phiL_rusanov = phiL
!~         end if
        
!~          ! --- left state ---!
!~         if(hL + zL < zR  .AND. hR < hL ) then
!~             phiR_rusanov = phiL
!~         else
!~             phiR_rusanov = phiR
!~         end if
!~ !    end if
!~  write(*,*) "anduse reconstruction ok"       
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !    Wave speed calculation for rusanov like scheme
!~    !-------------------------------------------------------------------------------------------------------------------!
!~            c1 =  sqrt(uL(1)**2 + vL(1)**2) + sqrt(g * hL_rusanov) 
!~            c2 =  sqrt(uR(1)**2 + vR(1)**2) + sqrt(g * hR_rusanov) 
!~            lambda_rusanov = max( c1,   c2)
           
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !    Weight given to each scheme calculated for edge ie
!~    !-------------------------------------------------------------------------------------------------------------------!
         
!~          weight_rusanov = max(all_weight_rusanov(iL), all_weight_rusanov(iR) )
!~          ! if bc, all_weight_rusanov(iL) = 1 and   all_weight_rusanov(iR) = 0
         
         
         
!~    !=============================================================================================================!
!~    !  Calculate FLUXES
!~    !=============================================================================================================!
!~ write(*,*) "calc flux"
!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  RUSANOV LIKE SCHEME
!~    !-------------------------------------------------------------------------------------------------------------------!
  
!~         F_rusanov = demi * ( &
!~                             hL_rusanov * (uL(1) * normal_x + vL(1) * normal_y ) &
!~                           + hR_rusanov * (uR(1) * normal_x + vR(1) * normal_y ) &
!~                             ) &
!~                    - lambda_rusanov * (hR_rusanov - hL_rusanov)
                   
!~         G_rusanov(1) = demi * (hL_rusanov * (uL(1) * uL(1) * normal_x + uL(1) * vL(1) * normal_y ) &
!~         		      + hR_rusanov * (uR(1) * uR(1) * normal_x + uR(1) * vR(1) * normal_y )  ) & 
!~                    - lambda_rusanov * (hR_rusanov * uR(1) - hL_rusanov * uL(1))  
                   
!~         G_rusanov(2) = demi * (hL_rusanov * (uL(1) * vL(1) * normal_x + vL(1) * vL(1) * normal_y ) &
!~         		      + hR_rusanov * (uR(1) * vR(1) * normal_x + vR(1) * vR(1) * normal_y )  ) & 
!~                    - lambda_rusanov * (hR_rusanov * vR(1) - hL_rusanov * vL(1))    
                   
!~         phi_rusanov = demi * (phiR_rusanov + phiL_rusanov )
!~       !                   write(*,*) "RUSANOV fluxes : ok"
      
!~ ! force flux to zer0 in case of wall edge     
      
!~        if (mesh%edge(ie)%boundary .AND. mesh%edgeb( mesh%edge(ie)%lim )%typlim .EQ. 'wall' ) then
!~        F_rusanov = 0
!~        G_rusanov(:) = 0
!~        phi_rusanov = 0
!~        end if
      

!~    !-------------------------------------------------------------------------------------------------------------------!
!~    !  ENTROPY DISSPATIVE SCHEME
!~    !-------------------------------------------------------------------------------------------------------------------!

!~   ! if the edge is at a boundary, the entropy disspative flux doesnt need to be calculated (because 1-weight_rusanov = 0)
!~   if (mesh%edge(ie)%boundary) then
!~   F_lf = 0._rp
!~   G_lf(1) = 0._rp
!~   G_lf(2) = 0._rp
!~   phi_lf = 0._rp
!~   else
!~ !         write(*,*) normal_x
!~ !         write(*,*) normal_y
!~ !         write(*,*) "tmp = ok", tmp
!~ !         write(*,*) "gamma_lf", gamma_lf
!~ !         write(*,*) "periL",  periL
!~ !         write(*,*) "periR",  periR  
!~ !         write(*,*) "surfL",  surfL
!~ !         write(*,*) "surfR",  surfR         
!~ !         write(*,*) "ratio surface", (periL/surfL + periR/surfR )
!~ !   
  
!~         F_lf = demi * (hL * (uL(1) * normal_x + vL(1) * normal_y ) &
!~                      + hR * (uR(1) * normal_x + vR(1) * normal_y )  ) & 
!~                    - 0.125_rp * gamma_lf * (periL/surfL + periR/surfR ) * (phiR - phiL)
                 
!~        !  write(*,*) "F ok"        
!~         ! tmp necessary calculus
!~         Fmax_lf = max(0._rp, F_lf)
!~         Fmin_lf = min(0._rp, F_lf)
!~        ! write(*,*) "Fmax/min ok"
!~         !! G flux
        
!~         G_lf(1) = uL(1) * Fmax_lf + uR(1) * Fmin_lf
!~         G_lf(2) = vL(1) * Fmax_lf + vR(1) * Fmin_lf        
!~       !  write(*,*) "G ok"
                         
!~         phi_lf =   demi * (phiR + phiL )  &
!~                    - 0.25_rp * alpha_lf * g * (periL/surfL + periR/surfR ) * &
!~                    ( (uL(1) * normal_x + vL(1) * normal_y) - (uR(1) * normal_x + vR(1) * normal_y) )
                   
!~   end if       
      
!~   !  write(*,*) "ENTROPY DISSIPATIVE fluxes : ok"
  

!~ ! if (mesh%edge(ie)%boundary ) then
!~ !   write(*,*) " type of boundary ", mesh%edgeb( mesh%edge(ie)%lim )%typlim 
!~ ! 
!~ !   write(*,*) "F lf ", F_lf
!~ !   write(*,*) "F rusanov ", F_rusanov  
!~ !   
!~ !   write(*,*) "G lf -x",  G_lf(1)
!~ !   write(*,*) "G rusanov -x",  G_rusanov(1)   
!~ ! 
!~ !   write(*,*) "G lf -y", G_lf(2)
!~ !   write(*,*) "G rusanov -y",  G_rusanov(2)   
!~ !   
!~ !   write(*,*) "PHI lf ", phi_lf
!~ !   write(*,*) "PHI rusanov",  phi_rusanov
!~ !   read(*,*)
!~ !   end if 

!~ write(*,*) "calc flux v2 ok"
!~  ! ---- Combine fluxes
!~  ! adapt weight rusanov to each edge (not scalar)
!~     Ftot_edge(1) = ( 1 - weight_rusanov) * F_lf + weight_rusanov * F_rusanov
    
!~     Ftot_edge(2) = ( 1 - weight_rusanov) * G_lf(1) + weight_rusanov * G_rusanov(1)   &
!~             + ( 1 - weight_rusanov) * phi_lf + weight_rusanov * phi_rusanov
            
!~     Ftot_edge(3) = ( 1 - weight_rusanov) * G_lf(2) + weight_rusanov * G_rusanov(2) &
!~             + ( 1 - weight_rusanov) * phi_lf + weight_rusanov * phi_rusanov    
            
!~ write(*,*) "combine flux"
            
!~       !     write(*,*) "Ftot: ok"     
!~          !=============================================================================================================!
!~          !  Boundary post treatment :
!~          !    - Feedback control of bathy_cell in ghost cells to properly control the Qin imposed
!~          !    - Calculation of nflux sum for each inflow
!~          !=============================================================================================================!
!~ ! SOMETHING TO DO ? 
!~ !          if ( mesh%edge(ie)%boundary ) call boundary_post( Ftot_edge(1) , iR , mesh )

!~          !=============================================================================================================!
!~          !  Flux rotation and summation (as antisymmetric part to save time computation)
!~          !=============================================================================================================!
         


!~          Ftot_edge(1:3)  =  Ftot_edge(1:3)  *  mesh%edge(ie)%length

         
!~          tflux( 1 , iL )  =  tflux( 1 , iL )  +  Ftot_edge(1)
!~          tflux( 2 , iL )  =  tflux( 2 , iL )  +  Ftot_edge(2)
!~          tflux( 3 , iL )  =  tflux( 3 , iL )  +  Ftot_edge(3)

         
!~ !          tflux( 2 , iL )  =  tflux( 2 , iL ) !HYDROSTATIC RECTONSTRUCTION .||   +  mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( ( hL(1)**2 - hL(2)**2 ) )
!~ ! 
!~ !          tflux( 3 , iL )  =  tflux( 3 , iL ) !HYDROSTATIC RECTONSTRUCTION||     +  mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( ( hL(1)**2 - hL(2)**2 ) )

!~          if ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) then

!~             tflux( 1 , iR )  =  tflux( 1 , iR )  -  Ftot_edge(1)
!~             tflux( 2 , iR )  =  tflux( 2 , iR )  -  Ftot_edge(2)
!~             tflux( 3 , iR )  =  tflux( 3 , iR )  -  Ftot_edge(3)

!~ !             tflux( 2 , iR )  =  tflux( 2 , iR ) ! HYDROSTATIC RECTONSTRUCTION|| -  mesh%edge(ie)%normal%x * mesh%edge(ie)%length * 0.5_rp * g * ( ( hR(1)**2 - hR(2)**2 ) )
!~ !             tflux( 3 , iR )  =  tflux( 3 , iR ) ! HYDROSTATIC RECTONSTRUCTION|| -  mesh%edge(ie)%normal%y * mesh%edge(ie)%length * 0.5_rp * g * ( & ( hR(1)**2 - hR(2)**2 ) )

!~          end if ! end if  ( .not. mesh%edge(ie)%boundary .and. .not. mesh%edge(ie)%subdomain ) 

!~       end if !     end if ( hL(1) > heps .or. hR(1) > heps )

!~    end do
!~        write(*,*) "end loop"    
!~   ! write(*,*) "END loop on edges "
!~    !===================================================================================================================!
!~    !  Euler Time Step
!~    !===================================================================================================================!

!~    do i = 1,mesh%nc

!~       h  =  dof%h(i)
!~       u  =  dof%u(i)
!~       v  =  dof%v(i)

!~       dof%h(i)  =  max( 0._rp , h  -  dt * tflux(1,i) * mesh%cell(i)%invsurf )

!~ !      dof%h(i)  =  h  -  dt * tflux(1,i) * mesh%cell(i)%invsurf

!~ !      if ( dof%h(i) < 0._rp ) call spread_mass_added( dof , mesh )

!~       !================================================================================================================!
!~       !   Positivity cut-off
!~       !================================================================================================================!

!~       if ( dof%h(i) <= heps ) then

!~          dof%u(i)  =  0._rp
!~          dof%v(i)  =  0._rp

!~       else

!~          dof%u(i)  =  (  h * u  -  dt * ( tflux(2,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)
!~          dof%v(i)  =  (  h * v  -  dt * ( tflux(3,i) * mesh%cell(i)%invsurf )  )  /  dof%h(i)

!~          !=============================================================================================================!
!~          !   Semi-Implicit Treatment of Friction Source Term (Manning/Strickler Formula)
!~          !=============================================================================================================!

!~          if      ( friction == 1 ) then

!~             vel  =  sqrt( dof%u(i)**2 + dof%v(i)**2 )

!~             sfl  =  dof%h(i)**d2p3 + sqrt( dof%h(i)**d4p3 + 4._rp * dt * g * &
!~                    (manning(land(i))*dof%h(i)**manning_beta(land(i)))**2 * vel )

!~             sfl  =  2._rp * dof%h(i)**d2p3 / sfl

!~          else if ( friction == 2 ) then

!~             sfl  =  one - dt * manning( land(i) )

!~          else

!~             sfl  =  1._rp

!~          end if

!~          dof%u(i)  =  dof%u(i) * sfl
!~          dof%v(i)  =  dof%v(i) * sfl

!~       end if

!~    end do

!~    !===================================================================================================================!
!~    !  Calling MPI and filling ghost cells
!~    !===================================================================================================================!

!~    call com_dof( dof , mesh )

!~    call com_var_r( bathy_cell , mesh )                   ! Required MPI Communication due to inverse variable dependencyEND
!~ END  SUBROUTINE low_froude_time_step
