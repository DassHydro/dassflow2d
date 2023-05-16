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
!> \file geometry.f90
!! \brief This file contains subroutines to calculate basic mesh geometric properties, and get the id of the cell coresponding to given coordinates.
!! \details The file includes  m_common, m_mesh and m_mpi modules.

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calculate basic mesh geometric properties :
!
!  -  Cells surface, perimeter and gravity center
!  -  Edges length and normal
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief calculate basic mesh geometric properties
!! \details filling mesh variable of following informations
!! - Cells geometrical properties :
!!  - Surface
!!  - Gravity Center
!!  - Perimeter
!!  - Ghost Cells Symetric Gravity Center
!! -   Nodes geometrical properties :
!!  - Connected Cell
!!  - Connected Edges
!! -   Edges geometrical properties :
!!  - Length
!!  - Normal  Vector (oriented from 1 to 2)
!!  - Tangent Vector (oriented from 1 to 2)
!!  - Center Coordinates
!!  - linked Cells gravity center Vector

SUBROUTINE Mesh_Geometric_Properties( mesh )

   USE m_common
   USE m_mesh
   USE m_mpi
   USE m_model

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   type( msh ), intent(inout)  ::  mesh

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   type(vec2d)  ::  node_coord(4)

   integer(ip)  ::  nbed, connected_num_bc

   integer(ip), dimension(size(mesh%node))  ::  connected_to_node

!======================================================================================================================!
!
!  Cells geometrical properties :
!
!   - Surface
!   - Gravity Center
!
!  Nodes geometrical properties :
!
!   - Connected Cells
!
!======================================================================================================================!
WRITE(*,*) "INTO Mesh_Geometric_Properties"
   if ( maxed > 4 ) call Stopping_Program_Sub( 'Fix Mesh_Geometric_Properties for maxed > 4' )

   mesh%surf  =  0._rp

   connected_to_node(:)  =  0

   do i = 1,size(mesh%cell)

      nbed  =  mesh%cell(i)%nbed

      node_coord(1:nbed)  =  mesh%node( mesh%cell(i)%node(1:nbed) )%coord

      if      ( nbed == 3 ) then

         mesh%cell(i)%surf  =  0.5_rp * abs( ( node_coord(3)%x - node_coord(2)%x ) * node_coord(1)%y + &
                                             ( node_coord(1)%x - node_coord(3)%x ) * node_coord(2)%y + &
                                             ( node_coord(2)%x - node_coord(1)%x ) * node_coord(3)%y )

      else if ( nbed == 4 ) then

         mesh%cell(i)%surf  =  0.5_rp * abs( ( node_coord(4)%x - node_coord(2)%x ) * node_coord(1)%y + &
                                             ( node_coord(1)%x - node_coord(3)%x ) * node_coord(2)%y + &
                                             ( node_coord(2)%x - node_coord(4)%x ) * node_coord(3)%y + &
                                             ( node_coord(3)%x - node_coord(1)%x ) * node_coord(4)%y )

     else

         call Stopping_Program_Sub( 'nbed problem Mesh_Geometric_Properties' )

      end if

      mesh%surf  =  mesh%surf  +  mesh%cell(i)%surf

      mesh%cell(i)%invsurf  =  1._rp / mesh%cell(i)%surf

      mesh%cell(i)%grav%x  =  sum( node_coord(1:nbed)%x ) / real(nbed,rp)
      mesh%cell(i)%grav%y  =  sum( node_coord(1:nbed)%y ) / real(nbed,rp)

      connected_to_node( mesh%cell(i)%node(1:nbed) ) = &
      connected_to_node( mesh%cell(i)%node(1:nbed) ) + 1

!write(*,*) "     mesh%cell(i)%surf  =",      mesh%cell(i)%surf
!write(*,*) "     connected_to_node( mesh%cell(i)%node(1:nbed) )  =",    connected_to_node( : )


   end do

!WRITE(*,*) "INTO step2"
!WRITE(*,*) "mesh%neb=",mesh%neb

   do ie = 1,mesh%neb
      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) = &
      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) + 1
   end do


!WRITE(*,*) "Def nodes"

   do i = 1,size(mesh%node)

      allocate( mesh%node(i)%cell( connected_to_node(i) ) )

   end do

   connected_to_node(:)  =  0

!   WRITE(*,*) "Def ++++"
   do i = 1,size(mesh%cell)

      nbed  =  mesh%cell(i)%nbed

      connected_to_node( mesh%cell(i)%node(1:nbed) ) = &
      connected_to_node( mesh%cell(i)%node(1:nbed) ) + 1

      do j = 1,nbed

         mesh%node( mesh%cell(i)%node(j) )%cell( connected_to_node( mesh%cell(i)%node(j) ) ) = i

      end do

  enddo

do i = 1, mesh%ne

!WRITE(*,*) "check mesh%edge(i)%boundary"

      if ( mesh%edge(i)%boundary ) then

!WRITE(*,*) "i=", i
!WRITE(*,*) "mesh%edge(i)%lim", mesh%edge(i)%lim
!WRITE(*,*) "mesh%edge(i)%cell(:)", mesh%edge(i)%cell(:)
!WRITE(*,*) "mesh%edgeb(mesh%edge(i)%lim)%typlim", mesh%edgeb(mesh%edge(i)%lim)%typlim

        if ( mesh%edgeb(mesh%edge(i)%lim)%typlim == 'internal_2D' ) then !then change connectivity to connected 1D-like cell
                read(bc%typ (mesh%edgeb(mesh%edge(i)%lim)%group, 3 ),'(i3)') connected_num_bc !Get connectivity to 1D-like cell from bc.txt

                do ib = 1,mesh%neb
                    if ( mesh%edgeb(ib)%group  ==  connected_num_bc ) then

                        mesh%edge(i)%cell1D2D  =  mesh%edge( mesh%edgeb( ib )%ind )%cell(1) !Get id of the single 1D-like cell with interface in the connected bc number

                    endif
                enddo


        endif

!  WRITE(*,*) "-----------------------------"
!  WRITE(*,*) "i=", i
!  WRITE(*,*) "mesh%edge(i)%lim", mesh%edge(i)%lim
!  WRITE(*,*) "mesh%edge(i)%cell(:)", mesh%edge(i)%cell(:)
!  WRITE(*,*) "mesh%edgeb(mesh%edge(i)%lim)%typlim", mesh%edgeb(mesh%edge(i)%lim)%typlim
!  WRITE(*,*) "-----------------------------"
      endif

   end do

!WRITE(*,*) " ==> TREAT NODES"
   do ie = 1,mesh%neb

      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) = &
      connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(:) ) + 1

      if ( mesh%edge( mesh%edgeb(ie)%ind )%cell(2) > mesh%nc ) then

         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)
         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%cell(2)

      else

         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(1) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%lim + mesh%nc
         mesh%node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) )%cell( connected_to_node( mesh%edge( mesh%edgeb(ie)%ind )%node(2) ) ) = mesh%edge( mesh%edgeb(ie)%ind )%lim + mesh%nc

      end if

   end do




!WRITE(*,*) "NODES OK"
!======================================================================================================================!
!
!  Edges geometrical properties :
!
!   - Length
!   - Normal  Vector (oriented from 1 to 2)
!   - Tangent Vector (oriented from 1 to 2)
!   - Center Coordinates
!   - linked Cells gravity center Vector
!
!  Cells geometrical properties :
!
!   - Perimeter
!   - Ghost Cells Symetric Gravity Center
!
!  Nodes geometrical properties :
!
!   - Connected Edges
!
!======================================================================================================================!

   mesh%cell(:)%peri  =  0._rp

   connected_to_node(:)  =  0

!WRITE(*,*) "edge1-----"
   do ie = 1,size(mesh%edge)

      node_coord(1:2)  =  mesh%node( mesh%edge(ie)%node(1:2) )%coord
!WRITE(*,*) "  node_coord(1:2)  ",   node_coord(1:2)
      mesh%edge(ie)%center%x  =  demi * sum( node_coord(1:2)%x )
      mesh%edge(ie)%center%y  =  demi * sum( node_coord(1:2)%y )

      mesh%edge(ie)%length  =  sqrt( ( node_coord(2)%x - node_coord(1)%x )**2 + &
                                     ( node_coord(2)%y - node_coord(1)%y )**2 )
!WRITE(*,*) "    mesh%edge(ie)%length   ",     mesh%edge(ie)%length

      mesh%edge(ie)%normal%x   =  ( node_coord(1)%y - node_coord(2)%y ) / mesh%edge(ie)%length
      mesh%edge(ie)%normal%y   =  ( node_coord(2)%x - node_coord(1)%x ) / mesh%edge(ie)%length

      mesh%edge(ie)%tangent%x  =  ( node_coord(2)%x - node_coord(1)%x ) / mesh%edge(ie)%length
      mesh%edge(ie)%tangent%y  =  ( node_coord(2)%y - node_coord(1)%y ) / mesh%edge(ie)%length

      if ( .not. mesh%edge(ie)%boundary ) then
!WRITE(*,*) "   not boundary   "
         mesh%cell( mesh%edge(ie)%cell(1:2) )%peri  =  &
         mesh%cell( mesh%edge(ie)%cell(1:2) )%peri  +  mesh%edge(ie)%length
!WRITE(*,*) "    mesh%cell( mesh%edge(ie)%cell(1:2) )%peri   ",  mesh%cell( mesh%edge(ie)%cell(1:2) )%peri

         mesh%edge(ie)%vcell%x = mesh%cell( mesh%edge(ie)%cell(2) )%grav%x - &
                                 mesh%cell( mesh%edge(ie)%cell(1) )%grav%x

         mesh%edge(ie)%vcell%y = mesh%cell( mesh%edge(ie)%cell(2) )%grav%y - &
                                 mesh%cell( mesh%edge(ie)%cell(1) )%grav%y
!WRITE(*,*) "  mesh%edge(ie)%vcell%y   ",  mesh%edge(ie)%vcell%y

         if ( ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) < 0._rp ) then

            mesh%edge(ie)%normal%x  = - mesh%edge(ie)%normal%x
            mesh%edge(ie)%normal%y  = - mesh%edge(ie)%normal%y

         end if

         mesh%edge(ie)%v_edge_cell(1)%x = mesh%edge(ie)%center%x - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%x
!WRITE(*,*) "  mesh%edge(ie)%v_edge_cell(1)%x  ",mesh%edge(ie)%v_edge_cell(1)%x
         mesh%edge(ie)%v_edge_cell(1)%y = mesh%edge(ie)%center%y - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%y

         mesh%edge(ie)%v_edge_cell(2)%x = mesh%edge(ie)%center%x - &
                                          mesh%cell( mesh%edge(ie)%cell(2) )%grav%x

         mesh%edge(ie)%v_edge_cell(2)%y = mesh%edge(ie)%center%y - &
                                          mesh%cell( mesh%edge(ie)%cell(2) )%grav%y
!WRITE(*,*) "   end not boundary   "
      else

!WRITE(*,*) "  In boundary edge  "
         mesh%cell( mesh%edge(ie)%cell(1) )%boundary = .true.

         mesh%cell( mesh%edge(ie)%cell(1) )%peri  =  &
         mesh%cell( mesh%edge(ie)%cell(1) )%peri  +  mesh%edge(ie)%length

         mesh%edge(ie)%vcell%x = 2._rp * ( mesh%edge(ie)%center%x - mesh%cell( mesh%edge(ie)%cell(1) )%grav%x )
         mesh%edge(ie)%vcell%y = 2._rp * ( mesh%edge(ie)%center%y - mesh%cell( mesh%edge(ie)%cell(1) )%grav%y )

         if ( ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) < 0._rp ) then

            mesh%edge(ie)%normal%x  = - mesh%edge(ie)%normal%x
            mesh%edge(ie)%normal%y  = - mesh%edge(ie)%normal%y

         end if

         mesh%edge(ie)%v_edge_cell(1)%x = mesh%edge(ie)%center%x - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%x

         mesh%edge(ie)%v_edge_cell(1)%y = mesh%edge(ie)%center%y - &
                                          mesh%cell( mesh%edge(ie)%cell(1) )%grav%y

         ib = mesh%edge(ie)%lim

         mesh%cellb( ib )%grav%x  =  mesh%cell( mesh%edge(ie)%cell(1) )%grav%x + mesh%edge(ie)%vcell%x
         mesh%cellb( ib )%grav%y  =  mesh%cell( mesh%edge(ie)%cell(1) )%grav%y + mesh%edge(ie)%vcell%y

!WRITE(*,*) "  Out boundary edge  "
      end if

!WRITE(*,*) "  tmp ",  connected_to_node( :)
      connected_to_node( mesh%edge(ie)%node(1:2) ) = &
      connected_to_node( mesh%edge(ie)%node(1:2) ) + 1

!WRITE(*,*) "   connected_to_node( :)  ",  connected_to_node( :)
   end do

!WRITE(*,*) "edge2-----"
   do ib = 1,size(mesh%edgeb)
!WRITE(*,*) "mesh%edgeb(ib)%typlim "

      if ( mesh%edgeb(ib)%typlim == 'periodic' ) then

         mesh%edge( mesh%edgeb(ib)%ind )%v_edge_cell(2)%x = mesh%edge( mesh%edgeb(ib)%perio )%center%x - &
                                                            mesh%cell( mesh%cellb(ib)%cell )%grav%x

         mesh%edge( mesh%edgeb(ib)%ind )%v_edge_cell(2)%y = mesh%edge( mesh%edgeb(ib)%perio )%center%y - &
                                                            mesh%cell( mesh%cellb(ib)%cell )%grav%y

      end if

   end do


!WRITE(*,*) "edge3-----"

   do i = 1,size(mesh%node)
!~ 		if (proc == 2) print *, proc, i, connected_to_node(i)
      allocate( mesh%node(i)%edge( connected_to_node(i) ) )

   end do

   connected_to_node(:)  =  0

!~ 	print *, proc, " : edge info"
!~ 	 seg fault for proc 2 in this loop
   do ie = 1,size(mesh%edge)
!~ 		if (proc==2) print *, "ite : " ,ie ,"/",size(mesh%edge), " les noeuds appeles : ", mesh%edge(ie)%node(1:2), "/", size(connected_to_node), size(mesh%node)
      connected_to_node( mesh%edge(ie)%node(1:2) ) = &
      connected_to_node( mesh%edge(ie)%node(1:2) ) + 1
!~ 		print *, "ite : " , ie, "/", size(mesh%edge), connected_to_node( mesh%edge(ie)%node(1:2) ), "tailles des tableaux edge : ", &
!~ 		size(mesh%node( mesh%edge(ie)%node(1) )%edge),size(mesh%node( mesh%edge(ie)%node(2) )%edge)

!~ 		if (ie==size(mesh%edge)) print *, mesh%edge(ie)%node(1:2)
      mesh%node( mesh%edge(ie)%node(1) )%edge( connected_to_node( mesh%edge(ie)%node(1) ) ) = ie
      mesh%node( mesh%edge(ie)%node(2) )%edge( connected_to_node( mesh%edge(ie)%node(2) ) ) = ie
!~ 	print *, ie, "/", size(mesh%edge)
   end do

!~    print *, proc, "finishing edges info"

WRITE(*,*) "OUT Mesh_Geometric_Properties"
END SUBROUTINE Mesh_Geometric_Properties


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  MUSCL Interpolation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief calculate  MUSCL Interpolation
!! \details Calculate  MUSCL Interpolation - NEVER CALLED
real(kind=rp) FUNCTION muscl_interp( var , grad_var , point , i_cell , mesh )

   USE m_mesh

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE(msh), intent(in)  ::  mesh

   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   TYPE(vec2d), dimension( mesh%nc ), intent(in)  ::  grad_var

   TYPE(vec2d), intent(in)  ::  point

   integer(ip), intent(in)  ::  i_cell

!======================================================================================================================!
!  Begin Function
!======================================================================================================================!

   muscl_interp  =  var(i_cell) + ( grad_var(i_cell)%x * ( point%x - mesh%cell(i_cell)%grav%x ) + &
                                    grad_var(i_cell)%y * ( point%y - mesh%cell(i_cell)%grav%y ) )

END FUNCTION muscl_interp


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Radial Basis Functions Interpolation
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief - NEVER CALLED - calculate Shepard interpolation
!! \details Radial Basis Functions Interpolation
real(kind=rp) FUNCTION Shepard_interp( var , point , i_cell , mesh )

   USE m_mesh

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE( msh ), intent(in)  ::  mesh

   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   TYPE( vec2d ), intent(in)  ::  point

   integer(ip), intent(in)  ::  i_cell

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   real(rp), parameter  ::  shepard_coef  =  1.5_rp       ! Power Law Basis Function ( > 0.5 and <= 1.5 )

   real(rp)  ::  basis_func_val , sum_d , sum_n

   integer(ip)  ::  neib , i_celle , i_celln

!======================================================================================================================!
!  Begin Function
!======================================================================================================================!

   basis_func_val = one / ( zerom + ( point%x - mesh%cell(i_cell)%grav%x )**2 + &
                                    ( point%y - mesh%cell(i_cell)%grav%y )**2 )**shepard_coef

   sum_n = var(i_cell) * basis_func_val
   sum_d =               basis_func_val

   do i_celle = 1,maxed

      i_celln = mesh%cell(i_cell)%cell(i_celle)

      basis_func_val = one / ( zerom + ( point%x - mesh%cell(i_celln)%grav%x )**2 + &
                                       ( point%y - mesh%cell(i_celln)%grav%y )**2 )**shepard_coef

      sum_n = sum_n + var(i_celln) * basis_func_val
      sum_d = sum_d +                basis_func_val

   end do

   Shepard_interp = sum_n / sum_d

END FUNCTION Shepard_interp


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Testing if a given point is in a Cell by Ray-Casting Technique (borders and corners included)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief Testing if a given point is in a Cell (borders and corners included)
!! \details Testing if a given point is in a Cell calling Ray-Casting Technique
logical FUNCTION point_in_cell( point , cell_ind , mesh )

   USE m_mesh

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   TYPE( vec2d ), intent(in)  ::  point

   integer(ip), intent(in)  ::  cell_ind

   TYPE( msh ), intent(in)  ::  mesh

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   integer(ip)  ::  cnt , n1 , n2

   logical  ::  ray_inter

   type( vec2d )  ::  point_offset

!======================================================================================================================!
!  Begin Function
!======================================================================================================================!

   cnt = 0

   point_offset = point * ( 1._rp - 10._rp * zerom )

   do ie = 1,maxed

      n1 = mesh%cell( cell_ind )%node(      ie               )
      n2 = mesh%cell( cell_ind )%node( mod( ie , maxed ) + 1 )

      if ( point_offset%x == mesh%node(n1)%coord%x .and. point_offset%y == mesh%node(n1)%coord%y ) then

         point_in_cell = .true. ; return

      else if ( ray_inter( point_offset , mesh%node(n1)%coord , mesh%node(n2)%coord ) ) then

         cnt = cnt + 1

      end if

   end do

   if ( mod( cnt , 2 ) == 0 ) then

      point_in_cell = .false.

   else

      point_in_cell = .true.

   end if

END FUNCTION point_in_cell

!> \brief APPLYING Ray-Casting Technique
!! \details APPLYING Ray-Casting Technique
logical FUNCTION ray_inter( p , a0 , b0 )

   USE m_mesh

   implicit none

!======================================================================================================================!
!  Interface Variables
!======================================================================================================================!

   type(vec2d), intent(in)  ::  p , a0 , b0

!======================================================================================================================!
!  Local Variables
!======================================================================================================================!

   type(vec2d)  ::  a , b

   real(rp)  ::  m_red , m_blue

!======================================================================================================================!
!  Begin Function
!======================================================================================================================!

   ray_inter = .false.

   if ( ( p%y /= a0%y ) .and. ( p%y /= b0%y ) ) then

      if ( a0%y > b0%y ) then
         b = a0
         a = b0
      else
         a = a0
         b = b0
      end if

      if ( ( p%y > b%y ) .or. ( p%y < a%y ) ) return

      if ( p%x > max( a%x , b%x ) ) return

      if ( p%x < min( a%x , b%x ) ) then

         ray_inter = .true.

      else

         if ( abs( a%x - b%x ) > tiny( a%x ) ) then
            m_red = ( b%y - a%y ) / ( b%x - a%x )
         else
            m_red = huge( m_red )
         end if

         if ( abs( a%x - p%x ) > tiny( a%x ) ) then
            m_blue = ( p%y - a%y ) / ( p%x - a%x )
         else
            m_blue = huge( m_blue )
         end if

         if ( m_blue >= m_red ) then
            ray_inter = .true.
         else
            ray_inter = .false.
         end if

      end if

   else if ( ( p%x /= a0%x ) .and. ( p%x /= b0%x ) ) then

      if ( a0%x > b0%x ) then
         b = a0
         a = b0
      else
         a = a0
         b = b0
      end if

      if ( ( p%x > b%x ) .or. ( p%x < a%x ) ) return

      if ( p%y > max( a%y , b%y ) ) return

      if ( p%y < min( a%y , b%y ) ) then

         ray_inter = .true.

      else

         if ( abs( a%y - b%y ) > tiny( a%y ) ) then
            m_red = ( b%x - a%x ) / ( b%y - a%y )
         else
            m_red = huge( m_red )
         end if

         if ( abs( a%y - p%y ) > tiny( a%y ) ) then
            m_blue = ( p%x - a%x ) / ( p%y - a%y )
         else
            m_blue = huge( m_blue )
         end if

         if ( m_blue >= m_red ) then
            ray_inter = .true.
         else
            ray_inter = .false.
         end if

      end if

   end if

END FUNCTION ray_inter



!> \brief Get the ID of the cell containing given point
!! \details return the ID of the cell containing given point
FUNCTION search_cell_inc_point( mesh , point ) RESULT( cell )

   USE m_mesh

   implicit none

   type( msh ), intent(in)  ::  mesh

   type( vec2d ), intent(in)  ::  point

   integer(ip)  ::  cell , icell

   logical  ::  point_in_cell

   cell = -1

   do icell = 1,mesh%nc

      if ( point_in_cell( point , icell , mesh ) ) then

         cell  =  icell

         exit

      end if

   end do

END FUNCTION search_cell_inc_point
