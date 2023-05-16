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


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  VF4 method to compute the normal discrete gradient for a scalar
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


FUNCTION vf4_normal_discrete_gradient_scal( var , mesh )  RESULT( grad_var )

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   real(rp), dimension( mesh%ne )  ::  grad_var

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  iL , iR

   real(rp)  ::  dist

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   do ie = 1,mesh%ne

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)

      dist = .norm. mesh%edge(ie)%vcell

      grad_var(ie)  =  ( var( iR ) - var( iL ) ) / dist

   end do

END FUNCTION vf4_normal_discrete_gradient_scal


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  VF4 method to compute the normal discrete gradient for a vector (type vec2d)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


FUNCTION vf4_normal_discrete_gradient_vec2d( var , mesh )  RESULT( grad_var )

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   type( vec2d ), dimension( mesh%ne )  ::  grad_var

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  iL , iR

   real(rp)  ::  dist

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   do ie = 1,mesh%ne

      iL  =  mesh%edge(ie)%cell(1)
      iR  =  mesh%edge(ie)%cell(2)

      dist = .norm. mesh%edge(ie)%vcell

      grad_var(ie)  =  ( var( iR ) - var( iL ) ) / dist

   end do

END FUNCTION vf4_normal_discrete_gradient_vec2d


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Diamond method to compute the discrete gradient for a scalar
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


FUNCTION diamond_discrete_gradient_scal( var , mesh )  RESULT( grad_var )

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   type( vec2d ), dimension( mesh%ne )  ::  grad_var

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   real(rp), dimension( mesh%nn )  ::  var_node

   real(rp)  ::  dist(2) , correction , grad(2)

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   var_node  =  Least_Square_Cell_to_Node( var , mesh )

   do ie = 1,mesh%ne

      dist(1) = .norm. mesh%edge(ie)%vcell

      correction = ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) / dist(1)

      if ( correction + 10._rp * zerom >= one ) then

         correction = 1._rp

      else

         correction = acos( correction )

      end if

      dist(1) = dist(1) * correction

      dist(2) = mesh%edge(ie)%length * correction

      grad(1) = ( var     ( mesh%edge(ie)%cell(2) ) - var     ( mesh%edge(ie)%cell(1) ) ) / dist(1)
      grad(2) = ( var_node( mesh%edge(ie)%node(2) ) - var_node( mesh%edge(ie)%node(1) ) ) / dist(2)

      grad_var(ie)  =  grad(1) * mesh%edge(ie)%normal + &
                       grad(2) * mesh%edge(ie)%tangent

   end do

END FUNCTION diamond_discrete_gradient_scal


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Diamond method to compute the discrete gradient for a vector (type vec2d)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


FUNCTION diamond_discrete_gradient_vec2d( var , mesh )  RESULT( grad_var )

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   type( vec2d ), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  var

   type( tens2d ), dimension( mesh%ne )  ::  grad_var

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   type(vec2d), dimension( mesh%nn )  ::  var_node

   real(rp)  ::  dist(2) , correction

   type(vec2d)  ::  grad(2)

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   var_node%x  =  Least_Square_Cell_to_Node( var%x , mesh )
   var_node%y  =  Least_Square_Cell_to_Node( var%y , mesh )

   do ie = 1,mesh%ne

      dist(1) = .norm. mesh%edge(ie)%vcell

      correction = ( mesh%edge(ie)%vcell .dotprod. mesh%edge(ie)%normal ) / dist(1)

      if ( correction + 10._rp * zerom >= one ) then

         correction = 1._rp

      else

         correction = acos( correction )

      end if

      dist(1) = dist(1) * correction

      dist(2) = mesh%edge(ie)%length * correction

      grad(1) = ( var     ( mesh%edge(ie)%cell(2) ) - var     ( mesh%edge(ie)%cell(1) ) ) / dist(1)
      grad(2) = ( var_node( mesh%edge(ie)%node(2) ) - var_node( mesh%edge(ie)%node(1) ) ) / dist(2)

      grad_var(ie)  =  ( grad(1) .tensprod. mesh%edge(ie)%normal  ) + &
                       ( grad(2) .tensprod. mesh%edge(ie)%tangent )

   end do

END FUNCTION diamond_discrete_gradient_vec2d


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Interpolate cell unknows to node using a Least Square Method
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


FUNCTION Least_Square_Cell_to_Node( cell , mesh )  RESULT( node )

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   real(rp), dimension( mesh%nc + mesh%ncb ), intent(in)  ::  cell

   real(rp), dimension( mesh%nn )  ::  node

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   do k = 1,mesh%nn

      node(k) = sum( cell_to_node(k)%weights(:) * cell( mesh%node(k)%cell(:) ) )

   end do

END FUNCTION Least_Square_Cell_to_Node


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Interpolate cell unknows to node using a Least Square Method
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE Init_Cell_to_Node( mesh )

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  nb_cell

   real(rp)  ::  Rx , Ry , Ixx , Iyy , Ixy , D , Lx , Ly , xK , yK

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   allocate( cell_to_node( size( mesh%node ) ) )

   do k = 1,size(mesh%node)

      nb_cell = size( mesh%node(k)%cell(:) )

      allocate( cell_to_node(k)%weights( nb_cell ) )

      if ( mesh%node(k)%boundary .and. nb_cell == 3 ) then

         cell_to_node(k)%weights(1) = 1._rp
         cell_to_node(k)%weights(2) = 0._rp
         cell_to_node(k)%weights(3) = 0._rp

      else

         Rx   =  0._rp
         Ry   =  0._rp
         Ixx  =  0._rp
         Iyy  =  0._rp
         Ixy  =  0._rp

         do i = 1,nb_cell

            iK  =  mesh%node(k)%cell(i)

            if ( iK <= mesh%nc ) then

               xK  =  mesh%node(k)%coord%x - mesh%cell(iK)%grav%x
               yK  =  mesh%node(k)%coord%y - mesh%cell(iK)%grav%y

            else

               xK  =  mesh%node(k)%coord%x - mesh%cellb(iK-mesh%nc)%grav%x
               yK  =  mesh%node(k)%coord%y - mesh%cellb(iK-mesh%nc)%grav%y

            end if

            Rx   =  Rx   +  xK
            Ry   =  Ry   +  yK
            Ixx  =  Ixx  +  xK**2
            Iyy  =  Iyy  +  yK**2
            Ixy  =  Ixy  +  xK * yK

         end do

         if ( abs( Rx  ) < zerom )  Rx   =  0._rp
         if ( abs( Ry  ) < zerom )  Ry   =  0._rp
         if ( abs( Ixy ) < zerom )  Ixy  =  0._rp

         D  =  Ixx * Iyy - Ixy**2

!         if ( D < zerom ) STOP 'Warning Init_Cell_to_Node : D < zerom'

         Lx  =  ( Ixy * Ry - Iyy * Rx ) / D
         Ly  =  ( Ixy * Rx - Ixx * Ry ) / D

!         write(6,'(I6,I3,20ES15.7)') k , nb_cell , Ixx , Iyy, Ixy , Rx , Ry , Lx , Ly , D

         do i = 1,nb_cell

            iK  =  mesh%node(k)%cell(i)

            if ( iK <= mesh%nc ) then

               xK  =  mesh%node(k)%coord%x - mesh%cell(iK)%grav%x
               yK  =  mesh%node(k)%coord%y - mesh%cell(iK)%grav%y

            else

               xK  =  mesh%node(k)%coord%x - mesh%cellb(iK-mesh%nc)%grav%x
               yK  =  mesh%node(k)%coord%y - mesh%cellb(iK-mesh%nc)%grav%y

!               write(6,'(3I4,20ES15.7)') k , i , iK , xK , yK

            end if

            cell_to_node(k)%weights(i)  =  ( 1._rp + Lx * xK + Ly * yK ) / ( real(nb_cell,rp) + Lx * Rx + Ly * Ry )

         end do

!         write(6,'(I6,I3,20ES15.7)') k , nb_cell , cell_to_node(k)%weights(:)

      end if

   end do! ; STOP

END SUBROUTINE Init_Cell_to_Node
