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
!> \file m_obs.f90
!! \brief This file includes m_obs module.
!! \details The file includes only m_obs module (see doc m_obs module).

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Module using Tapenade generated Output Files in /tap directory
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> Module m_obs.
!!
!! \details Module using Tapenade generated Output Files in /tap directory.
MODULE m_obs

   USE m_common
   USE m_mpi

   #if defined USE_SW_MONO || USE_SW_MULTI || USE_NS_MULTIFLUID

      USE m_model

   #else

      USE m_user_test

   #endif

   implicit none

   TYPE innovation_obs

      integer(ip) ::  nb_dt , nb_dx , ind_t

      real(rp), dimension(:), allocatable  ::  diff

   END TYPE

   type( innovation_obs ), dimension(:), allocatable  ::  innovation, innovW, innovQ

   integer(ip)  ::  nb_obs , nb_grp , iobs

   real(rp) :: w_mean


CONTAINS


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calculation of the cost function using model produced innovation vector
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief   Calculation of the cost function using model produced innovation vector
!! \details calculation of the cost function using inovation vector and regularization parameters
   SUBROUTINE calc_cost_function( cost , mesh )

      USE m_numeric

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( msh ), intent(in)  ::  mesh

      real(rp), intent(inout)  ::  cost

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      real(rp)  ::  cost_part(3) , filtered(4)

      integer(ip)  ::  idiff

      type( vec2d ), dimension( mesh%nc + mesh%ncb )  ::  grad_var

      !================================================================================================================!
      !  Loop on observations and quadratic norm of innovation vector / Regularization terms
      !================================================================================================================!

      #ifdef USE_SW_MONO

         if ( use_obs == 1 ) then

            !==========================================================================================================!
            !  Initialisation to zero of each part
            !==========================================================================================================!

            cost_part(:)  =  0._rp

	
            !==========================================================================================================!
            !  Loop on observations in Time/Space
            !==========================================================================================================!

            do iobs = 1,size( station )
               do idiff = 1,size( innovation( iobs )%diff )
                  cost_part(1)  =  cost_part(1)   +  station( iobs )%weight * innovation ( iobs )%diff( idiff )**2 
               end do

            end do

            !!!!!!!!!!!!!
            ! InnovW
            !!!!!!!!!!!!!            
            
            !do iobs = 1,size( station )
            !   do idiff = 1,size( innovW( iobs )%diff )
            !      cost_part(1)  =  cost_part(1)   +  station( iobs )%weight * innovW ( iobs )%diff( idiff )**2  
            !   end do
            !end do
            
            !!!!!!!!!!!!!            
            ! InnovQ with RMSE 
            !!!!!!!!!!!!!
            
            if (use_Qobs == 1 .or. use_Qobs_gr4 == 1) then
                do iobs = 1,size( stationQ )
!                 write(*,*) iobs, size( stationQ ), size( innovQ( iobs )%diff )
                    do idiff = 1,size( innovQ( iobs )%diff )
!                     write(*,*) iobs,  idiff, stationQ( iobs )%weight, innovQ ( iobs )%diff( idiff )
                        cost_part(1)  =  cost_part(1)   +   stationQ( iobs )%weight * innovQ ( iobs )%diff( idiff )**2

                    end do
                end do
            endif
            
            !!!!!!!!!!!!!            
            ! InnovQ with NSE 
            !!!!!!!!!!!!!
            
            
!             if (use_Qobs_gr4 == 1) then !Nash, temp use of cost part 2 and 3
!         
!                 do iobs = 1,size( stationQ )
! !                 write(*,*) iobs
!                     cost_part(3) = 0._rp
!                     cost_part(2) = 0._rp
!                     do idiff = 1,size( innovQ( iobs )%diff )
! !                     write(*,*) idiff
!                         cost_part(2)  =  cost_part(2)   +   (innovQ ( iobs                  )%diff( idiff ))**2
!                     end do
!                     do idiff = 1,size( innovQ( iobs )%diff )
! !                     write(*,*) idiff
!                         cost_part(3)  =  cost_part(3)   +   (innovQ ( iobs + size(stationQ) )%diff( idiff ))**2
!                     end do !FOR NSE
!                     cost_part(1)  =  cost_part(1)   +  cost_part(2) / cost_part(3)
!                     write(*,*)  'cost_part(1) in obs gr4 Nash',  cost_part(2) / cost_part(3) 
!                 end do
!             endif


             call mpi_sum_r( cost_part(1) )


!~             !==========================================================================================================!
!~             !  Bathymetry Regularization Term
!~             !==========================================================================================================!

             call FV_Cell_Grad( grad_var , bathy_cell , mesh )

            do i = 1,mesh%nc
               cost_part(2) = cost_part(2) + ( ( grad_var(i)%x )**2 +  ( grad_var(i)%y )**2 )
            end do

            call mpi_sum_r( cost_part(2) )

            cost_part(2) = cost_part(2) * regul_bathy

!~             !==========================================================================================================!
!~             !  Hydograph Regularization Term
!~             !==========================================================================================================!

             do k = 1,bc%nb_in

                filtered(1) = bc%hyd(k)%q(1)

                do i = 2,size( bc%hyd(k)%q(:) )-3

                   filtered(2) = filtered(1) + 0.2_rp * ( bc%hyd(k)%q(i  ) - filtered(1) )
                   filtered(3) = filtered(2) + 0.2_rp * ( bc%hyd(k)%q(i+1) - filtered(2) )
                   filtered(4) = filtered(3) + 0.2_rp * ( bc%hyd(k)%q(i+2) - filtered(3) )

                   cost_part(3) = cost_part(3) + ( bc%hyd(k)%q(i) - filtered(4) )**2

                   filtered(1) = filtered(4)

                end do

             end do

             cost_part(3) = cost_part(3) * regul_hydrograph

!~             !==========================================================================================================!
!~             !  Sum of each part
!~             !==========================================================================================================!

           cost = sum( cost_part )

         end if

      #endif

      !================================================================================================================!
      !  Fake operation for Tapenade Automatic Differentiation (Last operation ...)
      !================================================================================================================!

      cost = sqrt( cost**2 )

   END SUBROUTINE calc_cost_function


!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!  Calculation of the cost function directly updating its value during simulation (no observations)
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!

!> \brief   Calculation of the cost function directly updating its value during simulation (no observations)
!! \details sum of dof%h² times the weight given to the mesure point (station) ????
   SUBROUTINE update_cost_function( dof , cost )

      implicit none

      !================================================================================================================!
      !  Interface Variables
      !================================================================================================================!

      type( unk ), intent(in)  ::  dof

      real(rp), intent(inout)  ::  cost

      !================================================================================================================!
      !  Local Variables
      !================================================================================================================!

      integer(ip)  ::  cell , pt

      real(rp)  ::  h_mean

      !================================================================================================================!
      !  Begin
      !================================================================================================================!

      #ifdef USE_SW_MONO
      
        if (.not. allocated( station ) )  return

         do iobs = 1,size( station )

            if ( .not. test_dt_just_after( station( iobs )%dt ) ) cycle

            h_mean = 0._rp

            do pt = 1,size( station( iobs )%pt )

               cell = station( iobs )%pt( pt )%cell

               if ( cell < 0 ) cycle

               h_mean = h_mean + dof%h( cell )

            end do

            call mpi_sum_r( h_mean )

            h_mean = h_mean / real( size( station( iobs )%pt ) , 8 )

            cost = cost + station( iobs )%weight * h_mean**2

         end do

      #endif

   END SUBROUTINE update_cost_function


END MODULE m_obs
