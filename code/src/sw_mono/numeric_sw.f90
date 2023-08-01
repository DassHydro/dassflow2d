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
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE boundary_post( mass_flux , index_ghost , mesh )

   USE m_common
   USE m_mesh
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   real(rp), intent(in)  ::  mass_flux

   integer(ip), intent(in)  ::  index_ghost

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   integer(ip)  ::  group

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   ib = mesh%edge(ie)%lim

   group = mesh%edgeb(ib)%group
    
   if ( mesh%edgeb(ib)%typlim(1:8) == 'discharg' .or. &
        mesh%edgeb(ib)%typlim(1:3) == 'gr4' ) then !.or. &
!    mesh%edgeb(ib)%typlim(1:17) == 'internal_discharg') then !2nd item = TODO, 3rd item DEPRECATED

      bc%sum_mass_flux( group )  =  bc%sum_mass_flux( group )  -  mass_flux * mesh%edge(ie)%length

      if ( feedback_inflow == 1 ) then

         bathy_cell( index_ghost ) = bathy_cell( index_ghost ) + &

         coef_feedback * ( mass_flux - bc%inflow( mesh%neb + ib ) )

      end if

   end if

   if ( mesh%edgeb(ib)%typlim(1:6) == 'transm'   .or. &
        mesh%edgeb(ib)%typlim(1:7) == 'neumann'   .or. &
        mesh%edgeb(ib)%typlim(1:8) == 'ratcurve' .or. &
        mesh%edgeb(ib)%typlim(1:7) == 'zspresc'  .or. &
        mesh%edgeb(ib)%typlim(1:6) == 'hpresc' ) then
        !         mesh%edgeb(ib)%typlim(1:17)== 'internal_ratcurve'  .or. & !DEPRECATED

      bc%sum_mass_flux( group )  =  bc%sum_mass_flux( group )  +  mass_flux * mesh%edge(ie)%length

   end if
   
   if ( mesh%edgeb(ib)%typlim(1:11) == 'internal_2D') then
        
        bc%sum_mass_flux( group )  =  bc%sum_mass_flux( group )  +  mass_flux * mesh%edge(ie)%length
        
   endif

END SUBROUTINE boundary_post
                                                                                                                 !<NOADJ

!**********************************************************************************************************************!
!**********************************************************************************************************************!
!
!
!
!**********************************************************************************************************************!
!**********************************************************************************************************************!


SUBROUTINE spread_mass_added( dof , mesh )

   USE m_common
   USE m_mesh
   USE m_model

   implicit none

   !===================================================================================================================!
   !  Interface Variables
   !===================================================================================================================!

   type( msh ), intent(in)  ::  mesh

   type( unk ), intent(inout)  ::  dof

   !===================================================================================================================!
   !  Local Variables
   !===================================================================================================================!

   real(rp)  ::  coef_mass_balance

   integer(ip), dimension( maxed )  ::  icut

   !===================================================================================================================!
   !  Begin Subroutine
   !===================================================================================================================!

   k          =  count(                        mesh%cell(i)%cell(:) >= 1 .and. mesh%cell(i)%cell(:) <= mesh%nc )

   icut(1:k)  =  pack ( mesh%cell(i)%cell(:) , mesh%cell(i)%cell(:) >= 1 .and. mesh%cell(i)%cell(:) <= mesh%nc )

   coef_mass_balance  =  sum( dof%h( icut(1:k) ) * mesh%cell( icut(1:k) )%surf )

   if ( coef_mass_balance > - dof%h(i) * mesh%cell(i)%surf ) then

      coef_mass_balance  =  1._rp  +  dof%h(i) * mesh%cell(i)%surf / coef_mass_balance

      dof%h( icut(1:k) )  =  coef_mass_balance * dof%h( icut(1:k) )

   else

      mass_cut  =  mass_cut  -  dof%h(i) * mesh%cell(i)%surf  -  coef_mass_balance

      dof%h( icut(1:k) )  =  0._rp

   end if

   dof%h(i)  =  0._rp

END SUBROUTINE spread_mass_added                                                                                 !>NOADJ
