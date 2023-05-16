!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of sw_post_treatment in reverse (adjoint) mode (with options fixinterface noISIZE):
!   gradient     of useful results: *(bc.sum_mass_flux)
!   with respect to varying inputs: *(bc.sum_mass_flux)
!   Plus diff mem management of: bc.sum_mass_flux:in
SUBROUTINE SW_POST_TREATMENT_BACK(dof, mesh)
  USE M_COMMON ! Replaced by Perl Script
  USE M_MESH ! Replaced by Perl Script
  USE M_MODEL ! Replaced by Perl Script

  USE M_TAP_VARS ! Added by Perl Script -> Need to be filled !!!

  IMPLICIT NONE
  TYPE(UNK), INTENT(IN) :: dof
  TYPE(MSH), INTENT(IN) :: mesh
  INTEGER(ip) :: num_bc
  INTEGER :: branch
  DO num_bc=1,bc%nb
    IF (temp_scheme(1:2) .EQ. 'rk' .OR. temp_scheme(1:4) .EQ. 'imex') &
&   THEN
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END DO
  DO num_bc=bc%nb,1,-1
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) bc_back%sum_mass_flux(num_bc) = 0.5_rp*bc_back%&
&       sum_mass_flux(num_bc)
  END DO
END SUBROUTINE SW_POST_TREATMENT_BACK

