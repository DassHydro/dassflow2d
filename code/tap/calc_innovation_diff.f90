!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of calc_innovation in forward (tangent) mode (with options fixinterface noISIZE):
!   variations   of useful results: *(*innovation.diff)
!   with respect to varying inputs: *bathy_cell *(*innovation.diff)
!                *(dof.h)
!   Plus diff mem management of: bathy_cell:in innovation:in *innovation.diff:in
!                dof.h:in
SUBROUTINE CALC_INNOVATION_DIFF(dof, dof_diff, mesh)
  USE M_COMMON ! Replaced by Perl Script
  USE M_MODEL ! Replaced by Perl Script
  USE M_OBS ! Replaced by Perl Script
  USE M_OBS_DIFF

  USE M_TAP_VARS ! Added by Perl Script -> Need to be filled !!!

  IMPLICIT NONE
  TYPE(UNK), INTENT(IN) :: dof
  TYPE(UNK), INTENT(IN) :: dof_diff
  TYPE(MSH), INTENT(IN) :: mesh
  INTEGER(ip) :: cell, searched_time, pt, n_average
  REAL(rp) :: h_mean, s_total
  REAL(rp) :: h_mean_diff
  INTRINSIC SIZE
  DO iobs=1,SIZE(station)
    searched_time = innovation(iobs)%ind_t
    IF (searched_time .LE. innovation(iobs)%nb_dt) THEN
      IF (tc .GE. station(iobs)%t(searched_time)) THEN
        h_mean = 0._rp
        s_total = 0._rp
        h_mean_diff = 0.0_8
        DO pt=1,SIZE(station(iobs)%pt)
          cell = station(iobs)%pt(pt)%cell
          IF (cell .GE. 0) THEN
            IF (dof%h(cell) .GT. 0) THEN
!test on water presence determining if cell is used for calculating h_average
              h_mean_diff = h_mean_diff + mesh%cell(cell)%surf*(dof_diff&
&               %h(cell)+bathy_cell_diff(cell))
              h_mean = h_mean + (dof%h(cell)+bathy_cell(cell))*mesh%cell&
&               (cell)%surf
              s_total = s_total + mesh%cell(cell)%surf
            END IF
          END IF
        END DO
        IF (s_total .GT. 0) THEN
          h_mean_diff = h_mean_diff/s_total
          h_mean = h_mean/s_total
        END IF
        innovation_diff(iobs)%diff(searched_time) = h_mean_diff
        innovation(iobs)%diff(searched_time) = h_mean - station(iobs)%h(&
&         searched_time)
        innovation(iobs)%ind_t = innovation(iobs)%ind_t + 1
      END IF
    END IF
  END DO
END SUBROUTINE CALC_INNOVATION_DIFF

