!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of newton_qin in forward (tangent) mode (with options fixinterface noISIZE):
!   variations   of useful results: zs
!   with respect to varying inputs: *bathy_cell *(dof.h) *(dof.u)
!                *(dof.v) qin
!   Plus diff mem management of: bathy_cell:in dof.h:in dof.u:in
!                dof.v:in
SUBROUTINE NEWTON_QIN_DIFF(qin, qin_diff, dof, dof_diff, mesh, zs, &
& zs_diff)
  USE M_MODEL ! Replaced by Perl Script

  USE M_TAP_VARS ! Added by Perl Script -> Need to be filled !!!

  IMPLICIT NONE
  TYPE(MSH), INTENT(IN) :: mesh
  TYPE(UNK), INTENT(IN) :: dof
  TYPE(UNK), INTENT(IN) :: dof_diff
  REAL(rp), INTENT(IN) :: qin
  REAL(rp), INTENT(IN) :: qin_diff
  REAL(rp), INTENT(OUT) :: zs
  REAL(rp), INTENT(OUT) :: zs_diff
  INTEGER(ip) :: av
  REAL(rp) :: r, c, z, dzs, s1, s2
  REAL(rp) :: r_diff, z_diff, s1_diff, s2_diff
  INTRINSIC SQRT
  INTRINSIC REAL
  INTRINSIC MAX
  INTRINSIC ABS
  REAL(rp) :: result1
  REAL(rp) :: result1_diff
  REAL(rp) :: temp
  result1 = SQRT(g)
  c = two*result1
  zs = zero
  av = 0
  zs_diff = 0.0_8
  DO ie=1,mesh%neb
    IF (mesh%edgeb(ie)%typlim(1:8) .EQ. 'discharg') THEN
      i = mesh%edge(mesh%edgeb(ie)%ind)%cell(1)
      IF (dof%h(i) .GT. heps) THEN
        av = av + 1
        zs_diff = zs_diff + bathy_cell_diff(i) + dof_diff%h(i)
        zs = zs + bathy_cell(i) + dof%h(i)
      END IF
    END IF
  END DO
  IF (av .GT. 0) THEN
    temp = REAL(av, rp)
    zs_diff = zs_diff/temp
    zs = zs/temp
  ELSE
    zs = zero
    zs_diff = 0.0_8
  END IF
  dzs = one
  k = 0
  DO WHILE (k .LE. 50 .AND. dzs .GT. 10._rp*zerom)
    dzs = zs
    s1 = zero
    s2 = zero
    s1_diff = 0.0_8
    s2_diff = 0.0_8
    DO ie=1,mesh%neb
      IF (mesh%edgeb(ie)%typlim(1:8) .EQ. 'discharg') THEN
        i = mesh%edge(mesh%edgeb(ie)%ind)%cell(1)
        IF (dof%h(i) .GT. heps) THEN
          temp = SQRT(dof%h(i))
          IF (dof%h(i) .EQ. 0.0) THEN
            result1_diff = 0.0_8
          ELSE
            result1_diff = dof_diff%h(i)/(2.0*temp)
          END IF
          result1 = temp
          r_diff = mesh%edge(mesh%edgeb(ie)%ind)%normal%x*dof_diff%u(i) &
&           + mesh%edge(mesh%edgeb(ie)%ind)%normal%y*dof_diff%v(i) + c*&
&           result1_diff
          r = dof%u(i)*mesh%edge(mesh%edgeb(ie)%ind)%normal%x + dof%v(i)&
&           *mesh%edge(mesh%edgeb(ie)%ind)%normal%y + c*result1
          IF (zero .LT. zs - bathy_cell(i)) THEN
            z_diff = zs_diff - bathy_cell_diff(i)
            z = zs - bathy_cell(i)
          ELSE
            z = zero
            z_diff = 0.0_8
          END IF
          temp = SQRT(z)
          IF (z .EQ. 0.0) THEN
            result1_diff = 0.0_8
          ELSE
            result1_diff = z_diff/(2.0*temp)
          END IF
          result1 = temp
          s1_diff = s1_diff - mesh%edge(mesh%edgeb(ie)%ind)%length*((r-c&
&           *result1)*z_diff+z*(r_diff-c*result1_diff))
          s1 = s1 - z*(r-c*result1)*mesh%edge(mesh%edgeb(ie)%ind)%length
          temp = SQRT(z)
          IF (z .EQ. 0.0) THEN
            result1_diff = 0.0_8
          ELSE
            result1_diff = z_diff/(2.0*temp)
          END IF
          result1 = temp
          s2_diff = s2_diff - mesh%edge(mesh%edgeb(ie)%ind)%length*(&
&           r_diff-d3p2*c*result1_diff)
          s2 = s2 - (r-d3p2*c*result1)*mesh%edge(mesh%edgeb(ie)%ind)%&
&           length
        END IF
      END IF
    END DO
    temp = (s1-qin)/s2
    zs_diff = zs_diff - (s1_diff-qin_diff-temp*s2_diff)/s2
    zs = zs - temp
    IF (1._rp - dzs/zs .GE. 0.) THEN
      dzs = 1._rp - dzs/zs
    ELSE
      dzs = -(1._rp-dzs/zs)
    END IF
    k = k + 1
  END DO
END SUBROUTINE NEWTON_QIN_DIFF

