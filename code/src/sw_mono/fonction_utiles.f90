
FUNCTION calculate_wetted_area(H_K, y_coords, b_coords) RESULT(area_K)
    !=======================================================================
    ! Computes the wetted area of a cross-section using the trapezoid rule.
    ! The size of the input arrays is determined at runtime using SIZE().
    ! The division by 2 is factored out and applied at the end.
    !
    ! INPUTS:
    !   H_K        : REAL(rp), INTENT(IN) :: Free surface elevation for the cell
    !   y_coords   : REAL(rp), DIMENSION(:), INTENT(IN) :: Array of the profile's y-positions
    !   b_coords   : REAL(rp), DIMENSION(:), INTENT(IN) :: Array of the bed elevations
    !
    ! OUTPUT:
    !   area_K     : REAL(rp) :: The computed wetted area
    !=======================================================================
    
    IMPLICIT NONE

    ! --- Argument Declarations ---
    REAL(rp), INTENT(IN) :: H_K
    REAL(rp), DIMENSION(:), INTENT(IN) :: y_coords, b_coords

    ! --- Result Declaration ---
    REAL(rp) :: area_K

    ! --- Local Variable Declarations ---
    INTEGER  :: n_points
    INTEGER  :: j
    REAL(rp) :: h_j, h_j1, delta_y

    ! --- Start of Algorithm ---

    ! Get the size of the input arrays
    n_points = SIZE(y_coords)

    ! Initialize the area
    area_K = 0.0_rp

    ! Safety check: at least 2 points are needed to form a trapezoid.
    IF (n_points < 2) THEN
        RETURN
    END IF

    ! DO loop over the n_points-1 segments of the profile.
    DO j = 1, n_points - 1
        ! Local water depths at the ends of segment j.
        h_j     = MAX(0.0_rp, H_K - b_coords(j))
        h_j1    = MAX(0.0_rp, H_K - b_coords(j+1))

        ! Width of segment j.
        delta_y = y_coords(j+1) - y_coords(j)

        ! Add the (width) x (sum of heights) term to the total sum.
        area_K  = area_K + delta_y * (h_j + h_j1)
    END DO

    ! Final division by 2
    area_K = area_K / 2.0_rp

END FUNCTION calculate_wetted_area
