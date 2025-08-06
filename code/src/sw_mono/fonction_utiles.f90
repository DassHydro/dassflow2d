FUNCTION calculate_porosity(H_k, y_coords, b_coords, h) RESULT(poro)
	IMPLICIT NONE
	REAL, INTENT(IN) :: H_k			! surface elevation in the cell
	REAL, INTENT(IN) :: y_coords(:)		! y coordinate of the cell
	REAL, INTENT(IN) :: b_coords(:)		! bathymetric porfile of the cell (upstream face)
	REAL, INTENT(IN) :: h			! water depth in the cell
	REAL :: poro
	REAL :: yi 
	REAL :: yj 

	yi = y_coords(1)
	yj = y_coords(SIZE(y_coords)) 
	poro = calculate_wetted_area(H_k, y_coords, b_coords) / (yj-yi) / h
