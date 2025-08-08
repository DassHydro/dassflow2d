FUNCTION calculate_wetted_area_parabolic(H_k, y1, y_min, yN, b_min) RESULT(area_k)
    !=======================================================================
    ! Calcule l'aire mouillée en utilisant une parabole DÉPENDANTE DU NIVEAU D'EAU,
    ! selon les équations fournies.
    ! ATTENTION : Cette méthode a des incohérences physiques.
    !
    ! ENTRÉES:
    !   H_k        : REAL(rp), INTENT(IN) :: Cote de la surface libre
    !   y1         : REAL(rp), INTENT(IN) :: Coordonnée y de la rive gauche
    !   y_min      : REAL(rp), INTENT(IN) :: Coordonnée y du point le plus bas
    !   yN         : REAL(rp), INTENT(IN) :: Coordonnée y de la rive droite
    !   b_min      : REAL(rp), INTENT(IN) :: Altitude du fond au point le plus bas
    !
    ! SORTIE:
    !   area_k     : REAL(rp) :: L'aire mouillée calculée
    !=======================================================================
    IMPLICIT NONE

    ! --- Arguments ---
    REAL(rp), INTENT(IN) :: H_k, y1, y_min, yN, b_min
    
    ! --- Résultat ---
    REAL(rp) :: area_k

    ! --- Variables Locales ---
    REAL(rp) :: a, b, c, den
    REAL(rp) :: y_start, y_end

    ! Si le niveau d'eau est sous le fond, l'aire est nulle.
    IF (H_k <= b_min) THEN
        area_k = 0.0_rp
        RETURN
    END IF

    ! ======================================================================
    ! PARTIE 1 : CALCUL DES COEFFICIENTS a, b, c SELON VOTRE FORMULE
    ! ======================================================================
    
    ! Dénominateur commun
    den = (y1 - y_min) * (yN - y_min)
    IF (ABS(den) < 1.0E-9_rp) THEN ! Évite la division par zéro
        area_k = 0.0_rp
        RETURN
    END IF
    
    ! Coefficient 'a'
    a = (b_min - H_k) / den
    
    ! Coefficient 'b' (en supposant la symétrie, comme dans votre formule)
    b = -a * (y1 + yN)
    
    ! Coefficient 'c'
    c = H_k - a * y1**2 - b * y1

    ! ======================================================================
    ! PARTIE 2 : CALCULER L'AIRE MOUILLÉE
    ! ======================================================================

    ! Dans ce modèle, par définition, la parabole coupe la surface de l'eau
    ! aux points y1 et yN. Ce sont donc les bornes de l'intégration.
    y_start = MIN(y1, yN)
    y_end   = MAX(y1, yN)

    ! Calculer l'intégrale exacte de h(y) = H_k - (ay^2+by+c) entre y_start et y_end
    area_k = (H_k - c) * (y_end - y_start) - &
             (b / 2.0_rp) * (y_end**2 - y_start**2) - &
             (a / 3.0_rp) * (y_end**3 - y_start**3)
    
    area_k = MAX(0.0_rp, area_k) ! Assurer que l'aire est positive

END FUNCTION calculate_wetted_area_parabolic 


FUNCTION calculate_porosity(h_b, y1, b1, y_min, b_min, yN, bN) RESULT(poro)
    !=======================================================================
    ! Computes the porosity phi for a given cell based on its asymmetric
    ! parabolic sub-grid cross-section.
    ! H_k is calculated internally from h_b and b_min.
    !
    ! INPUTS:
    !   h_b        : REAL(rp), INTENT(IN) :: Macroscopic water depth in the cell
    !   y1, b1     : REAL(rp), INTENT(IN) :: Left bank point coordinates
    !   y_min, b_min : REAL(rp), INTENT(IN) :: Lowest point coordinates (b_min is also used as b_K)
    !   yN, bN     : REAL(rp), INTENT(IN) :: Right bank point coordinates
    !
    ! OUTPUT:
    !   poro       : REAL(rp) :: The computed porosity for the cell
    !=======================================================================

    IMPLICIT NONE

    ! --- Argument Declarations ---
    REAL(rp), INTENT(IN) :: h_b, y1, b1, y_min, b_min, yN, bN

    ! --- Result Declaration ---
    REAL(rp) :: poro

    ! --- Local Variable Declarations ---
    REAL(rp) :: H_k
    REAL(rp) :: wetted_area
    REAL(rp) :: macro_area
    REAL(rp) :: total_width
    
    ! --- External Function Declaration ---
    REAL(rp) :: calculate_wetted_area_parabolic

    ! --- Start of Algorithm ---

    ! Safety check: for a dry cell, porosity is conventionally 1.0
    IF (h_b < 1.0E-6_rp) THEN
        poro = 1.0_rp
        RETURN
    END IF

    ! 1. Calculate the free surface elevation inside the function.
    !    Assumption: The macroscopic bed elevation (b_K) is the lowest point (b_min).
    H_k = h_b + b_min

    ! 2. Calculate the real wetted area (the numerator)
    wetted_area = calculate_wetted_area_parabolic(H_k, y1, b1, y_min, b_min, yN, bN)

    ! 3. Calculate the macroscopic rectangular area (the denominator)
    total_width = yN - y1
    macro_area = total_width * h_b

    ! 4. Compute the final porosity
    IF (macro_area > 1.0E-9_rp) THEN
        poro = wetted_area / macro_area
    ELSE
        poro = 1.0_rp
    END IF

END FUNCTION calculate_porosity
