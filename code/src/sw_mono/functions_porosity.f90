MODULE fonctions_porosite_mod
    USE m_common
    USE m_mesh
    USE m_model

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: update_all_porosities

CONTAINS
	
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

    
	SUBROUTINE find_section(mesh, target_x, y1, yN)
		!=======================================================================
		! Analyse une section pour trouver les positions y des berges (y1, yN).
		!=======================================================================

		! --- Arguments ---
		TYPE(msh), INTENT(IN) :: mesh
		REAL(rp), INTENT(IN)       :: target_x
		REAL(rp), INTENT(OUT)      :: y1, yN

		! --- Variables Locales ---
		INTEGER  :: inode
		LOGICAL  :: first_point_found = .FALSE.
		REAL(rp), PARAMETER :: tolerance = 1.0E-6_rp

		! --- Boucle unique sur tous les nœuds du maillage ---
		DO inode = 1, mesh%nn
		    ! On ne considère que les nœuds qui appartiennent à la section
		    IF (ABS(mesh%node(inode)%coord%x - target_x) < tolerance) THEN

		        ! Si c'est le premier point qu'on trouve pour cette section
		        IF (.NOT. first_point_found) THEN
		            y1 = mesh%node(inode)%coord%y
		            yN = y1
		            first_point_found = .TRUE.
		        END IF

		        ! Mettre à jour les limites y des berges
		        y1 = MIN(y1, mesh%node(inode)%coord%y)
		        yN = MAX(yN, mesh%node(inode)%coord%y)
		    END IF
		END DO

	END SUBROUTINE find_section




	SUBROUTINE update_all_porosities(dof, mesh)
		!=======================================================================
		! Orchestre la mise à jour de la porosité pour toutes les cellules 1D-like.
		!=======================================================================
		
		! --- Arguments ---
		TYPE(unk), INTENT(IN)    :: dof
		TYPE(msh), INTENT(IN) :: mesh

		! --- Variables Locales ---
		INTEGER  :: inode, icell
		REAL(rp) :: h_b, b_min, H_k, phi_K_new, wetted_area
		REAL(rp) :: y1, y_min, yN, total_width, macro_area
		REAL(rp) :: min_dist_to_b_min
		REAL(rp), PARAMETER :: tolerance = 1.0E-6_rp

		! --- Boucle principale sur toutes les cellules/sections ---
		DO icell = 1, mesh%nc

		    ! 1. Récupérer les données macroscopiques et DÉFINIR b_min
		    h_b = dof%h(icell)
		    b_min = bathy_cell(icell)
		    H_k = h_b + b_min

		    ! 2. Trouver y1, yN, et le y_min correspondant à b_min pour cette section
		    CALL find_section(mesh, mesh%cell(icell)%grav%x, y1, yN)

		    ! Boucle supplémentaire pour trouver le y_min associé à b_min
		    min_dist_to_b_min = HUGE(0.0_rp)
		    y_min = (y1 + yN) / 2.0_rp ! Valeur par défaut au centre
		    DO inode = 1, mesh%nn
		        IF (ABS(mesh%node(inode)%coord%x - mesh%cell(icell)%grav%x) < tolerance) THEN
		            IF (ABS(bathy_node(inode) - b_min) < min_dist_to_b_min) THEN
		                min_dist_to_b_min = ABS(bathy_node(inode) - b_min)
		                y_min = mesh%node(inode)%coord%y
		            END IF
		        END IF
		    END DO

		    ! 3. Calculer l'aire mouillée avec le modèle parabolique
		    wetted_area = calculate_wetted_area_parabolic(H_k, y1, y_min, yN, b_min)

		    ! 4. Calculer la porosité finale
		    total_width = yN - y1
		    macro_area = total_width * h_b

		    IF (macro_area > 1.0E-9_rp) THEN
		        phi_K_new = wetted_area / macro_area
		    ELSE
		        phi_K_new = 1.0_rp 
		    END IF

		    ! 5. Stocker la nouvelle porosité dans le tableau global
		    SPorosity%phi(icell) = phi_K_new

		END DO

	END SUBROUTINE update_all_porosities

END MODULE fonctions_porosite_mod
