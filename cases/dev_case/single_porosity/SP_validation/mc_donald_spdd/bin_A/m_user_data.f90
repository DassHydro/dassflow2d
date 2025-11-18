real(kind=8) function bathy_user(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y

    ! Exemple : lit un relief simple en pente
    bathy_user = 0.01 * x   ! pente linéaire en x
    ! Tu peux remplacer par une expression ou interpolation à partir d'un fichier
end function bathy_user

real(kind=8) function manning_user(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y
    manning_user = 0.0d0
end function manning_user

real(kind=8) function zs0_user(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y

    ! η = h + zb, exemple : eau plate initiale à 1 m au-dessus du lit
    zs0_user = bathy_user(x, y) + 1.0d0
end function zs0_user

real(kind=8) function u0_user(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y

    ! Exemple : eau initialement au repos
    u0_user = 0.0d0
end function u0_user

real(kind=8) function v0_user(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y

    ! Eau initialement au repos
    v0_user = 0.0d0
end function v0_user

real(kind=8) function inflow_user(t, x, y)
    implicit none
    real(kind=8), intent(in) :: t, x, y

    ! Exemple : débit constant Q = 1 m^3/s sur toute la frontière
    ! Si tu connais la largeur L, tu peux calculer vitesse u = Q / L
    inflow_user = 2.0d0   ! valeur à adapter
end function inflow_user

  !---------------------------------------------------
  ! Exact solution (optionnel, pour validation)
  !---------------------------------------------------
  real(kind=8) function zs_exact(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y

    zs_exact = bathy_user(x, y) + 1.0d0  ! si eau plate
  end function zs_exact

  real(kind=8) function u_exact(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y

    u_exact = 0.0d0
  end function u_exact

  real(kind=8) function v_exact(x, y)
    implicit none
    real(kind=8), intent(in) :: x, y

    v_exact = 0.0d0
  end function v_exact