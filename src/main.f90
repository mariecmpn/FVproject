program scalaire_trafic
    use numerics
    use initialisation_sauvegarde
    use schemas
    IMPLICIT NONE
    real(rp) :: x_deb, x_fin
    integer :: Ns
    real(rp) :: CFL, T_fin, date, dt, dx
    real(rp) :: a, Delta
    integer :: i
    real(rp), dimension(:), allocatable :: U_O
    real(rp), dimension(:), allocatable :: U_N
    real(rp), dimension(:), allocatable :: Flux
    !integer :: condition, schema
    character(len = 1) :: condition
    character(len = 2) :: schema
    integer :: Nb_iter = 0

    ! lecture des donnees du fichier donnees.dat
    call lecture_donnees('donnees.dat', x_deb, x_fin, Ns, CFL, T_fin, condition, schema)

    dx = (x_fin - x_deb)/Ns

    ! allocation memoire des tableaux
    allocate(U_O(1:Ns), U_N(1:Ns), Flux(1:(Ns-1)))

    ! initialisation pour t = 0
    call initialisation(U_O, Ns, x_deb, x_fin)

    ! boucle en temps
    date = 0._rp
    do while (date < T_fin)
        ! CFL
        a = U_O(1)
        do i = 2,Ns
            a = max(a, U_O(i))
        end do
        dt = dx*CFL/a
        dt = min(dt, T_fin - date)
        date = date + dt

        ! calcul des flux
        if (schema == 'LF') then
            call flux_LF(Ns, Flux, U_O, dt, dx)
        else if (schema == 'MR') then
            call flux_MR(Ns, Flux, U_O)
        else if (schema == 'GD') then
            call flux_GD(Ns, Flux, U_O)
        else if (schema == 'LW') then
            call flux_LW(Ns, Flux, U_O, dt, dx)
        else if (schema == 'HL') then
            call flux_HLL(Ns, Flux, U_O)
        end if

        ! update calcul de u_i^{n+1}
        Delta = dt/dx
        do i = 2,(Ns-1)
            U_N(i) = U_O(i) - Delta*(Flux(i) - Flux(i-1))
        end do

        ! Conditions aux limites
        if (condition == 'D') then 
            ! Dirichlet
            U_N(1) = U_O(1)
            U_N(Ns) = U_O(Ns)
        else ! par defaut on prend des conditions de Neumann
            U_N(1) = U_N(2)
            U_N(Ns) = U_N(Ns-1)
        end if
        
        !mise a jour
        U_O(1:Ns) = U_N(1:Ns)

        Nb_iter = Nb_iter + 1
    end do

    write(6,*) 'Nombre d iterations', Nb_iter
    call sauvegarde('solution.dat', U_O, Ns, x_deb, x_fin)

    deallocate(U_O, U_N, Flux)

end program scalaire_trafic