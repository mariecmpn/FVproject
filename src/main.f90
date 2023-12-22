program scalaire_trafic
    use numerics
    use initialisation_sauvegarde
    use schemas
    IMPLICIT NONE
    real(rp) :: x_deb, x_fin, x
    integer :: Ns
    real(rp) :: CFL, T_fin, date, dt, dx
    real(rp) :: a, Delta
    integer :: i
    real(rp), dimension(:), allocatable :: U_O
    real(rp), dimension(:), allocatable :: U_N
    real(rp), dimension(:), allocatable :: Flux
    real(rp), dimension(:), allocatable :: U_ex, Err
    !integer :: condition, schema
    character(len = 1) :: condition
    character(len = 2) :: schema
    integer :: fonc
    integer :: Nb_iter = 0
    !character(len = 14) :: namefileout
    !character(len = 68) :: commandline

    write(6,*) '------------------------------------------'
    write(6,*) '--------- Modele trafic routier ----------'
    write(6,*) '---------------- Scalaire ----------------'
    write(6,*) '------------------------------------------'

    ! lecture des donnees du fichier donnees.dat
    call lecture_donnees('donnees.dat', x_deb, x_fin, Ns, CFL, T_fin, condition, schema, fonc)

    ! on calcule le pas en espace
    dx = (x_fin - x_deb)/Ns 

    ! allocation memoire des tableaux
    allocate(U_O(1:Ns), U_N(1:Ns), Flux(1:(Ns-1)), U_ex(1:Ns), Err(1:Ns))

    ! initialisation pour t = 0
    call initialisation(U_O, Ns, x_deb, x_fin, fonc)

    write(6,*) 'Norme L2 de la solution initiale: ', normeL2(U_O,Ns)

    ! boucle en temps
    date = 0._rp ! on initialise la date a 0
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
        else if (condition == 'P') then
            U_N(1) = U_O(1) - Delta*(Flux(1) - Flux(Ns-1))
            U_N(Ns) = U_N(1)
        else ! par defaut on prend des conditions de Neumann
            U_N(1) = U_N(2)
            U_N(Ns) = U_N(Ns-1)
        end if
        !mise a jour
        U_O(1:Ns) = U_N(1:Ns)
        ! On compte le nombre d'iterations
        Nb_iter = Nb_iter + 1
    end do

    ! Calcul solution exacte et de l'erreur
    do i = 1,Ns
        x = x_deb + i*(x_fin-x_deb)/Ns
        U_ex(i) = sol_ex(x, date, fonc)
        Err(i) = U_ex(i) - U_O(i)
    end do

    write(6,*) 'Nombre d iterations', Nb_iter
    write(6,*)
    write(6,*) 'Nombre de cellules: ', Ns
    write(6,*) 'pas en x: ', dx
    write(6,*) 'Erreur L2 entre solution approchee et solution exacte: ', normeL2(Err, Ns)
    write(6,*) 'Norme L2 de la solution approchee: ', normeL2(U_O,Ns)
    write(6,*)

    ! Enregistrement de la solution exacte et de la solution approchee
    !namefileout(1:8) = 'solution'
    !namefileout(9:10) = schema
    !namefileout(11:14) = '.dat'

    !write(6,*) 'Sauvegarde dans le fichier: ', namefileout

    call sauvegarde('solution.dat', U_O, Ns, x_deb, x_fin)
    call sauvegarde('solution_ex.dat', U_ex, Ns, x_deb, x_fin)

    !commandline(1:13) = "gnuplot plot "
    !commandline(14:27) = namefileout
    !commandline(28:68) = " with lines, 'solution_ex.dat' with lines"

    !call execute_command_line(commandline)

    deallocate(U_O, U_N, Flux, U_ex, Err)

end program scalaire_trafic