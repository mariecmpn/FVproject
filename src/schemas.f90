module schemas
    use initialisation_sauvegarde
    IMPLICIT NONE

    contains

    subroutine flux_LF(Ns, Flux, U_O, dt, dx)
    ! FLUX POUR SCHEMA DE LAX-FRIEDRICHS
        IMPLICIT NONE
        real(rp), intent(in) :: dt, dx
        integer, intent(in) :: Ns
        real(rp), dimension(Ns), intent(inout) :: Flux
        real(rp), dimension(Ns), intent(in) :: U_O
        real(rp) :: Delta
        integer :: i

        ! calcul flux pour Lax-Friedrichs
        Delta = (dx / dt) * 0.5_rp
        do i = 1,(Ns-1)
            Flux(i) = 0.5_rp*(f(U_O(i)) + f(U_O(i+1))) - Delta * (U_O(i+1) - U_O(i))
        end do
    end subroutine flux_LF

    subroutine flux_MR(Ns, Flux, U_O)
    ! FLUX POUR SCHEMA DE MURMAN-ROE
        integer, intent(in) :: Ns
        real(rp), dimension(Ns), intent(inout) :: Flux
        real(rp), dimension(Ns), intent(in) :: U_O
        real(rp) :: vitesse
        integer :: i

        do i = 1,(Ns-1)
            if (U_O(i) == U_O(i+1)) then
                vitesse = a_f(U_O(i))
            else
                vitesse = (f(U_O(i+1))- f(U_O(i)))/(U_O(i+1) - U_O(i))
            end if

            Flux(i) = 0.5_rp*(f(U_O(i)) + f(U_O(i+1))) - 0.5_rp*abs(vitesse)*(U_O(i+1) - U_O(i))
        end do
    end subroutine flux_MR

    !subroutine flux_GD (Ns, Flux, U_O)
    ! FLUX POUR SCHEMA DE GODUNOV
    !    integer, intent(in) :: Ns
    !    real(rp), dimension(Ns), intent(inout) :: Flux
    !    real(rp), dimension(Ns), intent(in) :: U_O
    !    real(rp) :: vitesse
    !    integer :: i

    !    do i = 1,(Ns-1)
    !        if (U_O(i) < U_O(i+1)) then ! cas d'une detente pour f convexe
    !            if (a_f(U_O(i)) > 0._rp) then
    !                Flux(i) = f(U_O(i))
    !            else if (a_f(U_O(i+1)) < 0._rp) then
    !                Flux(i) = f(U_O(i+1))
    !            else
    !                Flux(i) = a_inv(0._rp)
    !            end if
    !        else
    !            vitesse = (f(U_O(i+1)) - f(U_O(i))) / (U_O(i+1) - U_O(i))
    !            if (vitesse > 0._rp) then
    !                Flux(i) = f(U_O(i))
    !            else
    !                Flux(i) = f(U_O(i+1))
    !            end if
    !        end if ! cas d'un choc pour f convexe
    !    end do
    !end subroutine flux_GD


    subroutine flux_GD (Ns, Flux, U_O)
        ! FLUX POUR SCHEMA DE GODUNOV
            integer, intent(in) :: Ns
            real(rp), dimension(Ns), intent(inout) :: Flux
            real(rp), dimension(Ns), intent(in) :: U_O
            real(rp) :: vitesse
            integer :: i
    
            do i = 1,(Ns-1)
                if (U_O(i) > U_O(i+1)) then ! cas d'une detente pour f concave
                    if (a_f(U_O(i)) > 0._rp) then
                        Flux(i) = f(U_O(i))
                    else if (a_f(U_O(i+1)) < 0._rp) then
                        Flux(i) = f(U_O(i+1))
                    else
                        Flux(i) = a_inv(0._rp)
                    end if
                else
                    vitesse = (f(U_O(i+1)) - f(U_O(i))) / (U_O(i+1) - U_O(i))
                    if (vitesse > 0._rp) then
                        Flux(i) = f(U_O(i))
                    else
                        Flux(i) = f(U_O(i+1))
                    end if
                end if 
            end do
        end subroutine flux_GD

    subroutine flux_LW(Ns, Flux, U_O, dt, dx)
    ! FLUX POUR SCHEMA LAX-WENDROFF
        IMPLICIT NONE
        real(rp), intent(in) :: dt, dx
        integer, intent(in) :: Ns
        real(rp), dimension(Ns), intent(inout) :: Flux
        real(rp), dimension(Ns), intent(in) :: U_O
        real(rp) :: Delta
        integer :: i
        real(rp) :: mid

        Delta = (dx / dt) * 0.5_rp
        do i = 1,(Ns-1)
            mid = 0.5_rp*(U_O(i) + U_O(i+1))
            Flux(i) = 0.5_rp*(f(U_O(i)) + f(U_O(i+1))) - Delta * a_f(mid) * (U_O(i+1) - U_O(i))
        end do
    end subroutine flux_LW

end module schemas