module schemas
    use numerics
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


    subroutine flux_GD (Ns, Flux, U_O)
    ! FLUX POUR SCHEMA DE GODUNOV (si f est concave !)
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
            else ! cas d'un choc pour f concave
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
        real(rp), intent(in) :: dt, dx
        integer, intent(in) :: Ns
        real(rp), dimension(Ns), intent(inout) :: Flux
        real(rp), dimension(Ns), intent(in) :: U_O
        real(rp) :: Delta, coef
        integer :: i
        real(rp) :: mid

        Delta = (dx / dt) * 0.5_rp ! on calcule dx/2dt une fois 
        do i = 1,(Ns-1)
            mid = 0.5_rp*(U_O(i) + U_O(i+1))
            if (a_f(mid) /= 0.) then
                coef = a_f(mid)
            else ! points soniques
                coef = (f(U_O(i+1))- f(U_O(i)))/(U_O(i+1) - U_O(i))
            end if
            Flux(i) = 0.5_rp*(f(U_O(i)) + f(U_O(i+1))) - Delta * coef * (f(U_O(i+1)) - f(U_O(i)))
        end do
    end subroutine flux_LW

    subroutine flux_HLL(Ns, Flux, U_O)
    ! FLUX POUR SCHEMA HLL
        integer, intent(in) :: Ns
        real(rp), dimension(Ns), intent(inout) :: Flux
        real(rp), dimension(Ns), intent(in) :: U_O
        integer :: i
        real(rp) :: bl, br

        do i=1,Ns-1
            !bl = min(0._rp, a_f(U_O(i)) - U_O(i))
            !br = max(0._rp, U_O(i+1)+a_f(U_O(i+1)))
            bl = min(a_f(U_O(i)) - U_O(i), a_f(U_O(i+1)) - U_O(i+1))
            br = max(a_f(U_O(i)) + U_O(i), a_f(U_O(i+1)) + U_O(i+1))
            if (bl >= 0) then
                Flux(i) = f(U_O(i))
            else if (bl<0 .AND. br>=0) then
                Flux(i) = (br*f(U_O(i))-bl*f(U_O(i+1))+(bl*br)*(U_O(i+1)-U_O(i)))/(br-bl)
            else if (br < 0) then
                Flux(i) = f(U_O(i+1))
            else
                write(6,*) 'Probleme de calcul de vitesse pour flux HLL'
            end if
        end do
    end subroutine flux_HLL

    real(rp) function minmod(sigma_l, sigma_r)
        real(rp) :: sigma_l, sigma_r
        if (sigma_l >= 0 .AND. sigma_r >= 0.) then
            minmod = min(sigma_l, sigma_r)
        else if (sigma_l < 0 .AND. sigma_r < 0.) then
            minmod = max(sigma_l, sigma_r)
        else
            minmod = 0._rp
        end if
    end function minmod

    subroutine MUSCL(Ns, U_O, U_plus, U_moins, dx)
        integer, intent(in) :: Ns
        real(rp), intent(in) :: dx
        real(rp), dimension(Ns), intent(in) :: U_O 
        real(rp), dimension(Ns), intent(out) :: U_plus
        real(rp), dimension(Ns), intent(out) :: U_moins
        integer :: i
        real(rp) :: sigmal, sigmar, Delta, sigma

        Delta = 1._rp/dx
        do i =2,Ns-1
            sigmal = (U_O(i)-U_O(i-1))*Delta
            sigmar = (U_O(i+1)-U_O(i))*Delta
            sigma = minmod(sigmal, sigmar)
            U_plus(i) = U_O(i) + sigma*0.5_rp*dx
            U_moins(i) = U_O(i) + sigma*0.5_rp*dx
        end do 
    end subroutine MUSCL

end module schemas