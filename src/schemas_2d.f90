module schemas_2d
    use initialisation_sauvegarde
    IMPLICIT NONE

    contains

    subroutine flux_LF(Ns, Flux, W_O, dt, dx, v_max, rho_max)
    ! FLUX POUR SCHEMA DE LAX-FRIEDRICHS
        IMPLICIT NONE
        real(rp), intent(in) :: dt, dx, v_max, rho_max
        integer, intent(in) :: Ns
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp), dimension(2) :: F_ex1, F_ex
        real(rp) :: Delta
        integer :: i

        ! calcul flux pour Lax-Friedrichs
        Delta = (dx / dt) * 0.5_rp
        call F(F_ex1, W_O(:,1), v_max, rho_max)
        do i = 1,(Ns-1)
            call F(F_ex, W_O(:,(i+1)), v_max, rho_max)
            Flux(1,i) = 0.5_rp*(F_ex(1) + F_ex1(1)) - Delta * (F_ex(1) - F_ex1(1))
            Flux(2,i) = 0.5_rp*(F_ex(2) + F_ex1(2)) - Delta * (F_ex(2) - F_ex1(2))
            F_ex1(:) = F_ex(:)
        end do
    end subroutine flux_LF


    subroutine flux_GD (Ns, Flux, W_O)
    ! FLUX POUR SCHEMA DE GODUNOV
        integer, intent(in) :: Ns
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp) :: v
        integer :: i,j

            !do j = 1,2
            !    do i = 1,(Ns-1)
            !        if (U_O(i) > U_O(i+1)) then ! cas d'une detente pour f concave
            !        if (a_f(U_O(i)) > 0._rp) then
            !                Flux(i) = f(U_O(i))
            !            else if (a_f(U_O(i+1)) < 0._rp) then
            !                Flux(i) = f(U_O(i+1))
            !            else
            !                Flux(i) = a_inv(0._rp)
            !            end if
            !        else
            !            v = (f(U_O(i+1)) - f(U_O(i))) / (U_O(i+1) - U_O(i))
            !            if (v > 0._rp) then
            !                Flux(i) = f(U_O(i))
            !            else
            !                Flux(i) = f(U_O(i+1))
            !            end if
            !        end if 
            !    end do
            !end do
    end subroutine flux_GD

    subroutine flux_RS (Ns, Flux, W_O, v_max, rho_max)
    ! FLUX POUR LE SCHEMA DE RUSANOV
        IMPLICIT NONE
        real(rp), intent(in) :: v_max, rho_max
        integer, intent(in) :: Ns
        real(rp), dimension(2,Ns), intent(inout) :: Flux
        real(rp), dimension(2,Ns), intent(in) :: W_O
        real(rp), dimension(2) :: F_ex1, F_ex
        real(rp) :: Delta
        integer :: i

        call F(F_ex1, W_O(:,1), v_max, rho_max)
        do i = 1,(Ns-1)
            call F(F_ex, W_O(:,(i+1)), v_max, rho_max)
            ! max des valeurs propres
            Delta = max(abs(W_O(2,i)/W_O(1,i)-W_O(1,i)*p(W_O(1,i), v_max, rho_max)), abs(W_O(2,i)/W_O(1,i)-W_O(1,i)*p(W_O(1,i), v_max, rho_max) - W_O(1,i)*p_prime(W_O(1,i), v_max, rho_max)))
            Delta = max(Delta, abs(W_O(2,i+1)/W_O(1,i+1)-W_O(1,i+1)*p(W_O(1,i+1), v_max, rho_max)), abs(W_O(2,i+1)/W_O(1,i+1)-W_O(1,i+1)*p(W_O(1,i+1), v_max, rho_max) - W_O(1,i+1)*p_prime(W_O(1,i+1), v_max, rho_max)))
            Flux(1,i) = 0.5_rp*(F_ex(1) + F_ex1(1)) - 0.5_rp*Delta * (F_ex(1) - F_ex1(1))
            Flux(2,i) = 0.5_rp*(F_ex(2) + F_ex1(2)) - Delta * (F_ex(2) - F_ex1(2))
            F_ex1(:) = F_ex(:)
        end do

    end subroutine flux_RS


end module schemas_2d