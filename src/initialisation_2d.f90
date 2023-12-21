module initialisation_sauvegarde_2d
    use numerics
    IMPLICIT NONE

    contains 

    real(rp) function initial_u(x)
        IMPLICIT NONE
        real(rp) :: x
        if (x<0.) then
            initial_u = 3._rp
        else 
            initial_u = 4._rp
        end if
    end function initial_u


    real(rp) function initial_rho(x)
        IMPLICIT NONE
        real(rp) :: x
        if (x<0.) then
            initial_rho = 2._rp
        else 
            initial_rho = 1._rp
        end if
    end function initial_rho


    subroutine lecture_donnees_syst(file_name, x_deb, x_fin, Ns, CFL, T_fin, condition, schema, v_max, rho_max)
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name
        integer, intent(inout) :: Ns
        real(rp), intent(inout) :: x_deb, x_fin
        real(rp), intent(inout) :: CFL
        real(rp), intent(inout) :: T_fin
        character(len = 1), intent(inout) :: condition
        !integer, intent(inout) :: condition, schema
        character(len = 2), intent(inout) :: schema
        real(rp), intent(inout) :: v_max, rho_max

        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')
        
        read(my_unit, *) x_deb, x_fin
        read(my_unit, *) Ns
        read(my_unit, *) CFL
        read(my_unit, *) T_fin
        read(my_unit, *) condition
        read(my_unit, *) schema
        read(my_unit, *) rho_max
        read(my_unit, *) v_max

        !if (condition_lim == 'D') then
        !    condition = 0
        !else if (condition_lim == 'N') then
        !    condition = 1
        !else
        !    condition = -1
        !end if

        if (schema == 'LF') then ! Lax-Friedrichs
            write(6,*) "Schema utilise: Lax_Friedrichs"
        else if (schema == 'RS') then ! Rusanov
            write(6,*) "Schema utilise: Rusanov"
        !else if (schema_use == 'GD') then ! Godunov
        !    schema = 2
        !else if (schema_use == 'LW') then ! Lax-Wendroff
        !    schema = 3
        !else
        !    schema = 0
        end if

        close(my_unit)
    end subroutine lecture_donnees_syst


    subroutine initialisation_syst(W_O, Ns, x_deb, x_fin)
        IMPLICIT NONE
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        real(rp), intent(in) :: x_deb, x_fin
        real(rp) :: x
        integer :: i
        real(rp) :: deltax

        deltax = (x_fin-x_deb)/Ns
        do i = 1,Ns
            x = x_deb + i*deltax
            W_O(1,i) = initial_rho(x)
            W_O(2,i) = initial_u(x)
        end do
    end subroutine initialisation_syst


    subroutine sauvegarde_syst(file_name_u, file_name_rho, W_O, Ns, x_deb, x_fin)
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name_u, file_name_rho
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(in) :: W_O
        real(rp), intent(in) :: x_deb, x_fin
        real(rp) :: x
        integer :: i
        integer :: my_unit_1 = 60, my_unit_2 = 70

        open(my_unit_1, file = file_name_rho, action = 'write', form = 'formatted', status = 'unknown')
        open(my_unit_2, file = file_name_u, action = 'write', form = 'formatted', status = 'unknown')

        do i = 1,Ns
            x = x_deb + i*(x_fin-x_deb)/Ns
            write(my_unit_1, *) x, W_O(1,i)
            write(my_unit_2, *) x, W_O(2,i)
        end do

        close(my_unit_1)
        close(my_unit_2)
    end subroutine sauvegarde_syst


    real(rp) function p(rho, v_max, rho_max)
        ! fonction p
        real(rp) :: rho, v_max, rho_max
        p = v_max*(rho**2/rho_max)
    end function p


    real(rp) function p_prime(rho, v_max, rho_max)
        ! fonction derivee de p
        real(rp) :: rho, v_max, rho_max
        p_prime = 2._rp*v_max*(rho/rho_max)
    end function p_prime


    real(rp) function vitesse(W, v_max, rho_max)
        ! fonction pour la vitesse
        ! au cas ou rho = 0
        real(rp) :: eps = 10.E-14
        real(rp), dimension(2), intent(in) :: W
        real(rp) :: rho, u
        real(rp) :: v_max, rho_max

        rho = W(1)
        if (rho > eps) then 
            u = (W(2) / rho) - p(rho, v_max, rho_max)
        else
            u = 0
        end if
        vitesse = u 
    end function vitesse


    subroutine conserv_to_prim(W_O, Ns, v_max, rho_max)
        ! Subroutine pour passer des variables conservatives aux variables primitives
        ! Q = rho*u + rho*p(rho) devient u
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        real(rp), intent(in) :: v_max, rho_max
        integer :: i
        do i = 1,Ns
            W_O(2,i) = vitesse(W_O(1:2,i), v_max, rho_max)
        end do
    end subroutine conserv_to_prim


    subroutine prim_to_conserv(W_O, Ns, v_max, rho_max)
        ! Subroutine pour passer des variables primitives aux variables conservatives
        ! u devient Q = rho*u + rho*p(rho)
        integer, intent(in) :: Ns
        real(rp), dimension(2,1:Ns), intent(inout) :: W_O
        real(rp), intent(in) :: v_max, rho_max
        integer :: i
        do i = 1,Ns
            W_O(2,i) = W_O(1,i)*W_O(2,i) + W_O(1,i)*p(W_O(1,i), v_max, rho_max)
        end do
    end subroutine prim_to_conserv


    real(rp) function lambda_1(W, v_max, rho_max)
        ! fonction pour la valeur propre 1 du systeme
        real(rp), dimension(2) :: W
        real(rp) :: v_max, rho_max
        lambda_1 = W(2) - W(1)*p_prime(W(1), v_max, rho_max)
    end function lambda_1


    real(rp) function lambda_2(W, v_max, rho_max)
        ! fonction pour la valeur propre 1 du systeme
        real(rp), dimension(2) :: W 
        real(rp) :: v_max, rho_max
        lambda_2 = W(2)
    end function lambda_2


    function mat_A(u, rho, v_max, rho_max)
        ! fonction qui calcule la matrice A
        real(rp), intent(in) :: u, rho, v_max, rho_max
        real(rp), dimension(2,2) :: mat_A
        mat_A(1,1) = u
        mat_A(1,2) = rho
        mat_A(2,1) = 0._rp
        mat_A(2,2) = u - rho*p_prime(rho, v_max, rho_max)
    end function mat_A


    subroutine F(F_ex, W, v_max, rho_max)
        ! subroutine pour calculer la fonction F flux exact
        ! W en variables conservatives
        real(rp), dimension(2), intent(in) :: W
        real(rp), dimension(2), intent(inout) :: F_ex
        real(rp) ,intent(in) ::  v_max, rho_max
        real(rp) :: u
        u = vitesse(W(1:2), v_max, rho_max) ! on calcule u car on est en variables conservatives
        F_ex(:) = u * W(:) ! on multiplie par u car u est 
    end subroutine F


end module initialisation_sauvegarde_2d