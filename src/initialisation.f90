module initialisation_sauvegarde
    use numerics
    IMPLICIT NONE

    contains 

    real(rp) function choc_det(x)
        real(rp) :: x
        if (x>=0. .AND. x<100.) then
            choc_det = 1._rp
        else if (x>=100. .AND. x<500.) then
            choc_det = 2._rp
        else 
            choc_det = 3._rp
        end if 
    end function choc_det

    real(rp) function initial(x)
        IMPLICIT NONE
        real(rp) :: x
        if (x>=0.) then
            initial = 1._rp
        !else if (x>=1. .AND. x<3.) then
        !    initial = 2._rp
        else 
            initial = 0.5_rp
        end if 
    end function initial


    subroutine lecture_donnees(file_name, x_deb, x_fin, Ns, CFL, T_fin, condition, schema)
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name
        integer, intent(inout) :: Ns
        real(rp), intent(inout) :: x_deb, x_fin
        real(rp), intent(inout) :: CFL
        real(rp), intent(inout) :: T_fin
        character(len = 1), intent(inout) :: condition
        character(len = 2), intent(inout) :: schema

        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')
        
        read(my_unit, *) x_deb, x_fin
        read(my_unit, *) Ns
        read(my_unit, *) CFL
        read(my_unit, *) T_fin
        read(my_unit, *) condition
        read(my_unit, *) schema

        if (schema == 'LF') then ! Lax-Friedrichs
            write(6,*) 'Resolution par le schema de Lax-Friedrichs'
        else if (schema == 'MR') then ! Murman-Roe
            write(6,*) 'Resolution par le schema de Murman-Roe'
        else if (schema == 'GD') then ! Godunov
            write(6,*) 'Resolution par le schema de Godunov'
        else if (schema == 'LW') then ! Lax-Wendroff
            write(6,*) 'Resolution par le schema de Lax-Wendroff'
        else if (schema == 'HL') then ! HLL
            write(6,*) 'Resolution par le schema HLL'
        end if
        close(my_unit)
    end subroutine lecture_donnees

    subroutine initialisation(U_O, Ns, x_deb, x_fin)
        IMPLICIT NONE
        integer, intent(in) :: Ns
        real(rp), dimension(1:Ns), intent(inout) :: U_O
        real(rp), intent(in) :: x_deb, x_fin
        real(rp) :: x
        integer :: i
        real(rp) :: delta

        delta = (x_fin-x_deb)/Ns
        do i = 1,Ns
            x = x_deb + i*delta
            U_O(i) = initial(x)
        end do
    end subroutine initialisation

    subroutine sauvegarde(file_name, U_O, Ns, x_deb, x_fin)
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name
        integer, intent(in) :: Ns
        real(rp), dimension(1:Ns), intent(in) :: U_O
        real(rp), intent(in) :: x_deb, x_fin
        real(rp) :: x
        integer :: i
        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'write', form = 'formatted', status = 'unknown')

        do i = 1,Ns
            x = x_deb + i*(x_fin-x_deb)/Ns
            write(my_unit, *) x, U_O(i)
        end do

        close(my_unit)
    end subroutine sauvegarde

    real(rp) function f(x) ! f pour trafic routier
        real(rp) :: x
        f = x*(1-x)
    end function f

    real(rp) function a_f(x) ! a = f' pour trafic routier
        real(rp) :: x
        a_f = 1._rp - 2._rp*x
    end function a_f

    real(rp) function a_inv(x) ! a^{-1} pour trafic routier
        real(rp) :: x
        a_inv = 0.5_rp*(1-x)
    end function a_inv

end module initialisation_sauvegarde