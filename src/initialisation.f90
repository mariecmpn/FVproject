module initialisation_sauvegarde
    use numerics
    IMPLICIT NONE

    contains 

    real(rp) function choc_det0(x)
        real(rp) :: x
        if (x<1.) then
            choc_det0 = 1._rp
        else if (x>=1. .AND. x<2.) then
            choc_det0 = 3._rp
        else 
            choc_det0 = 2._rp
        end if 
    end function choc_det0

    real(rp) function choc_det(x,t)
        real(rp) :: x,t

        if (t< 0.5) then 
            if (x < -3.*t+1) then
                choc_det = 1._rp
            else if (x>=-3.*t+1 .AND. x<-5.*t+2) then
                choc_det = 3._rp
            else if (x>=-5.*t+2 .AND. x<-3.*t+2) then 
                choc_det = 0.5_rp*(1 - (x-2)/t)
            else
                choc_det = 2._rp
            end if 
        else if (t>=0.5 .AND. t<0.25*(3*sqrt(5.)+7)) then
            if (x < -3.*t+1) then
                choc_det = 1._rp
            else if (x>=-3.*t+1 .AND. x<t-3*sqrt(2*t)+2) then
                choc_det = 0.5_rp*(1 - (x-2)/t)
            else
                choc_det = 2._rp
            end if 
        else 
            if (x<-2.*t+2) then
                choc_det = 1._rp
            else 
                choc_det = 2._rp
            end if 
        end if
    end function choc_det

    real(rp) function det_choc0(x)
        real(rp) :: x
        if (x<0.) then 
            det_choc0 = 3._rp
        else if (x>=0. .AND. x<1.) then
            det_choc0 = 1._rp
        else 
            det_choc0 = 2._rp
        end if 
    end function det_choc0

    real(rp) function det_choc(x,t)
        real(rp) :: x,t

        if (t<1.) then
            if (x<-5.*t) then
                det_choc = 3._rp
            else if (x>=-5.*t .AND. x<-t) then 
                det_choc = 0.5_rp*(1-(x/t))
            else if (x>=-t .AND. x<-2.*t+1) then
                det_choc = 1._rp
            else
                det_choc = 2._rp
            end if 
        else
            if (x<-5.*t) then
                det_choc = 3._rp
            else if (x>=-5.*t .AND. x<-3.*t+4.*sqrt(t)) then
                det_choc = 0.5_rp*(1-(x/t))
            else 
                det_choc = 2._rp
            end if 
        end if 
    end function det_choc

    real(rp) function initial(x, fonc)
        real(rp) :: x
        integer :: fonc
        if (fonc == 1) then
            initial = choc_det0(x)
        else if (fonc == 2) then
            initial = det_choc0(x)
        end if 
    end function initial

    real(rp) function sol_ex(x, t, fonc)
        real(rp) :: x,t
        integer :: fonc

        if (fonc == 1) then
            sol_ex = choc_det(x,t)
        else if (fonc == 2) then
            sol_ex = det_choc(x,t)
        end if 
    end function sol_ex


    subroutine lecture_donnees(file_name, x_deb, x_fin, Ns, CFL, T_fin, condition, schema, fonc)
        IMPLICIT NONE
        character(len = *), intent(in) :: file_name
        integer, intent(inout) :: Ns
        real(rp), intent(inout) :: x_deb, x_fin
        real(rp), intent(inout) :: CFL
        real(rp), intent(inout) :: T_fin
        character(len = 1), intent(inout) :: condition
        character(len = 2), intent(inout) :: schema
        integer, intent(inout) :: fonc

        integer :: my_unit

        open(newunit = my_unit, file = file_name, action = 'read', form = 'formatted', status = 'old')
        
        read(my_unit, *) x_deb, x_fin
        read(my_unit, *) Ns
        read(my_unit, *) CFL
        read(my_unit, *) T_fin
        read(my_unit, *) condition
        read(my_unit, *) schema
        read(my_unit, *) fonc

        if (fonc == 1) then 
            write(6,*) 'Fonction initiale choisie : '
            write(6,*) '1 si x<1'
            write(6,*) '3 si 1<x<2'
            write(6,*) '2 si x>2'
        else if (fonc == 2) then
            write(6,*) 'Fonction initiale choisie : '
            write(6,*) '3 si x<0'
            write(6,*) '1 si 0<x<1'
            write(6,*) '2 si x>1'
        end if 

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

    subroutine initialisation(U_O, Ns, x_deb, x_fin, fonc)
        IMPLICIT NONE
        integer, intent(in) :: Ns
        real(rp), dimension(1:Ns), intent(inout) :: U_O
        real(rp), intent(in) :: x_deb, x_fin
        integer, intent(in) :: fonc
        real(rp) :: x
        integer :: i
        real(rp) :: delta

        delta = (x_fin-x_deb)/Ns
        do i = 1,Ns
            x = x_deb + i*delta
            U_O(i) = initial(x, fonc)
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