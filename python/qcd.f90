! qcd.f90
!
! by Adma Freese
! part of the package tiktaalik
!
! wrappers for f2py to access

module dummy
  use alpha_qcd

  implicit none
  public

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Running coupling

    subroutine alpha_wrap(nQ2, Q2, alpha)
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nQ2
        real(dp), intent(in)  :: Q2(nQ2)
        real(dp), intent(out) :: alpha(nQ2)
        !
        integer :: i
        do i=1, nQ2, 1
          alpha(i) = get_alpha_qcd(Q2(i))
        end do
    end subroutine alpha_wrap

end module dummy
