! testing.f90
!
! by Adma Freese
! part of the package tiktaalik
!
! wrappers for f2py to access
!
! These routines are just for me to test things,
! and will probably be removed in a public release.
! Some of what's here may be replaced by a validation module.

module dummy
  use continuum_evolution
  use continuum_wilson
  use gk
  use pixelation

  implicit none
  public

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Test CFFs

    subroutine test_cff_q(nxi, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nxi
        real(dp), intent(in)  :: xi(nxi), Q2
        logical,  intent(in)  :: nlo
        complex(dp), intent(out) :: v(nxi)
        !
        integer :: ixi
        !$OMP PARALLEL DO
        do ixi=1, nxi, 1
          v(ixi) = continuum_cff_q(proxy_function, xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = 4./9.*(Hu(z, xi(ixi), 0.0_dp) - Hu(-z, xi(ixi), 0.0_dp)) &
              & + 1./9.*(Hd(z, xi(ixi), 0.0_dp) - Hd(-z, xi(ixi), 0.0_dp)) &
              & + 1./9.*(Hs(z, xi(ixi), 0.0_dp) - Hs(-z, xi(ixi), 0.0_dp))
          end function proxy_function
    end subroutine test_cff_q

    subroutine test_cff_g(nxi, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nxi
        real(dp), intent(in)  :: xi(nxi), Q2
        logical,  intent(in)  :: nlo
        complex(dp), intent(out) :: v(nxi)
        !
        integer :: ixi
        !$OMP PARALLEL DO
        do ixi=1, nxi, 1
          v(ixi) = continuum_cff_g(proxy_function, xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hg(z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_cff_g

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Test shifts (V-type)

    subroutine test_shift_cNS(nx, nxi, x, xi, Q2, nlo, nstype, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nstype
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_QQ(proxy_function, x(ix), xi(ixi), Q2, nlo, nstype)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hu(z, xi(ixi), 0.0_dp) - Hu(-z, xi(ixi), 0.0_dp) &
              & - Hd(z, xi(ixi), 0.0_dp) + Hd(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cNS

    subroutine test_shift_cQQ(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_QQ(proxy_function, x(ix), xi(ixi), Q2, nlo, 0)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hu(z, xi(ixi), 0.0_dp) - Hu(-z, xi(ixi), 0.0_dp) &
              & + Hd(z, xi(ixi), 0.0_dp) - Hd(-z, xi(ixi), 0.0_dp) &
              & + Hs(z, xi(ixi), 0.0_dp) - Hs(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cQQ

    subroutine test_shift_cQG(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_QG(proxy_function, x(ix), xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hg(z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cQG

    subroutine test_shift_cGQ(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_GQ(proxy_function, x(ix), xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hu(z, xi(ixi), 0.0_dp) - Hu(-z, xi(ixi), 0.0_dp) &
              & + Hd(z, xi(ixi), 0.0_dp) - Hd(-z, xi(ixi), 0.0_dp) &
              & + Hs(z, xi(ixi), 0.0_dp) - Hs(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cGQ

    subroutine test_shift_cGG(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_GG(proxy_function, x(ix), xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hg(z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cGG

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Test shifts (A-type)

    subroutine test_shift_cNS_tilde(nx, nxi, x, xi, Q2, nlo, nstype, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi, nstype
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_QQ_tilde(proxy_function, x(ix), xi(ixi), Q2, nlo, nstype)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hu_tilde(z, xi(ixi), 0.0_dp) + Hu_tilde(-z, xi(ixi), 0.0_dp) &
              & - Hd_tilde(z, xi(ixi), 0.0_dp) - Hd_tilde(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cNS_tilde

    subroutine test_shift_cQQ_tilde(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_QQ_tilde(proxy_function, x(ix), xi(ixi), Q2, nlo, 0)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hu_tilde(z, xi(ixi), 0.0_dp) - Hu_tilde(-z, xi(ixi), 0.0_dp) &
              & + Hd_tilde(z, xi(ixi), 0.0_dp) - Hd_tilde(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cQQ_tilde

    subroutine test_shift_cQG_tilde(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_QG_tilde(proxy_function, x(ix), xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              ! need a stand-in since Hg_tilde=0 in GK model
              ! the following has the right symmetry
              f = Hu_tilde(z, xi(ixi), 0.0_dp) + Hu_tilde(-z, xi(ixi), 0.0_dp) &
              & + Hd_tilde(z, xi(ixi), 0.0_dp) + Hd_tilde(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cQG_tilde

    subroutine test_shift_cGQ_tilde(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_GQ_tilde(proxy_function, x(ix), xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              f = Hu_tilde(z, xi(ixi), 0.0_dp) - Hu_tilde(-z, xi(ixi), 0.0_dp) &
              & + Hd_tilde(z, xi(ixi), 0.0_dp) - Hd_tilde(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cGQ_tilde

    subroutine test_shift_cGG_tilde(nx, nxi, x, xi, Q2, nlo, v)
        ! Test how the shift operates on a proxy function
        integer,  parameter   :: dp = kind(1d0)
        integer,  intent(in)  :: nx, nxi
        real(dp), intent(in)  :: x(nx), xi(nxi), Q2
        logical,  intent(in)  :: nlo
        real(dp), intent(out) :: v(nx,nxi)
        !
        integer :: ix, ixi
        do ixi=1, nxi, 1
        !$OMP PARALLEL DO
        do ix=1, nx, 1
          v(ix,ixi) = continuum_shift_GG_tilde(proxy_function, x(ix), xi(ixi), Q2, nlo)
        end do
        !$OMP END PARALLEL DO
        end do
        contains
          function proxy_function(z) result(f)
              real(dp), intent(in) :: z
              real(dp) :: f
              ! need a stand-in since Hg_tilde=0 in GK model
              ! the following has the right symmetry
              f = Hu_tilde(z, xi(ixi), 0.0_dp) + Hu_tilde(-z, xi(ixi), 0.0_dp) &
              & + Hd_tilde(z, xi(ixi), 0.0_dp) + Hd_tilde(-z, xi(ixi), 0.0_dp)
          end function proxy_function
    end subroutine test_shift_cGG_tilde

end module dummy
