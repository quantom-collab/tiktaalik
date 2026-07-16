! kvlo4tests.f90
!
! by Adam Freese
! part of the package tiktaalik for GPD evolution
!
! Created August 8, 2024.
!
! All kernels are taken from:
!   Belitsky, Fruend and Mueller
!   Nuclear Physics B 574 (2000) 347-406
!   Belitsky:1999hf
!   arxiv:hep-ph/9912379
!
! This module is only used for accuracy benchmarks.
! It uses a fixed nfl=4.

module kvlo4tests
  use constants, only: CF, CA, TF, pi
  use kernels_lo
  use kernels_nlo

  implicit none
  private

  integer,  parameter, private :: dp = kind(1d0)
  integer,  parameter, private :: nfl = 4

  public :: &
      ! V-type
      & KV0_QQ_pls, KV0_QQ_cst, &
      & KV0_QG_reg, KV0_GQ_reg, &
      & KV0_GG_reg, KV0_GG_pls, KV0_GG_cst, &
      & tKV1_QQ_reg, tKV1_NSp_pls, tKV1_NSp_cst, tKV1_NSm_pls, tKV1_NSm_cst, &
      & tKV1_QG_reg, tKV1_GQ_reg, tKV1_GG_pls, tKV1_GG_cst, &
      & KV0_QG_cst, & ! NEW
      ! A-type
      & tKA0_QQ_pls, tKA0_QQ_cst, &
      & tKA0_QG_reg, tKA0_GQ_reg, &
      & tKA0_GG_pls, tKA0_GG_cst, &
      & tKA1_QQ_reg, tKA1_NSp_pls, tKA1_NSp_cst, tKA1_NSm_pls, tKA1_NSm_cst, &
      & tKA1_QG_reg, tKA1_GQ_reg, tKA1_GG_pls, tKA1_GG_cst, &
      & tKA0_QG_cst ! NEW

  contains

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! QQ kernel pieces

    ! LO, V-type

    function KV0_QQ_pls(x, y, xi) result(K)
        ! redundant
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = K0_qq_pls(x, y, xi)
    end function KV0_QQ_pls

    function KV0_QQ_cst(x, xi) result(K)
        ! redundant
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = K0_qq_cst(x, xi)
    end function KV0_QQ_cst

    ! LO, A-type

    function tKA0_QQ_pls(x, y, xi) result(K)
        ! redundant
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = K0_qq_pls(x, y, xi)
    end function tKA0_QQ_pls

    function tKA0_QQ_cst(x, xi) result(K)
        ! redundant
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = K0_qq_cst(x, xi)
    end function tKA0_QQ_cst

    ! NLO, V-type

    function tKV1_QQ_reg(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KV1_qq_reg(x, y, xi)*real(nfl)
    end function tKV1_QQ_reg

    function tKV1_NSp_pls(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KV1_NSp_pls(x, y, xi) + real(nfl)*KV1_NSp_pls_nfl(x, y, xi)
    end function tKV1_NSp_pls

    function tKV1_NSp_cst(x, xi) result(K)
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KV1_NSp_cst(x, xi) + real(nfl)*KV1_NSp_cst_nfl(x, xi)
    end function tKV1_NSp_cst

    function tKV1_NSm_pls(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KV1_NSm_pls(x, y, xi) + real(nfl)*KV1_NSm_pls_nfl(x, y, xi)
    end function tKV1_NSm_pls

    function tKV1_NSm_cst(x, xi) result(K)
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KV1_NSm_cst(x, xi) + real(nfl)*KV1_NSm_cst_nfl(x, xi)
    end function tKV1_NSm_cst

    ! NLO, A-type

    function tKA1_QQ_reg(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KA1_qq_reg(x, y, xi)*real(nfl)
    end function tKA1_QQ_reg

    function tKA1_NSp_pls(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KV1_NSm_pls(x, y, xi) + real(nfl)*KV1_NSm_pls_nfl(x, y, xi)
    end function tKA1_NSp_pls

    function tKA1_NSp_cst(x, xi) result(K)
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KV1_NSm_cst(x, xi) + real(nfl)*KV1_NSm_cst_nfl(x, xi)
    end function tKA1_NSp_cst

    function tKA1_NSm_pls(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KV1_NSp_pls(x, y, xi) + real(nfl)*KV1_NSp_pls_nfl(x, y, xi)
    end function tKA1_NSm_pls

    function tKA1_NSm_cst(x, xi) result(K)
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KV1_NSp_cst(x, xi) + real(nfl)*KV1_NSp_cst_nfl(x, xi)
    end function tKA1_NSm_cst

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! QG kernel pieces

    ! LO, V-type

    function KV0_QG_reg(x, y, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = real(nfl) * ( KA0_qG_reg(x, y, xi) + KVmA0_qG_reg(x, y, xi) )
    end function KV0_QG_reg

    function KV0_QG_cst(x, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = real(nfl) * KA0_qG_cst(x, xi)
    end function KV0_QG_cst

    ! LO, A-type

    function tKA0_QG_reg(x, y, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = real(nfl) * KA0_qG_reg(x, y, xi)
    end function tKA0_QG_reg

    function tKA0_QG_cst(x, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = real(nfl) * KA0_qG_cst(x, xi)
    end function tKA0_QG_cst

    ! NLO, V-type

    function tKV1_QG_reg(x, y, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = real(nfl) * KV1_qG_reg(x, y, xi)
    end function tKV1_QG_reg

    ! NLO, A-type

    function tKA1_QG_reg(x, y, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = real(nfl) * KA1_qG_reg(x, y, xi)
    end function tKA1_QG_reg

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GQ kernel pieces

    ! LO, V-type

    function KV0_GQ_reg(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KA0_Gq_reg(x, y, xi) + KVmA0_Gq_reg(x, y, xi)
    end function KV0_GQ_reg

    ! LO, A-type

    function tKA0_GQ_reg(x, y, xi) result(K)
        ! kinda redundant ... maybe drop?
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KA0_Gq_reg(x, y, xi)
    end function tKA0_GQ_reg

    ! NLO, V-type

    function tKV1_GQ_reg(x, y, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KV1_Gq_reg(x, y, xi) + real(nfl) * KV1_Gq_reg_nfl(x, y, xi)
    end function tKV1_GQ_reg

    ! NLO, A-type

    function tKA1_GQ_reg(x, y, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KA1_Gq_reg(x, y, xi) + real(nfl) * KA1_Gq_reg_nfl(x, y, xi)
    end function tKA1_GQ_reg

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GG kernel pieces

    ! LO, V-type

    function KV0_GG_pls(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KA0_GG_pls(x, y, xi)
    end function KV0_GG_pls

    function KV0_GG_reg(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KVmA0_GG_reg(x, y, xi)
    end function KV0_GG_reg

    function KV0_GG_cst(x, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KA0_GG_cst(x, xi) + real(nfl)*KA0_GG_nfl(x, xi)
    end function KV0_GG_cst

    ! LO, A-type

    function tKA0_GG_pls(x, y, xi) result(K)
        ! redundant
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KA0_GG_pls(x, y, xi)
    end function tKA0_GG_pls

    function tKA0_GG_reg(x, y, xi) result(K)
        ! redundant
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = 0.0_dp
    end function tKA0_GG_reg

    function tKA0_GG_cst(x, xi) result(K)
        ! NOTICE:
        ! Assumes nfl=4
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KA0_GG_cst(x, xi) + real(nfl)*KA0_GG_nfl(x, xi)
    end function tKA0_GG_cst

    ! NLO, V-type

    function tKV1_GG_pls(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KV1_GG_pls(x, y, xi) + real(nfl)*KV1_GG_pls_nfl(x, y, xi)
    end function tKV1_GG_pls

    function tKV1_GG_cst(x, xi) result(K)
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KV1_GG_cst(x, xi) + real(nfl)*KV1_GG_cst_nfl(x, xi)
    end function tKV1_GG_cst

    ! NLO, A-type

    function tKA1_GG_pls(x, y, xi) result(K)
        real(dp), intent(in) :: x, y, xi
        real(dp) :: K
        !
        K = KA1_GG_pls(x, y, xi) + real(nfl)*KA1_GG_pls_nfl(x, y, xi)
    end function tKA1_GG_pls

    function tKA1_GG_cst(x, xi) result(K)
        real(dp), intent(in) :: x, xi
        real(dp) :: K
        !
        K = KA1_GG_cst(x, xi) + real(nfl)*KA1_GG_cst_nfl(x, xi)
    end function tKA1_GG_cst

end module kvlo4tests
