      SUBROUTINE CALC_V_FALL(P, a, g, rhod, T, V_FALL)
        ! Equations and constants taken from Parmentier et al 2013
        IMPLICIT NONE
        REAL, INTENT(IN)  :: P, a, g, rhod, T
        REAL, INTENT(OUT) :: V_FALL
        REAL :: lambda, rhog, eta

        CALL SED_GAS(P, T, lambda, rhog, eta)
        CALL V_FALL_GRAIN(a, g, rhod, lambda, rhog, eta, V_FALL)
      END SUBROUTINE CALC_V_FALL

      SUBROUTINE SED_GAS(P, T, lambda, rhog, eta)
!       Grain-independent fall velocity stuff: mean free path,
!       density and viscosity of the H2 gas at P (Pa) and T (K). Split
!       out so a column can compute them once per level and share them
!       across every grain size.
        IMPLICIT NONE
        REAL, INTENT(IN)  :: P, T
        REAL, INTENT(OUT) :: lambda, rhog, eta
!       H2 gas properties, SI. eps is the L-J well depth over kB, in K.
        REAL, PARAMETER :: kB  = 1.380649E-23
        REAL, PARAMETER :: d   = 2.827E-10
        REAL, PARAMETER :: mH2 = 3.35E-27
        REAL, PARAMETER :: eps = 59.7
        REAL, PARAMETER :: pi  = 3.14159265358979
        REAL, PARAMETER :: mu  = 2.3 * 1.67E-27

        lambda = kB * T / (sqrt(2.0) * pi * d*d * P)
        rhog   = P / (kB * T) * mu
!       Ackerman & Marley (2001), eq. B2
        eta    = 5.0 / 16.0 * sqrt(pi * mH2 * kB * T) / (pi * d*d)
     &           * (T / eps)**0.16 / 1.22
      END SUBROUTINE SED_GAS

      SUBROUTINE V_FALL_GRAIN(a, g, rhod, lambda, rhog, eta, V_FALL)
!       The grain-dependent half: fall speed of a grain of radius a (m)
!       and density rhod.
        IMPLICIT NONE
        REAL, INTENT(IN)  :: a, g, rhod, lambda, rhog, eta
        REAL, INTENT(OUT) :: V_FALL
        REAL :: KN, Beta

        KN     = lambda / a
        Beta   = 1.0 + KN * (1.256 + 0.4 * EXP(-1.1 / KN))
!       Stokes law with Cunningham slip; valid at low Reynolds number
        V_FALL = MIN(2.0 * Beta * a*a * g * (rhod - rhog) / (9.0 * eta), 100.0)
      END SUBROUTINE V_FALL_GRAIN
