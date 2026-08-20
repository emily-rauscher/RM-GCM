      SUBROUTINE CALC_V_FALL(P, a, g, rhod, T, V_FALL)
        ! Equations and constants taken from Parmentier et al 2013
        IMPLICIT NONE
        REAL, INTENT(IN)  :: P, a, g, rhod, T
        REAL, INTENT(OUT) :: V_FALL
!       H2 gas properties, SI. eps is the L-J well depth over kB, in K.
        REAL, PARAMETER :: kB  = 1.380649E-23
        REAL, PARAMETER :: d   = 2.827E-10
        REAL, PARAMETER :: mH2 = 3.35E-27
        REAL, PARAMETER :: eps = 59.7
        REAL, PARAMETER :: pi  = 3.14159265358979
        REAL, PARAMETER :: mu  = 2.3 * 1.67E-27
        REAL :: lambda, KN, Beta, rhog, eta

        lambda = kB * T / (sqrt(2.0) * pi * d*d * P)
        KN     = lambda / a
        Beta   = 1.0 + KN * (1.256 + 0.4 * EXP(-1.1 / KN))
        rhog   = P / (kB * T) * mu
!       Ackerman & Marley (2001), eq. B2
        eta    = 5.0 / 16.0 * sqrt(pi * mH2 * kB * T) / (pi * d*d)
     &           * (T / eps)**0.16 / 1.22
!       Stokes law with Cunningham slip; valid at low Reynolds number
        V_FALL = MIN(2.0 * Beta * a*a * g * (rhod - rhog) / (9.0 * eta), 100.0)
      END SUBROUTINE CALC_V_FALL