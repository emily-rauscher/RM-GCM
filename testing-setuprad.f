      subroutine opacity_wrapper(t_pass, tau_IRe, tau_Ve, Beta_V, Beta_IR, gravity_SI, incident_starlight_fraction)
          include 'rcommons.h'

          integer :: NLAYER, J, k
          real :: mu_0, Tirr, Tint, gravity_SI, incident_starlight_fraction

          real, dimension(NIRP) :: Beta_IR
          real, dimension(NSOLP) :: Beta_V

          real, dimension(NIRP,NLAYER+1) :: tau_IRe
          real, dimension(NSOLP,NLAYER+1) :: tau_Ve
          real, dimension(NLAYER) :: dpe, Pl, Tl, t_pass

          ! Malsky check whether this is NLAYER+1
          real, dimension(NLAYER) :: pe

          mu_0 = incident_starlight_fraction

          Tint = (FBASEFLUX / 5.670367E-8) ** 0.25
          Tirr = (SOLC_IN   / 5.670367E-8) ** 0.25

          do J = 1, NLAYER
             pe(J) = press(J) / 10! convert to pascals
          end do

          DO J = 1, NLAYER-1
              dpe(J) = pe(J+1) - pe(J)
              pl(J) = dpe(J) / log(pe(J+1)/pe(J))
          END DO


          if (tt(1) .ge. 100.0) then
              DO J = 1, NLAYER
                  Tl(J) = tt(J)
              END DO
          else
              DO J = 1, NLAYER
                  Tl(J) = t_pass(J)
              END DO
          end if

          dpe(NLAYER) = 10.0 ** (LOG10(dpe(NLAYER-1)) + (LOG10(dpe(NLAYER-1)) - LOG10(dpe(NLAYER-2))))
          pl(NLAYER)  = 10.0 ** (LOG10(pl(NLAYER-1))  + (LOG10(pl(NLAYER-1))  - LOG10(pl(NLAYER-2))))
          Tl(NLAYER)  = Tl(NLAYER-1) + ABS(Tl(NLAYER-1) - Tl(NLAYER-2)) / 2.0




          CALL calculate_opacities(NLAYER, NSOLP, NIRP, mu_0,Tirr, Tint, Tl, Pl, dpe, tau_IRe,tau_Ve, Beta_V,
     &                            Beta_IR,gravity_SI, with_TiO_and_VO, pe)


      end subroutine opacity_wrapper

      subroutine calculate_opacities(NLAYER, NSOLP, NIRP, mu_0, Tirr, Tint, Tl, Pl, dpe, tau_IRe,tau_Ve,Beta_V,
     &                               Beta_IR,gravity_SI, with_TiO_and_VO, pe)
        ! Input:
        ! Teff - Effective temperature [K] (See Parmentier papers for various ways to calculate this)
        ! for non-irradiated atmosphere Teff = Tint
        ! table_num - Table selection from Parmentier et al. (2015): 1 = w. TiO/VO, 2 = w.o. TiO/VO

        ! Output:
        ! gam_V(3) - gamma ratio for 3 visual bands (gam_V = kV_Ross/kIR_Ross)
        ! beta_V(3) - fraction of total incident stellar flux in band (1/3 for Parmentier values)
        ! Beta - equilvalent bandwidth for picket fence IR model
        ! gam_1 - gamma ratio for IR band 1 (gam_1 = kIR_1/kIR_Ross)
        ! gam_2 - gamma ratio for IR band 2 (gam_2 = kIR_2/kIR_Ross)
        ! gam_P - gamma ratio for Planck mean (gam_P = kIR_Planck/kIR_Ross)
        ! tau_lim - tau limit variable (usually for IC system)

        implicit none
        real :: gam_1, gam_2, tau_lim, gam_P
        real, dimension(NSOLP) :: Beta_V, gam_V
        real, dimension(NIRP) :: Beta_IR
        integer :: k, NLAYER, J, NSOLP, NIRP, i
        real :: Teff, Tint, Tirr, mu_0

        real :: R,gravity_SI
        real :: aP, bP, cP
        real :: aV1, bV1, aV2, bV2, aV3, bV3
        real :: aB, bB
        real :: l10T, l10T2, RT

        real, dimension(NLAYER) :: dpe, Pl, Tl, pe

        real, dimension(NIRP, NLAYER) :: k_IRl
        real, dimension(NSOLP,NLAYER) :: k_Vl

        real, dimension(NIRP,NLAYER+1) :: tau_IRe
        real, dimension(NSOLP,NLAYER+1) :: tau_Ve
        real :: grav
        real :: with_TiO_and_VO
        real :: Bond_Albedo

        grav = gravity_SI

        !! Parmentier opacity profile parameters - first get Bond albedo
        Teff = ((Tint * Tint * Tint * Tint) + (1.0 / sqrt(3.0)) *
     &          (Tirr * Tirr * Tirr * Tirr)) ** (0.25)
        call Bond_Parmentier(Teff, grav, Bond_Albedo)

        !! Recalculate Teff and then find parameters
        Teff = ((Tint * Tint * Tint * Tint) + (1.0 - Bond_Albedo) * mu_0 *
     &          (Tirr * Tirr * Tirr * Tirr)) ** (0.25)

        ! Log 10 T_eff variables
        l10T = log10(Teff)
        l10T2 = l10T * l10T

        if (with_TiO_and_VO .eq. 1) then
          ! First table in Parmentier et al. (2015) w. TiO/VO
          ! Start large if statements with visual band and Beta coefficents
          if (Teff <= 200.0) then
            aV1 = -5.51 ; bV1 = 2.48
            aV2 = -7.37 ; bV2 = 2.53
            aV3 = -3.03 ; bV3 = -0.20
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 200.0) .and. (Teff <= 300.0)) then
            aV1 = 1.23 ; bV1 = -0.45
            aV2 = 13.99 ; bV2 = -6.75
            aV3 = -13.87 ; bV3 = 4.51
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 300.0) .and. (Teff <= 600.0)) then
            aV1 = 8.65 ; bV1 = -3.45
            aV2 = -15.18 ; bV2 = 5.02
            aV3 = -11.95 ; bV3 = 3.74
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 600.0) .and. (Teff <= 1400.0)) then
            aV1 = -12.96 ; bV1 = 4.33
            aV2 = -10.41 ; bV2 = 3.31
            aV3 = -6.97 ; bV3 = 1.94
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 1400.0) .and. (Teff < 2000.0)) then
            aV1 = -23.75 ; bV1 = 7.76
            aV2 = -19.95 ; bV2 = 6.34
            aV3 = -3.65 ; bV3 = 0.89
            aB = 0.84  ; bB = 0.0
          else if (Teff >= 2000.0) then
            aV1 = 12.65 ; bV1 = -3.27
            aV2 = 13.56 ; bV2 = -3.81
            aV3 = -6.02 ; bV3 = 1.61
            aB = 6.21  ; bB = -1.63
          end if

          ! gam_P coefficents
          aP = -2.36
          bP = 13.92
          cP = -19.38

        else if (with_TiO_and_VO .eq. 2)then
          ! Appendix table from Parmentier et al. (2015) - without TiO and VO
          if (Teff <= 200.0) then
            aV1 = -5.51 ; bV1 = 2.48
            aV2 = -7.37 ; bV2 = 2.53
            aV3 = -3.03 ; bV3 = -0.20
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 200.0) .and. (Teff <= 300.0)) then
            aV1 = 1.23 ; bV1 = -0.45
            aV2 = 13.99 ; bV2 = -6.75
            aV3 = -13.87 ; bV3 = 4.51
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 300.0) .and. (Teff <= 600.0)) then
            aV1 = 8.65 ; bV1 = -3.45
            aV2 = -15.18 ; bV2 = 5.02
            aV3 = -11.95 ; bV3 = 3.74
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 600.0) .and. (Teff <= 1400.0)) then
            aV1 = -12.96 ; bV1 = 4.33
            aV2 = -10.41 ; bV2 = 3.31
            aV3 = -6.97 ; bV3 = 1.94
            aB = 0.84  ; bB = 0.0
          else if ((Teff > 1400.0) .and. (Teff < 2000.0)) then
            aV1 = -1.68 ; bV1 = 0.75
            aV2 = 6.96 ; bV2 = -2.21
            aV3 = 0.02 ; bV3 = -0.28
            aB = 3.0  ; bB = -0.69
          else if (Teff >= 2000.0) then
            aV1 = 10.37 ; bV1 = -2.91
            aV2 = -2.4 ; bV2 = 0.62
            aV3 = -16.54 ; bV3 = 4.74
            aB = 3.0  ; bB = -0.69
          end if

          !gam_P coefficents
          if (Teff <= 1400.0) then
            aP = -2.36
            bP = 13.92
            cP = -19.38
          else
            aP = -12.45
            bP = 82.25
            cP = -134.42
          end if
        end if

        ! Calculation of all values
        ! Visual band gamma
        gam_V(1) = 10.0**(aV1 + bV1 * l10T)
        gam_V(2) = 10.0**(aV2 + bV2 * l10T)
        gam_V(3) = 10.0**(aV3 + bV3 * l10T)

        ! Visual band fractions
        Beta_V(:) = 1.0/3.0

        ! gamma_Planck - if < 1 then make it grey approximation (k_Planck = k_Ross, gam_P = 1)
        gam_P = 10.0**(aP * l10T2 + bP * l10T + cP)
        if (gam_P < 1.0000001) then
          gam_P = 1.0000001
        end if

        ! equivalent bandwidth value
        Beta_IR(1) = aB + bB * l10T
        Beta_IR(2) = 1.0 - Beta_IR(1)

        ! IR band kappa1/kappa2 ratio - Eq. 96 from Parmentier & Menou (2014)
        RT = (gam_P - 1.0)/(2.0*Beta_IR(1)*Beta_IR(2))
        R = 1.0 + RT + sqrt(RT * RT + RT)

        ! gam_1 and gam_2 values - Eq. 92, 93 from Parmentier & Menou (2014)
        gam_1 = Beta_IR(1) + R - Beta_IR(1)*R
        gam_2 = gam_1 / R



        tau_Ve(:,1) = 0.0
        tau_IRe(:,1) = 0.0


        do k = 1, NLAYER
          call k_Ross_Freedman(Tl(k), pl(k), 0.0, k_IRl(1,k))
          !call k_Ross_Valencia(Tl(k), pl(k), 0.0, k_IRl(1,k))

          k_Vl(1,k) = k_IRl(1,k) * gam_V(1)
          k_Vl(2,k) = k_IRl(1,k) * gam_V(2)
          k_Vl(3,k) = k_IRl(1,k) * gam_V(3)

          k_IRl(2,k) = k_IRl(1,k) * gam_2
          k_IRl(1,k) = k_IRl(1,k) * gam_1

          tau_Ve(:,k)  = ((k_Vl(:,k)  * dpe(k)) / grav)
          tau_IRe(:,k) = ((k_IRl(:,k) * dpe(k)) / grav)
        end do

      end subroutine calculate_opacities



      subroutine k_Ross_Freedman(Tin, Pin, met, k_IR)
        implicit none

        !! Coefficent parameters for Freedman et al. (2014) table fit
        real, parameter :: c1 = 10.602
        real, parameter :: c2 = 2.882
        real, parameter :: c3 = 6.09e-15
        real, parameter :: c4 = 2.954
        real, parameter :: c5 = -2.526
        real, parameter :: c6 = 0.843
        real, parameter :: c7 = -5.490
        real, parameter :: c8_l = -14.051, c8_h = 82.241
        real, parameter :: c9_l = 3.055, c9_h = -55.456
        real, parameter :: c10_l = 0.024, c10_h = 8.754
        real, parameter :: c11_l = 1.877, c11_h = 0.7048
        real, parameter :: c12_l = -0.445, c12_h = -0.0414
        real, parameter :: c13_l = 0.8321, c13_h = 0.8321

        ! Input:
        ! T - Local gas temperature [K]
        ! P - Local gas pressure [pa]
        ! met - Local metallicity [M/H] (log10 from solar, solar [M/H] = 0.0)

        ! Output:
        ! k_IR - IR band Rosseland mean opacity [m2 kg-1]

        real, intent(in) :: Tin, Pin, met
        real, intent(out) :: k_IR

        real :: k_lowP, k_hiP
        real :: T, P, Tl10, Pl10

        k_IR = 0.0

        T = Tin
        P = Pin * 10.0 ! CoNLAYER to dyne cm-2

        Tl10 = log10(T)
        Pl10 = log10(P)


        ! Low pressure expression
        k_lowP = c1*atan(Tl10 - c2) -
     &    (c3/(Pl10 + c4))*exp((Tl10 - c5)**2) +
     &    c6*met + c7

      ! Temperature split for coefficents = 800 K
        if (T <= 800.0) then
          k_hiP = c8_l + c9_l*Tl10 +
     &    c10_l*Tl10**2 + Pl10*(c11_l + c12_l*Tl10) +
     &    c13_l * met * (0.5 + 0.31830988*atan((Tl10 - 2.5) / 0.2))
        else
          k_hiP = c8_h + c9_h*Tl10 +
     &    c10_h*Tl10**2 + Pl10*(c11_h + c12_h*Tl10) +
     &    c13_h * met * (0.5 + 0.31830988*atan((Tl10 - 2.5) / 0.2))
        end if
        ! Total Rosseland mean opacity - coverted to m2 kg-1
        k_IR = (10.0**k_lowP + 10.0**k_hiP) / 10.0

        ! Avoid divergence in fit for large values
        k_IR = min(k_IR,1.0e30)
      end subroutine k_Ross_Freedman


      subroutine Bond_Parmentier(Teff0, grav,  Bond_Albedo)
        implicit none

        ! Input:
        ! Teff0 - Atmospheric profile effective temperature [K] with zero albedo
        ! grav - Surface gravity of planet [m s-2]

        ! Output:
        ! Bond_Albedo - Bond_Albedo

        real, intent(in) :: Teff0, grav
        real, intent(out) :: Bond_Albedo

        real :: a, b

        ! a and b cofficents dependent on T_eff and grav
        if (Teff0 <= 250.0) then
          a = -0.335 * grav**(0.070)
          b = 0.0
        else if ((Teff0 > 250.0) .and. (Teff0 <= 750.0)) then
          a = -0.335 * grav**(0.070) + 2.149 * grav**(0.135)
          b = -0.896 * grav**(0.135)
        else if ((Teff0 > 750.0) .and. (Teff0 < 1250.0)) then
          a = -0.335 * grav**(0.070) -  0.428 * grav**(0.135)
          b = 0.0
        else if (Teff0 >= 1250.0) then
          a = 16.947 - 3.174 * grav**(0.070) - 4.051 * grav**(0.135)
          b = -5.472 + 0.917 * grav**(0.070) + 1.170 * grav**(0.135)
        end if

        ! Final Bond Albedo expression
        Bond_Albedo = 10.0**(a + b * log10(Teff0))
      end subroutine Bond_Parmentier




      !! Calculates the IR band Rosseland mean opacity (local T) according to the
      !! Valencia et al. (2013) fit and coefficents
      subroutine k_Ross_Valencia(Tin, Pin, met, k_IR)
        implicit none

        ! Input:
        ! T - Local gas temperature [K]
        ! P - Local gas pressure [pa]
        ! met - Local metallicity [M/H] (log10 from solar, solar [M/H] = 0.0)

        ! Output:
        ! k_IR - IR band Rosseland mean opacity [m2 kg-1]

        real Tin, Pin, met
        real k_IR

        real k_lowP, k_hiP
        real Tl10, Pl10
        real T, P

        !! Coefficents parameters for the Valencia et al. (2013) table fit
        real :: c1_v = -37.50
        real :: c2_v = 0.00105
        real :: c3_v = 3.2610
        real :: c4_v = 0.84315
        real :: c5_v = -2.339
        real :: c6_vl = -14.051, c6_vh = 82.241
        real :: c7_vl = 3.055, c7_vh = -55.456
        real :: c8_vl = 0.024, c8_vh = 8.754
        real :: c9_vl = 1.877, c9_vh = 0.7048
        real :: c10_vl = -0.445, c10_vh = -0.0414
        real :: c11_vl = 0.8321, c11_vh = 0.8321
        real :: onedivpi = 1.0 / 3.141592653


        T = Tin
        P = Pin * 10.0 ! Convert to dyne

        Tl10 = log10(T)
        Pl10 = log10(P)

        k_lowP = c1_v * (Tl10-c2_v*Pl10-c3_v)**2 + (c4_v*met + c5_v)

        ! Temperature split for coefficents = 800 K
        if (T <= 800.0) then
          k_hiP = (c6_vl+c7_vl*Tl10+c8_vl*Tl10**2)
     &     + Pl10*(c9_vl+c10_vl*Tl10)
     &     + met*c11_vl*(0.5 + onedivpi*atan((Tl10-2.5)/0.2))
        else
          k_hiP = (c6_vh+c7_vh*Tl10+c8_vh*Tl10**2)
     &     + Pl10*(c9_vh+c10_vh*Tl10)
     &     + met*c11_vh*(0.5 + onedivpi*atan((Tl10-2.5)/0.2))

        end if

        ! Total Rosseland mean opacity - converted to m2 kg-1
        k_IR = (10.0**k_lowP + 10.0**k_hiP) / 10.0

        ! Avoid divergence in fit for large values
        k_IR = min(k_IR,1.0e30)
      end subroutine k_Ross_Valencia
