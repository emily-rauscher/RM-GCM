      subroutine opacity_wrapper(t_pass, tau_IRe, tau_Ve, Beta_V, Beta_IR, gravity_SI)
          include 'rcommons.h'

          integer :: NLAYER, J, k
          real :: mu_0, Tirr, Tint, gravity_SI

          real, dimension(NIRP) :: Beta_IR
          real, dimension(NSOLP) :: Beta_V

          real, dimension(NIRP,NLAYER+1) :: tau_IRe
          real, dimension(NSOLP,NLAYER+1) :: tau_Ve
          real, dimension(NLAYER) :: dpe, Pl, Tl, t_pass

          real, dimension(NLAYER + 1) :: pe
          !real, dimension(NIRP,54) :: tau_IRe
          !real, dimension(NSOLP,54) :: tau_Ve
          !real, dimension(54) :: dpe, Tl, Pl
          !real, dimension(55) :: pe

          ! This is to calculate the incident fraction of the starlight
          mu_0 = COS(ALAT * PI / 180.0) * COS(ALON * PI / 180.0)

          IF (mu_0 .le. 0) THEN
              mu_0 = 0
          END IF

          ! DOUBLE CHECK THESE VALUES
          Tint = (FBASEFLUX / 5.670367E-5) ** 0.25
          Tirr = (SOLC_IN   / 5.670367E-5) ** 0.25

          do J = 1, NLAYER
             pe(J) = press(J) / 10 ! convert to pascals
          end do
          pe(NLAYER + 1) = pe(NLAYER) + ABS(pe(NLAYER) - pe(NLAYER-1))

          DO J = 1, NLAYER
              dpe(J) = pe(J+1) - pe(J)
          END DO

          if (tt(1) .ge. 100.0) then
              DO J = 1, NZ
                  Tl(J) = tt(J)
              END DO
          else
              DO J = 1, NLAYER
                  Tl(J) = t_pass(J)
              END DO
          end if

!         Tl = (/718.78183337991857, 718.78196091148016,718.78214656096566,718.78239803706924,718.78273476562697,
!     &           718.78319220410356,718.78383172427198,718.78475875694073,718.78615568539703,718.78834193123851,
!     &           718.79188473388683,718.79780493261580,718.80796128319946,718.82577081906561,718.85756232909739,
!     &           718.91512323959023,719.02049599727252,719.21501129307808,719.57628412807946,720.25010577821195,
!     &           721.50993516436063,723.86654051086100,728.26518647629769,736.42299633394953,751.34799183336133,
!     &           777.95722188917193,823.47593814010361,876.10034134591615,935.84511854567438,1002.8896470966265,
!     &           1077.2723727565547,1158.1148589946081,1241.5667790992175,1316.5397902652905,1364.5635331171977,
!     &           1385.1226548601949,1413.9046140350752,1465.7626232167763,1524.0924098449020,1567.7408353719700,
!     &           1584.8527305053142,1587.3708376723982,1587.8967927077658,1588.7533923968615,1590.3067227028648,
!     &           1593.1331676256652,1598.2831166374394,1607.6522157406928,1624.5976694257308,1654.8784620830322,
!     &           1707.9304181965206,1798.5183231978438,1949.8101424451597,2157.3218077666902/)

!         dpe = (/0.20859733251215448,0.22044416870542050,0.23805758824042167,0.26424454167295941,0.30317829492372228,
!     &            0.36106349894703027,0.44712498968589842,0.57507789774188556,0.76531335078162543,1.0481480916489416,
!     &            1.4686558853501106,2.0938506687815250,3.0233661809425669,4.4053340255971900,6.4599905629787742,
!     &            9.5147747266870759,14.056510105357699,20.808987250855296,30.848310064465466,45.774388359017564,
!     &            67.965906391795158,100.95940025406642,150.01285552611540,222.94364255254089,331.37432631977185,
!     &            492.58488857477539,732.26654029583096,1088.6159839561274,1618.4226048850774,2406.1186728147713,
!     &            3577.2348327903010,5318.4052105484952,7907.1102322996994,11755.897103041276,17478.125109680550,
!     &            25985.713204227170,38634.467343784374,57440.147579469718,85399.708278087812,126968.90570473982,
!     &            188772.38657207199,280659.42033206439,417273.52126489556,620386.10211536102,922366.06039229129,
!     &            1371338.2255782168,2038852.7396602808,3031287.5965643590,4506801.4983895328,6700538.7087804656,
!     &            9962102.6576534621,14811270.252659660,22020825.770522781,32739705.094676733/)

         do k = 1, NLAYER
             pl(k) = dpe(k) / log(pe(k+1)/pe(k))
         end do


         CALL calculate_opacities(NLAYER, NSOLP, NIRP,
     &                             mu_0,Tirr, Tint, Tl, Pl, dpe, tau_IRe,tau_Ve, Beta_V, Beta_IR,gravity_SI)

      end subroutine opacity_wrapper

      subroutine calculate_opacities(NLAYER, NSOLP, NIRP,
     &                               mu_0, Tirr, Tint, Tl, Pl, dpe, tau_IRe,tau_Ve,Beta_V, Beta_IR,gravity_SI)
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

        real, dimension(NLAYER) :: dpe, Pl, Tl

        real, dimension(NIRP,NLAYER+1) :: k_IRl, tau_IRe
        real, dimension(NSOLP,NLAYER+1) :: k_Vl, tau_Ve
        real :: grav

        grav = gravity_SI

        !! Parmentier opacity profile parameters - first get Bond albedo
        Teff = ((Tint * Tint * Tint * Tint) + (1.0 / sqrt(3.0)) *
     &          (Tirr * Tirr * Tirr * Tirr)) ** (0.25)

        call Bond_Parmentier(Teff, grav, AB)

        !! Recalculate Teff and then find parameters
        Teff = ((Tint * Tint * Tint * Tint) + (1.0 - AB) * mu_0 *
     &          (Tirr * Tirr * Tirr * Tirr)) ** (0.25)

        ! Log 10 T_eff variables
        l10T = log10(Teff)
        l10T2 = l10T * l10T

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

        ! Calculate tau_lim parameter
        tau_lim = 1.0/(gam_1*gam_2) * sqrt(gam_P/3.0)

        tau_Ve(:,1) = 0.0
        tau_IRe(:,1) = 0.0

        ! SOMETHING IS OFF HERE, WHY CAN"T I CALL THIS FOR NLAYERS????
        ! MALSKY STOP BEING LAZY! AND FIGURE THIS OUT
        do k = 1, NLAYER-1
          call k_Ross_Freedman(Tl(k), pl(k), 0.0, k_IRl(1,k))
          !call k_Ross_Valencia(Tl(k), pl(k), 0.0, k_IRl(1,k))

          k_Vl(:,k) = k_IRl(1,k) * gam_V(:)

          k_IRl(2,k) = k_IRl(1,k) * gam_2
          k_IRl(1,k) = k_IRl(1,k) * gam_1

          tau_Ve(:,k+1)  = (k_Vl(:,k) * dpe(k))  / grav
          tau_IRe(:,k+1) = (k_IRl(:,k) * dpe(k)) / grav
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
        ! Total Rosseland mean opacity - coNLAYERed to m2 kg-1
        k_IR = (10.0**k_lowP + 10.0**k_hiP) / 10.0

        ! Avoid divergence in fit for large values
        k_IR = min(k_IR,1.0e30)
      end subroutine k_Ross_Freedman





      subroutine Bond_Parmentier(Teff0, grav,  AB)
        implicit none

        ! Input:
        ! Teff0 - Atmospheric profile effective temperature [K] with zero albedo
        ! grav - Surface gravity of planet [m s-2]

        ! Output:
        ! AB - Bond albedo

        real, intent(in) :: Teff0, grav
        real, intent(out) :: AB

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
        AB = 10.0**(a + b * log10(Teff0))
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





