      SUBROUTINE SETUPRAD_SIMPLE(Beta_V, Beta_IR, t, pr, P_PASS,
     &                           incident_starlight_fraction,TAURAY,TAUL,TAUGAS,TAUAER,
     &                           solar_calculation_indexer, DPG,
     &             ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &             fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &             pbar, dpgsub, pbarsub,
     &           LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,
     &           EMISIR, EPSILON, HEATI, HEATS, HEAT, SOLNET, TPI, SQ3, SBK, AM, AVG, ALOS,
     &  SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
     &  GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
     &  WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
     &  uG0,uTAUL,W0,uW0,uopd,U1S,
     &  U1I,TOON_AK,B1,B2,EE1,EM1,
     &  EM2,EL1,EL2,GAMI,AF,
     &  BF,EF,SFCS,B3,CK1,CK2,
     &  CP,CPB,CM,CMB,DIRECT,EE3,
     &  EL3,FNET,TMI,AS,DF,
     &  DS,XK,DIREC,DIRECTU,DINTENT,
     &  UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
     &  tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
     &  fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &  qrad,alb_tomi,alb_toai, num_layers, SLOPE,
     &  dpe, Pl, Tl, pe,
     &  k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
     &  Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val, tau_IRe, tau_Ve, k_IRl, k_Vl)

!     *********************************************************
!     *  Purpose            :  Defines all constants, and     *
!     *                        calculates pressure averaged   *
!     *                        absorption coefficients.       *
!     *                                                       *
!     *       THIS IS A SIMPLER VERSION OF SETUPRAD           *
!     * If you are attempting to add complexity to the model, *
!     * Such as adding more wavelength bins or mie calculations
!     * then consult that program first.                      *
!     * *******************************************************
      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5),NPROB(5), SOL(5),RAYPERBAR(5),WEIGHT(5)
      REAL GOL(5,2*NL+2), WOL(5,2*NL+2), WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2), PTEMP(5,2*NL+2)
      REAL uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2), uopd(5,2*NL+2),  U1S( 5)
      REAL U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(  5,2*NL+2), EE1( 5,2*NL+2), EM1(5,2*NL+2)
      REAL EM2(5,2*NL+2), EL1( 5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2), EE3(5,2*NL+2)
      REAL EL3(5,2*NL+2), FNET(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4), DF(5,4*NL+4)
      REAL DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2), DINTENT(5,3,2*NL+2)
      REAL UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(2),fird(2),fsLu(3), fsLd(3),fsLn(3),alb_toa(3), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai, SLOPE(5,2*NL+2)

      REAL tau_IRe(2,NL+1), tau_Ve(3,NL+1)
      real, dimension(NL+1) :: dpe, Pl, Tl, pe
      real :: k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met
      real :: Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val

      real, dimension(2, NL+1) :: k_IRl
      real, dimension(2, 2*NL+2) :: k_irl_doubled

      real, dimension(3, NL+1) :: k_Vl

      ! New variables for calculating the IR absorbtion coefficient as a power law
      real, dimension(NLAYER) :: IR_ABS_COEFFICIENT

! **********************************************************************
!
!           LOCAL DECLARATIONS
!
! **********************************************************************
      integer :: L, J, K, solar_calculation_indexer, I
      REAL G,WVO, incident_starlight_fraction
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL) :: Beta_V
      dimension rup_1(NGROUP)
      dimension rhoi(NRAD), dbnds(NRAD+1)
      dimension zbnds(6), pbnds(6), rn2ds(NRAD,6)
      dimension tauem(5,NWAVE), ssam(5,NWAVE), asmm(5,NWAVE)
      dimension temparr(6,NWAVE)
      dimension pbndsm(6)
      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL, TAUGAS,TAUAER
      real dpg(nl+1), pbar(nl+1)
      real dpgsub(2*nl+2), pbarsub(2*nl+2)
      real t(NLAYER), pr(NLAYER)
      integer i1, i2, indorder(5)
      logical all_ok
      integer ifsetup
      real ibinm
      real rfluxes_aerad(2,2,2)
      real psol_aerad
      real heati_aerad(NL+1)
      real heats_aerad(NL+1)
      real fsl_up_aerad(NL+1)
      real fsl_dn_aerad(NL+1)
      real fir_up_aerad(NL+1)
      real fir_dn_aerad(NL+1)
      real fir_net_aerad(NL+1)
      real fsl_net_aerad(NL+1)

      ! For getting a doubled grid for the IR channels
      REAL :: LOG_START, LOG_END, LOG_STEP
      REAL, DIMENSION(NLAYER) :: P_PASS
      REAL, DIMENSION(2*NLAYER) :: P_PASS_SUB

! ******************************************
!            DEFINE CONSTANTS
! *****************************************
!
!     UNITS ARE (CM**2)/GM
!
!     ***********************
!     *  DATA FOR INFRARED  *
!     ***********************
!
!     GAUSS ANGLES AND GAUSS WEIGHTS FOR GAUSSIAN INTEGRATION
!     MOMENTS (USE FIRST MOMENT VALUES) N=3

      DATA PI     /3.14159265359/
!     ALOS   - LOCSHMIDT'S NUMBER (#/CM**3)
!     AM     - MOLECULAR WEIGHT OF AIR (G/MOL)
!     AVG    - AVAGODROS' NUMBER (#/MOL)
!     G      - GRAVITY (CM/S**2)
!     PI     - PI
!     RGAS   - UNIVERSAL GAS CONSTANT (ERG / MOL K)
!     SCDAY  - NUMBER OF SECONDS IN ONE DAY (S)

      AM= RGAS/R_AIR
      G = GA*100.
      AM= RGAS/R_AIR

      SQ3     =   SQRT(3.)
      JDBLE   =   2*NLAYER
      JDBLEDBLE = 2*JDBLE
      JN      =   JDBLE-1
      JN2     =   2*JDBLE-1
      TPI     =   2.*PI
      CPCON   =   GASCON/AKAP/1000.  ! Cp in J/gm/K
      FDEGDAY =   1.0E-4*G*SCDAY/CPCON

!     Get scalars from interface common block:
!     ISL        - do solar calculations when = 1
!     IR         - do infrared calculations when = 1
!     IRS        - do infrared scattering when = 1
!     FLXLIMDIF  - do flux limited diffusion correction = 1 ~MTR
!     UO         - SOLAR ZENITH ANGLE
!     EMISIR     - SURFACE IR EMISSIVITY
!     SFC_WIND   - wind speed at 10 m altitude (m/s)

      ISL          = 0
      IR           = 0
      IRS          = 0
      IF(DOSWRAD) THEN
        ISL = 1
      ENDIF
      IF(DOLWRAD) THEN
        IR  = 1
      ENDIF
      IF(LWSCAT) THEN
        IRS  = 1
      ENDIF

      ! SET WAVELENGTH LIMITS LLA AND LLS BASED ON VALUES OF ISL AND IR
      LLA = NTOTAL
      LLS = 1

      IF(ISL .EQ. 0) THEN
          LLS =  NSOLP+1
      ENDIF

      IF(IR .EQ. 0) THEN
          LLA =  NSOLP
      ENDIF

      EMISIR       = SURFEMIS

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!    SET UP THE LAYER TEMPS     !!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO J = 2, NVERT
          TT(J) = T(J-1) * (P_PASS(J)/pr(J-1)) ** (log(T(J)/T(J-1))/log(pr(J)/pr(J-1)))
      END DO

      TT(1)=((T(1)-TT(2))/log(pr(1)/P_PASS(2)))*log(P_PASS(1)/pr(1))+T(1)

      TT(NLAYER) = T(NVERT) * ((P_PASS(NLAYER)*10)/(pr(NVERT)*10.0)) **
     &             (log(T(NVERT)/T(NVERT-1))/log((pr(NVERT)*10.0)/(pr(NVERT-1)*10.0)))

    !     The layer thickness in pressure should be the difference between
    !     the top and bottom edge of the layer.  Eg. So for layer 3, the layer thickness
    !     is the pressure at the boundaries 4 minus at the boundary 3.
    !     These boundary pressures are given in p_full (NZ+1 of them).
    !     The mass of these layers is this pressure thickness divide by G.
    !
    !     (NOTE - THE TOP LAYER IS FROM PTOP TO 0, SO AVERAGE = PTOP/2)
    !     P_PASS=PRESSURE AT EDGES OF LAYERS (PASCALS)
    !     PRESS - PRESSURE AT EDGE OF LAYER (dyne/cm^2)
    !     DPG   - MASS OF LAYER (G / CM**2)!    D PBAR
    !     PBARS - THICKNESS OF LAYER IN PRESSURE (BARS)
      P_PASS(1) = 10.0 ** (LOG10(P_PASS(2)) - (LOG10(P_PASS(3)) - LOG10(P_PASS(2))))

      DO J  = 2,NLAYER
          PBAR(J)  = (P_PASS(J)-P_PASS(J-1))*1e-5
          DPG(J)   = ((P_PASS(J) * 10.0)-(P_PASS(J-1)*10.0)) / G
      END DO

      PBAR(1)  = 10.0 ** (LOG10(PBAR(2)) - (LOG10(PBAR(3)) - LOG10(PBAR(2))))
      DPG(1)   = 10.0 ** (LOG10(DPG(2))  - (LOG10(DPG(3))  - LOG10(DPG(2))))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! New interpolation to smooth the layers a bit !!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      LOG_START = LOG(P_PASS(1))
      LOG_END = LOG(P_PASS(NLAYER))
      LOG_STEP = (LOG_END - LOG_START) / (2*NLAYER - 1)

      ! Creating the new grid
      DO I = 1, 2*NLAYER
          P_PASS_SUB(I) = EXP(LOG_START + LOG_STEP * (I - 1))
      END DO

      ! Calculating PBARSUB and DPGSUB for indices 2 to 2*NLAYER
      DO J = 2, 2*NLAYER
          PBARSUB(J) = (P_PASS_SUB(J) - P_PASS_SUB(J-1)) * 1e-5
          DPGSUB(J)  = ((P_PASS_SUB(J) * 10.0) - (P_PASS_SUB(J-1) * 10.0)) / G
      END DO

      ! Extrapolating the first values using a logarithmic trend
      PBARSUB(1) = 10.0 ** (LOG10(PBARSUB(2)) - (LOG10(PBARSUB(3)) - LOG10(PBARSUB(2))))
      DPGSUB(1)  = 10.0 ** (LOG10(DPGSUB(2))  - (LOG10(DPGSUB(3))  - LOG10(DPGSUB(2))))




      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! Keeping Michael Roman's Old Regridding Just !!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (.FALSE.) THEN
          PBAR(1)                = p_pass(1) * 1e-5
          DPG(1)                 = (p_pass(1) * 10.0)/G

          DO J  = 2,NLAYER
              PBAR(J) = (p_pass(J)-p_pass(J-1))*1e-5
              DPG(J)  = ((p_pass(J) * 10.0)-(p_pass(J-1)*10.0)) / G
          END DO

          K = 1
          DO J  = 2, NDBL-1,2
             L  =  J
             PBARsub(L) =  ((pr(K)*10.0) - (p_pass(K) * 10.0))*1e-6
             DPGsub (L) =  ((pr(K)*10.0) - (p_pass(K) * 10.0)) / G
             K = K + 1
             L = L + 1
             PBARsub(L) =  ((p_pass(K) * 10.0) - (pr(K-1)*10.0))*1e-6
             DPGsub (L) =  ((p_pass(K) * 10.0) - (pr(K-1)*10.0)) / G
          END DO

          PBARsub(NDBL) = PBARsub(NDBL-1) + ABS(PBARsub(NDBL-2) - PBARsub(NDBL-3))
          DPGsub(NDBL)  = DPGsub(NDBL-1)  + ABS(DPGsub(NDBL-2)  - DPGsub(NDBL-3))

          PBARsub(1) = PBAR(1)
          DPGsub(1)  = DPG(1)
      END IF

      TAURAY(:,:) = 0.0
      TAUAER(:,:) = 0.0
      TAUGAS(:,:) = 0.0
      TAUL(:,:)   = 0.0

      WOL(:,:)    = 0.0
      GOL(:,:)    = 0.0

      IF (picket_fence_optical_depths) THEN
          IF (NIRP .EQ. 1) THEN
              WRITE(*,*) 'Stopping! Running picket fence gas optics with the wrong number of channels'
              STOP
          ENDIF
      ENDIF

      IF (picket_fence_optical_depths) THEN
        CALL opacity_wrapper(t, P_PASS, tau_IRe, tau_Ve, Beta_V, Beta_IR, GA, incident_starlight_fraction,
     &           LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,
     &           EMISIR, EPSILON, HEATI, HEATS, HEAT, SOLNET, TPI, SQ3, SBK, AM, AVG, ALOS,
     &  SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
     &  GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
     &  WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
     &  uG0,uTAUL,W0,uW0,uopd,U1S,
     &  U1I,TOON_AK,B1,B2,EE1,EM1,
     &  EM2,EL1,EL2,GAMI,AF,
     &  BF,EF,SFCS,B3,CK1,CK2,
     &  CP,CPB,CM,CMB,DIRECT,EE3,
     &  EL3,FNET,TMI,AS,DF,
     &  DS,XK,DIREC,DIRECTU,DINTENT,
     &  UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
     &  tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
     &  fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &  qrad,alb_tomi,alb_toai, num_layers,
     &  dpe, Pl, Tl, pe,
     &  k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
     &  Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val, k_IRl, k_Vl)

        DO L = solar_calculation_indexer,NSOL
          DO J     =   1,NLAYER
            TAUGAS(L,J) = k_VL(L,J)*10. * DPG(J) ! 10 converts from m^2/kg to cm^2/g
          END DO
        END DO
        k_irl_doubled = 0.0
        ! smooth out the IR opacities to twice the resolution (linear interpolation)
        DO L = NSOL+1,NTOTAL
          k = 1
          DO J     =   1,NDBL, 2
            k_irl_doubled(L-NSOL,J) = k_irl(L-NSOL, k)
            k_irl_doubled(L-NSOL,J+1) = k_irl(L-NSOL, k)+ ABS(k_irl(L-NSOL,k) - k_irl(L-NSOL,k+1)) / 2.0
            k = k + 1
          END DO
        END DO
        DO L  = NSOL+1,NTOTAL
          DO J     =   1,NDBL
            TAUGAS(L,J) = k_irl_doubled(L-NSOL,J)*10. * DPGSUB(J) ! 10 converts from m^2/kg to cm^2/g
          END DO
        END DO
    ! Removed in v5.2
    !     DO L = solar_calculation_indexer,NSOLP
    !       tau_Ve(L,NLAYER) = 10.0**(LOG10(tau_Ve(L,NLAYER-1))+(LOG10(tau_Ve(L,NLAYER-1)) - LOG10(tau_Ve(L,NLAYER-2))))
    !     END DO

    !     DO L = NSOLP+1, NTOTAL
    !       tau_IRe(L-NSOLP,NLAYER) = 10.0 ** (LOG10(tau_IRe(L-NSOLP,NLAYER-1))+
    !  &            (LOG10(tau_IRe(L-NSOLP,NLAYER-1))-LOG10(tau_IRe(L-NSOLP,NLAYER-2))))
    !     END DO

    !     DO L = solar_calculation_indexer,NSOLP
    !         DO J = 1,NLAYER
    !             TAUGAS(L,J) = tau_Ve(L,J)
    !         END DO
    !     END DO

    !     DO L = NSOLP+1, NTOTAL
    !         k  =  1
    !         DO  J = 1,NDBL,2
    !             TAUGAS(L, J)   = tau_IRe(L - NSOLP, k)
    !             TAUGAS(L, J+1) = tau_IRe(L - NSOLP, k)+ ABS(tau_IRe(L - NSOLP,k) - tau_IRe(L - NSOLP,k+1)) / 2.0
    !             k = k + 1
    !         END DO
    !     END DO
      ELSE
          if (NSOLP .gt. 1) then
              Beta_V(1) = 1.0
              Beta_V(2) = 0.0
              Beta_V(3) = 0.0

              Beta_IR(1) = 1.0
              Beta_IR(2) = 0.0
          else
              Beta_V(1)  = 1.0
              Beta_IR(1) = 1.0
          end if

          ! Calculate the opacity power law
          ! Split it up into if statements for better efficiency
          ! If the opacity power law is zero, resume using double gray


          IF (OPACIR_POWERLAW.eq.0) THEN
            DO J  = 1,(2*nl+2)
              IR_ABS_COEFFICIENT(J) = ABSLW
            END DO
          ELSE IF (OPACIR_POWERLAW.eq.1) THEN
            DO J  = 1,(2*nl+2)
              IR_ABS_COEFFICIENT(J) = ABSLW
     &        * MAX(1e-6, (P_PASS_SUB(J) / OPACIR_REFPRES))
            END DO
          ELSE IF (OPACIR_POWERLAW.eq.2) THEN
            DO J  = 1,(2*nl+2)
              IR_ABS_COEFFICIENT(J) = ABSLW
     &        * MAX(1e-6, (P_PASS_SUB(J) / OPACIR_REFPRES))
     &        * MAX(1e-6, (P_PASS_SUB(J) / OPACIR_REFPRES))
            END DO
          ELSE IF (OPACIR_POWERLAW.eq.3) THEN
            DO J  = 1,(2*nl+2)
              IR_ABS_COEFFICIENT(J) = ABSLW
     &        * MAX(1e-6, (P_PASS_SUB(J) / OPACIR_REFPRES))
     &        * MAX(1e-6, (P_PASS_SUB(J) / OPACIR_REFPRES))
     &        * MAX(1e-6, (P_PASS_SUB(J) / OPACIR_REFPRES))
            END DO
          ELSE
            DO J  = 1,(2*nl+2)
              IR_ABS_COEFFICIENT(J) = ABSLW * MAX(1e-6, ((P_PASS_SUB(J) / OPACIR_REFPRES) ** OPACIR_POWERLAW))
            END DO
          END IF

          ! Set the tau gas equal to the absorbtion coefficient times dpg
          DO L = solar_calculation_indexer,NSOLP
              DO J     =   1,NLAYER
                  TAUGAS(L,J) = ABSSW * DPG(J)
              END DO
          END DO

          DO L  = NSOLP+1,NTOTAL
             DO J     =   1,NDBL
                 TAUGAS(L,J)=IR_ABS_COEFFICIENT(J)*DPGsub(J)
             END DO
          END DO
      END IF

      FNET(:,:)   = 0.0
      TMI(:,:)    = 0.0
      DIRECT(:,:) = 0.0

!     @@@@@@@@@    RAYLEIGH SCATTERING CONDITIONAL    @@@@@@@@@@@@@@@
!     WAVE MUST BE IN MICRONS
!     CALCULATE RAYLEIGH OPTICAL DEPTH PARAMETERS.

      IF (RAYSCAT) THEN
        DO J = 1,NLAYER
          ! Calculate the rayleigh scattering
          DO L = 1,NTOTAL
            if( L .LE. NSOLP )then
              TAURAY(L,J) = RAYPERBARCONS(L) * PBAR(J)
            else
              TAURAY(L,J)= 0.0
            endif
          END DO
        END DO
      ELSE
        DO 320 J     = 1,NLAYER
          DO 325 L    = 1,NTOTAL
            TAURAY(L,J)= 0.0
325       CONTINUE
320     CONTINUE
      ENDIF

      DO 360 L   =   1,NSOLP
        SOL(L)  = PSOL_aerad
 360  CONTINUE


! *********************************************************************
!
!     COMPUTE PLANCK FUNCTION TABLE. WAVE IS IN UNITS OF MICRONS.
!
! **********************************************************************
!     Set <iblackbody_above> = 1 to include a source of radiation
!     at the top of the radiative transfer model domain

      iblackbody_above = 0.0

!     Set <ibeyond_spectrum> = 1 to include blackbody radiation at
!     wavelengths longer than WAVE(NWAVE+1) in PLANK(NWAVE+1-NSOL) and
!     at wavelengths shorter than WAVE(NSOL+1) in PLANK(1)

      ibeyond_spectrum = 0

      RETURN
      END
