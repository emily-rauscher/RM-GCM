      SUBROUTINE SETUPRAD_SIMPLE(Beta_V, Beta_IR, t_pass, incident_starlight_fraction,TAURAY,TAUL,TAUGAS,TAUAER,
     &                         solar_calculation_indexer)
!
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
!
      include 'rcommons.h'

! **********************************************************************
!
!           LOCAL DECLARATIONS
!
! **********************************************************************
      integer :: testing, L, J, solar_calculation_indexer
      REAL G,WVO,AM, incident_starlight_fraction
      real, dimension(NIR,NLAYER) :: tau_IRe
      real, dimension(NSOL,NLAYER) :: tau_Ve
      real, dimension(NWAVE) :: MALSKY_ABSCOEFF(NWAVE)
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL) :: Beta_V
      dimension rup_1(NGROUP)
      dimension rhoi(NRAD), dbnds(NRAD+1)
      dimension zbnds(6), pbnds(6), rn2ds(NRAD,6)
      dimension tauem(5,NWAVE), ssam(5,NWAVE), asmm(5,NWAVE)
      dimension temparr(6,NWAVE)
      dimension pbndsm(6)
      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL, TAUGAS,TAUAER


      real t_pass(NZ)
      integer i1, i2, indorder(5)
      logical all_ok
      integer :: malsky_switch


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
!
      DATA GANGLE  / 0.2123405382, 0.5905331356,0.9114120405/
      DATA GRATIO  / 0.4679139346, 0.3607615730, 0.1713244924/
      DATA GWEIGHT /  0.0698269799, 0.2292411064,0.2009319137 /

      DATA AVG    /6.02252E+23/
      DATA PI     /3.14159265359/
!     ALOS   - LOCSHMIDT'S NUMBER (#/CM**3)
!     AM     - MOLECULAR WEIGHT OF AIR (G/MOL)
!     AVG    - AVAGODROS' NUMBER (#/MOL)
!     G      - GRAVITY (CM/S**2)
!     PI     - PI
!     RGAS   - UNIVERSAL GAS CONSTANT (ERG / MOL K)
!     SCDAY  - NUMBER OF SECONDS IN ONE DAY (S)

      AM= RGAS/R_AIR
      DATA ALOS   / 2.68719E19   /
      DATA RGAS   / 8.31430E+07  /
      DATA SBK    / 5.6697E-8    /
      DATA SCDAY  / 86400.0      /
      DATA EPSILON / ALMOST_ZERO  /

      G = GA*100.
      AM= RGAS/R_AIR

      DO L = 1,NSOLP
          MALSKY_ABSCOEFF(L)=ABSSW
      END DO

      DO L = NSOLP+1,NTOTAL
          MALSKY_ABSCOEFF(L)=ABSLW
      END DO

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
!     PTOP       - PRESSURE AT TOP OF MODEL (DYNES/CM**2)
!     PBOT       - PRESSURE AT BOTTOM OF MODEL (DYNES/CM**2)
!     SFC_WIND   - wind speed at 10 m altitude (m/s)
!     SFC_ALB    - surface albedo when fixed

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
      PTOP         = p_aerad(1)*10.
      PBOT         = p_aerad(NL+1)*10.

      ALBEDO_SFC = ALBSW

      testing = 0
      if (testing .eq. 1) then
          p_aerad(1) = 10.0 ** (LOG10(p_aerad(2)) - (LOG10(p_aerad(3)) - LOG10(p_aerad(2))))

          PRESSMID = PLAYER*10.
          PRESS    = P_aerad*10.

          DO J  = 2,NLAYER
             PBAR(J)  = (p_aerad(J)-p_aerad(J-1))*1e-5
             DPG(J) = (PRESS(J)-PRESS(J-1)) / G
          END DO

          PBAR(1)  = 10.0 ** (LOG10(PBAR(2)) - (LOG10(PBAR(3)) - LOG10(PBAR(2))))
          DPG(1)   = 10.0 ** (LOG10(DPG(2))  - (LOG10(DPG(3))  - LOG10(DPG(2))))

           k = 1
           DO J = 1, NDBL-4, 2
               PBARsub(J)   = PBAR(k)
               PBARsub(J+1) = PBAR(k) + ABS(PBAR(k+2) - PBAR(k+1)) / 2.0

               DPGsub(J)   = DPG(k)
               DPGsub(J+1) = DPG(k) + ABS(DPG(k+2) - DPG(k+1)) / 2.0
               k = k + 1
           END DO

           PBARsub(NDBL-3) = PBARsub(NDBL-4) + ABS(PBARsub(NDBL-5) - PBARsub(NDBL-6))
           PBARsub(NDBL-2) = PBARsub(NDBL-3) + ABS(PBARsub(NDBL-4) - PBARsub(NDBL-5))
           PBARsub(NDBL-1) = PBARsub(NDBL-2) + ABS(PBARsub(NDBL-3) - PBARsub(NDBL-4))
           PBARsub(NDBL)   = PBARsub(NDBL-1) + ABS(PBARsub(NDBL-2) - PBARsub(NDBL-3))

           DPGsub(NDBL-3) = DPGsub(NDBL-4) + ABS(DPGsub(NDBL-5) - DPGsub(NDBL-6))
           DPGsub(NDBL-2) = DPGsub(NDBL-3) + ABS(DPGsub(NDBL-4) - DPGsub(NDBL-5))
           DPGsub(NDBL-1) = DPGsub(NDBL-2) + ABS(DPGsub(NDBL-3) - DPGsub(NDBL-4))
           DPGsub(NDBL)   = DPGsub(NDBL-1) + ABS(DPGsub(NDBL-2) - DPGsub(NDBL-3))

           PBARsub(1)  = 10.0 ** (LOG10(PBARsub(2)) - (LOG10(PBARsub(3)) - LOG10(PBARsub(2))))
           DPGsub(1)   = 10.0 ** (LOG10(DPGsub(2))  - (LOG10(DPGsub(3))  - LOG10(DPGsub(2))))
      else
    !     Get atmospheric pressure profile from interface common block
    !     [ dyne / cm^2 ]

          do k = 1,NZ
            press(k)=p_aerad(k)*10.
          enddo

    !     The layer thickness in pressure should be the difference between
    !     the top and bottom edge of the layer.  Eg. So for layer 3, the layer thickness
    !     is the pressure at the boundaries 4 minus at the boundary 3.
    !     These boundary pressures are given in p_full (NZ+1 of them).
    !     The mass of these layers is this pressure thickness divide by G.
    !
    !     (NOTE - THE TOP LAYER IS FROM PTOP TO 0, SO AVERAGE = PTOP/2)
    !     P_aerad=PRESSURE AT EDGES OF LAYERS (PASCALS)
    !     PRESS - PRESSURE AT EDGE OF LAYER (dyne/cm^2)
    !     DPG   - MASS OF LAYER (G / CM**2)!    D PBAR
    !     PBARS - THICKNESS OF LAYER IN PRESSURE (BARS)
    !     PRESSMID- PRESSURE AT CENTER OF LAYER (dyne/cm^2

          PRESSMID=PLAYER*10.
          PRESS=P_aerad*10.
          PBAR(1)  = P_AERAD(1)*1e-5
          DPG(1)= PRESS(1)/G

          DO J  = 2,NLAYER
             PBAR(J)  = (p_aerad(J)-p_aerad(J-1))*1e-5
             DPG(J) = (PRESS(J)-PRESS(J-1)) / G
          END DO
                     K  =  1
          DO J  = 2, NDBL,2
                     L  =  J
             PBARsub(L) =  (PRESSMID(K) - PRESS(K))*1e-6
             DPGsub (L) =  (PRESSMID(K) - PRESS(K)) / G
                     K  =  K+1
                     L  =  L+1
             PBARsub(L) =  (PRESS(K) - PRESSMID(K-1))*1e-6
             DPGsub (L) =  (PRESS(K) - PRESSMID(K-1)) / G
          END DO
             PBARsub(1) = PBAR(1)
             DPGsub(1)  = DPG(1)
      end if

      TAURAY(:,:) = 0.0
      TAUAER(:,:) = 0.0
      TAUGAS(:,:) = 0.0

      TAUCLD(:,:)  = 0.

      WCLD(:,:)    = 0.
      WOL(:,:)     = 0.
      GOL(:,:)     = 0.
      GCLD(:,:)    = 0.


      malsky_switch = 0
      IF (malsky_switch .gt. 0) THEN
        CALL opacity_wrapper(t_pass, tau_IRe, tau_Ve, Beta_V, Beta_IR, GA, incident_starlight_fraction)

        DO L = solar_calculation_indexer,NSOLP
          tau_Ve(L,NLAYER) = 10.0**(LOG10(tau_Ve(L,NLAYER-1))+(LOG10(tau_Ve(L,NLAYER-1)) - LOG10(tau_Ve(L,NLAYER-2))))
        END DO

        DO L = NSOLP+1, NTOTAL
          tau_IRe(L-NSOLP,NLAYER) = 10.0 ** (LOG10(tau_IRe(L-NSOLP,NLAYER-1))+
     &            (LOG10(tau_IRe(L-NSOLP,NLAYER-1))-LOG10(tau_IRe(L-NSOLP,NLAYER-2))))
        END DO

        DO L = solar_calculation_indexer,NSOLP
            DO J = 1,NLAYER
                TAUGAS(L,J) = tau_Ve(L,J)
            END DO
        END DO


        DO L = NSOLP+1, NTOTAL
            k  =  1
            DO  J = 1,NDBL,2
                TAUGAS(L, J)   = tau_IRe(L - NSOLP, k)
                TAUGAS(L, J+1) = tau_IRe(L - NSOLP, k)+ ABS(tau_IRe(L - NSOLP,k) - tau_IRe(L - NSOLP,k+1)) / 2.0
                k = k + 1
            END DO
        END DO


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

          DO L = 1,NSOLP
              DO J     =   1,NLAYER
                  PM          =   DPG(J)
                  TAUGAS(L,J) = MALSKY_ABSCOEFF(L)*PM
              END DO
          END DO

          DO L  = NSOLP+1,NTOTAL
             DO J     =   1,NDBL
                 PM          =   DPGsub(J)
                 TAUGAS(L,J)=MALSKY_ABSCOEFF(L)*PM
             END DO
          END DO
      END IF

      FNET(:,:)   = 0.0
      TMI(:,:)    = 0.0
      DIRECT(:,:) = 0.0



      DO  L = NSOLP+1,NTOTAL
          TAUCONST(L)=MALSKY_ABSCOEFF(L)/GA/100.
      ENDDO

!@@@@@@@@@@@@@@RAYLEIGH SCATTERING CONDITIONAL@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     WAVE MUST BE IN MICRONS
!     CALCULATE RAYLEIGH OPTICAL DEPTH PARAMETERS.

!     RAYLEIGH SCATTERING CONDITIONAL:
! NOTE: The conversion factor KMAMGperBAR is computed in insimprad
! KMAMGperBAR  = AVG*1.E5/(LO*GA*MWTOT)
!it is derived by the relation:
!km-amagats= Pressure * Avogadro's# * mole fraction/  (Loshchmidt's# *
!gravity* molecular weight)
!where pressure is in Pa, Avogadro's=6.0221367*10^23 mol^-1
!Lo=2.6867630*10^25 m^-3, gravity= 24.40 m s^-2
!molec. weight for mole fraction .86 H2 and .136 He (von Zahn)
! 2.27 *10^-3 kg/mole

      IF (RAYSCAT) THEN
        DO  L = 1,NTOTAL
          RAYPERBAR(L) = RAYPERBARCONS
        END DO

        DO 330 J          =   1,NLAYER
          DO 335 L         =   1,NTOTAL
            if( L .LE. NSOLP )then
              TAURAY(L,J) = RAYPERBAR(L)*PBAR(J) !PER BAR X LAYER THICKNESS IN BAR
            else
              TAURAY(L,J)= 0.0
            endif
335     CONTINUE
330   CONTINUE
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
!
!
!     Set <iblackbody_above> = 1 to include a source of radiation
!     at the top of the radiative transfer model domain
!
      iblackbody_above = ir_above_aerad
      t_above = tabove_aerad
!
!     Set <ibeyond_spectrum> = 1 to include blackbody radiation at
!     wavelengths longer than WAVE(NWAVE+1) in PLANK(NWAVE+1-NSOL) and
!     at wavelengths shorter than WAVE(NSOL+1) in PLANK(1)
!
      ibeyond_spectrum = 0


      RETURN
      END
