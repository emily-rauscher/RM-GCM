      SUBROUTINE SETUPRAD_SIMPLE(Beta_V, Beta_IR, t_pass)
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

      REAL p(NZ),ABSCOEFF(NWAVE),G,WVO,AM

      real, dimension(2,NLAYER) :: tau_IRe
      real, dimension(3,NLAYER) :: tau_Ve
      real, dimension(2) :: Beta_IR
      real, dimension(3) :: Beta_V
      dimension rup_1(NGROUP)
      dimension rhoi(NRAD), dbnds(NRAD+1)
      dimension zbnds(6), pbnds(6), rn2ds(NRAD,6)
      dimension tauem(5,NWAVE), ssam(5,NWAVE), asmm(5,NWAVE)
      dimension temparr(6,NWAVE)
      dimension pbndsm(6)

      real t_pass(NZ)


      integer i1, i2, indorder(5)
      logical all_ok
      DATA AVG    / 6.02252E+23  /
      DATA PI     /3.14159265359/

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
      DATA GANGLE  /  0.2123405382, 0.5905331356,     
     &               0.9114120405                       /
      DATA GRATIO/0.4679139346, 0.3607615730, 0.1713244924/
      DATA GWEIGHT /  0.0698269799, 0.2292411064,     
     &                 0.2009319137                       /
!
!
!     ALOS   - LOCSHMIDT'S NUMBER (#/CM**3)
!     AM     - MOLECULAR WEIGHT OF AIR (G/MOL)
!     AVG    - AVAGODROS' NUMBER (#/MOL)
!     G      - GRAVITY (CM/S**2)
!     PI     - PI
!     RGAS   - UNIVERSAL GAS CONSTANT (ERG / MOL K)
!     SCDAY  - NUMBER OF SECONDS IN ONE DAY (S)

      AM= RGAS/R_AIR
      DATA ALOS   / 2.68719E19   /
      G= GA*100.
      DATA RGAS   / 8.31430E+07  /
      DATA SBK    / 5.6697E-8    /
      DATA SCDAY  / 86400.0      /


      DATA EPSILON / ALMOST_ZERO  /

      AM= RGAS/R_AIR
      ABSCOEFF(1)=ABSSW
      ABSCOEFF(2)=ABSLW
      ABSCOEFF(3)=ABSLW
      ABSCOEFF(4)=ABSLW
      ABSCOEFF(5)=ABSLW

      SQ3     =   SQRT(3.)
      JDBLE   =   2*NLAYER
      JDBLEDBLE = 2*JDBLE
      JN      =   JDBLE-1
      JN2     =   2*JDBLE-1
      TPI     =   2.*PI
      CPCON   =   GASCON/AKAP/1000.  ! Cp in J/gm/K 
      FDEGDAY =   1.0E-4*G*SCDAY/CPCON


!     Get scalars from interface common block:
!
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
!

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
      LLA                   =  NTOTAL
      LLS                   =  1

      IF(ISL .EQ. 0) THEN
          LLS =  NSOLP+1
      ENDIF

      IF(IR .EQ. 0) THEN
          LLA =  NSOLP
      ENDIF

      EMISIR       = SURFEMIS
      PTOP         =p_aerad(1)*10.
      PBOT         =p_aerad(NL+1)*10.

      ALBEDO_SFC = ALBSW !SFC_ALB MTR

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

!        Here we are saying that the top layer has a thickness that
!        extends from the top of the atmosphere (pressure = 0) down to
!        the interface ABOVE the top sigma levels. The atmosphere
!        above the upper most flux boundary (Press(1)) is not explicitly
!        heated or emitting; it simply attenuates the incoming and
!        outgoing radiation by a very small amount.
!        In practice, the small amount of mass is non very signicant,
!        and the absorption coefficient at such rarified pressures and
!        temperatures is probably not very valid. 
!        DPG and PBAR include this mass above the model in index.  Since
!        there are N vertical levels and N+1 layers (i.e. boundaries,
!        interfaces), the arrays of DPG and PBAR have N+1 too account
!        for this cap above the model.  The attenuation through it must
!        still be calculated using these values.

!     CALCULATE CLOUD AND AEROSOL OPACITY.

      DO J             = 1,NDBL
          DO L           = 1,NTOTAL
               TAUAER(L,J)  = 0.0
               TAUCLD(L,J)  = 0.0
               WCLD(L,J)    = 0.0
               WOL(L,J)     = 0.0
               GOL(L,J)     = 0.0
               GCLD(L,J)    = 0.0
               TAURAY(L,J)  = 0.0
          END DO
      END DO


!     THIS IS DOUBLE-GRAY SPECIFIC. NOT YET GENERALIZED
!     SHORTWAVE:
      DO L = LLS,NSOLP
          DO J = 1,NLAYER
              PM          =   DPG(J)
              TAUGAS(L,J) = ABSCOEFF(L)*PM
              FNET(L,J)   = 0.0
              TMI(L,J)    = 0.0
              DIRECT(L,J) = 0.0
          END DO
      END DO


      malsky_switch = 0
      if (malsky_switch .gt. 0) then
          CALL opacity_wrapper(t_pass, tau_IRe, tau_Ve, Beta_V, Beta_IR, GA)

          !IF ((ANY(TT .ge. 1d4)) .or. (ANY(TT .lt. 0)) .or. (ANY(tau_IRe .gt. 1d10)) .or. (ANY(tau_Ve .gt. 1d10))) then
          !    DO J = 1, NLAYER
          !        write(*,*) TT(J), tau_IRe(1, J), tau_IRe(2, J)
          !    END DO
          !    STOP
          !END IF

          DO L = 1,NSOLP
              tau_Ve(L, 1)  = ABS(tau_Ve(L,3) - tau_Ve(L,2))
          END DO

          DO L = NSOLP+1, NTOTAL
              tau_IRe(L - NSOLP,1) = ABS(tau_IRe(L - NSOLP,3) - tau_IRe(L - NSOLP,2))
          END DO

          DO J = 1,NLAYER
              DO L = 1,NSOLP
                  TAUGAS(L,J) = tau_Ve(L,J)
              END DO
          END DO

          DO L = NSOLP+1, NTOTAL
              DO J = 1,NLAYER - 1
                  TAUGAS(L,2*J-1) = tau_IRe(L - NSOLP,J) / 2
                  TAUGAS(L,2*J)   = (tau_IRe(L - NSOLP,J) + tau_IRe(L - NSOLP,J+1)) / 4
              END DO

              TAUGAS(L,2*NLAYER-1) = tau_IRe(L-NSOLP, NLAYER)/2
              TAUGAS(L,2*NLAYER)   = tau_IRe(L-NSOLP,NLAYER)/2 + ABS(tau_IRe(L-NSOLP,NLAYER)-tau_IRe(L-NSOLP,NLAYER-1))/2
          END DO


          DO J = 1, NLAYER
              DO L = 1, NSOLP
                  FNET(L,J)   = 0.0
                  TMI(L,J)    = 0.0
                  DIRECT(L,J) = 0.0
              END DO
          END DO

          DO J = 1, NDBL
              DO L = NSOLP+1, NTOTAL
                  FNET(L,J)   = 0.0
                  TMI(L,J)    = 0.0
                  DIRECT(L,J) = 0.0
              END DO
          END DO
      end if






!     THIS IS DOUBLE-GRAY SPECIFIC. NOT YET GENERALIZED
!     SHORTWAVE:
      DO L = LLS,NSOLP
          DO J     =   1,NLAYER
              PM          =   DPG(J)
              TAUGAS(L,J) = ABSCOEFF(L)*PM
              FNET(L,J)   = 0.0
              TMI(L,J)    = 0.0
              DIRECT(L,J) = 0.0
          END DO
      END DO


      DO L  = NSOLP+1,NTOTAL
          DO J     =   1,NDBL
         PM          =   DPGsub(J)
!         write(*,*)'PM',PM
         IF (OPACIR_POWERLAW.eq.0) THEN !MTR Modif to avoid exponent
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
                     ELSE IF (OPACIR_POWERLAW.eq.1) THEN
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
                     ELSE IF (OPACIR_POWERLAW.eq.2) THEN
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
                     ELSE IF (OPACIR_POWERLAW.eq.3) THEN
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
                     ELSE
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     & *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES)**OPACIR_POWERLAW)
          ENDIF
         FNET(L,J)    =  0.0
         TMI(L,J)     =  0.0
         DIRECT(L,J)  =  0.0
          END DO
      END DO


      DO  L           =   NSOLP+1,NTOTAL
          TAUCONST(L)=ABSCOEFF(L)/GA/100.
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
               TAURAY(L,J) = 0.0
             endif
335             CONTINUE
330            CONTINUE
      ELSE
         DO 320 J     = 1,NLAYER
          DO 325 L    = 1,NTOTAL
           TAURAY(L,J)= 0.0
325       CONTINUE
320      CONTINUE
      ENDIF

       DO 360 L   =   1,NSOLP
          SOL(L)  = PSOL_aerad !  SOLC_IN  !SOLFX(NPROB(L)) * WEIGHT(L)
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
