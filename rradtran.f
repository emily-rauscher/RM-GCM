      SUBROUTINE RADTRAN(Beta_V,Beta_IR, incident_starlight_fraction,TAURAY,TAUL,TAUGAS,TAUAER,
     &                   solar_calculation_indexer, DPG, pr, t, p_pass,
     &             ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &             fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &             pbar, dpgsub, pbarsub,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
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
     &  qrad,alb_tomi,alb_toai, num_layers, SLOPE)

!
!     **************************************************************
!     Purpose:    Driver routine for radiative transfer model.
!
!     Input:      Temperature, vapor, and aerosol profiles and solar
!                 zenith angle are taken from interface common block.
!
!     Output:     Profiles of radiative fluxes, heating
!                 rates for air and particles; verticaly ! spelled wrong so I don't grep a certain word
!                 integrated optical depths; and albedos (which
!                 are all loaded into interface common block).
!     **************************************************************
!
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
      REAL qrad(NL+1),alb_tomi,alb_toai

      integer, parameter :: nwave_alb = NTOTAL
      real wavea(nwave_alb),albedoa(nwave_alb),t(NZ),p_cgs(NZ)
      real maxopd(nwave_alb)
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL) :: Beta_V
      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL,TAUGAS,TAUAER
      real incident_starlight_fraction
      integer jflip, solar_calculation_indexer

      real, dimension(NTOTAL,NDBL) :: SLOPE
      real pr(NLAYER), p_pass(NLAYER)

      real dpg(nl+1), pbar(nl+1)
      real dpgsub(2*nl+2), pbarsub(2*nl+2)

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

      integer L, J, K

!     Reset flag for computation of solar fluxes
      if (incident_starlight_fraction .gt. 1e-5) then
          ISL = 1
      else
          ISL = 0
      endif

      do k = 1,nvert
        p_cgs(k) = pr(k)*10.  ! to convert from Pa to Dyne cm^2
      enddo

      TT(1) = t(1) ! MALSKY ADDED
      DO J = 2, NVERT
          TT(J) = T(J-1) * ((p_pass(J)*10.0)/p_cgs(J-1)) ** (log(T(J)/T(J-1))/log(p_cgs(J)/p_cgs(J-1)))
      END DO

      TT(1)=((T(1)-TT(2))/log(p_cgs(1)/p_pass(2)))*log(p_pass(1)/p_cgs(1))+T(1)
      TT(NLAYER) = T(NVERT) * ((p_pass(NLAYER)*10)/p_cgs(NVERT)) **
     &             (log(T(NVERT)/T(NVERT-1))/log(p_cgs(NVERT)/p_cgs(NVERT-1)))


!     Solar zenith angle
      u0 = incident_starlight_fraction

!     SURFACE REFLECTIVITY AND EMISSIVITY
!     Hack: use spectrally dependent surface albedo
      DO 20 L =  1,NSOLP
         RSFX(L) = ALBSW
         EMIS(L) =  1.0 - RSFX(L)
 20   CONTINUE

!...Hack: specify EMIS based on RSFX rather than visa versa
      DO 30 L =  NSOLP+1,NTOTAL
         EMIS(L) =  EMISIR
         RSFX(L) = 1.0 - EMIS(L)

         if( wave(nprob(L)).gt.wavea(nwave_alb) ) then
             rsfx(L) = albedoa(nwave_alb)
         endif

         EMIS(L) = 1.0 - RSFX(L)
 30   CONTINUE

!     CALCULATE THE OPTICAL PROPERTIES
      IF(AEROSOLCOMP .EQ. 'All') THEN
          !CALL OPPRMULTI(TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG, p_pass,
!     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
!     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK)

          ! This one only works with 50 layers
          IF (NL .eq. 50) THEN
              CALL DOUBLEGRAY_OPPRMULTI(TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
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
     &  qrad,alb_tomi,alb_toai, num_layers)
          ELSE
              write(*,*) 'Youre doing the old cloud version with NL not equal to 50'
              stop
          END IF
      ELSE
          write(*,*) 'ERROR! Dont run without aerosols'
          STOP
      ENDIF





      SLOPE(:,:) = 0.0
      CALL OPPR1(TAUL, SLOPE, t,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, EMISIR,
     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
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
     &  qrad,alb_tomi,alb_toai, num_layers)


!     IF NO INFRARED SCATTERING THEN SET INDEX TO NUMBER OF SOLAR INTERVALS
      IF(IRS .EQ. 0) THEN
          LLA  =  NSOLP
          write(*,*) "Something funny is going on, why no IR scattering?"
          stop
      ENDIF

!     IF EITHER SOLAR OR INFRARED SCATTERING CALCULATIONS ARE REQUIRED
!     GET TWO STREAM CODE AND FIND THE SOLUTION
      IF(incident_starlight_fraction .gE. 0 .OR. IRS .NE. 0) THEN
          CALL TWOSTR(TAUL, solar_calculation_indexer,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
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
     &  qrad,alb_tomi,alb_toai, num_layers)

          CALL ADD(TAUL, solar_calculation_indexer, SLOPE,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
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
     &  qrad,alb_tomi,alb_toai, num_layers)
      ENDIF


!     IF INFRARED CALCULATIONS ARE REQUIRED THEN NEWFLUX1 FOR MORE ACCURATE SOLUTION

      IF(IR .NE. 0) THEN
          CALL NEWFLUX1(TAUL,SLOPE,LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &                  EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
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
     &  qrad,alb_tomi,alb_toai, num_layers)
      ENDIF

!     CLOUD FRACTION
!     NOW, IF WE ARE INCLUDING AEROSOLS, AND WE WOULD LIKE A  CLOUD
!     FRACTION LESS THAN UNITY, THEN RECOMPUTE THESE FLUXES FOR A
!     CLEAR SKY AND COMBINE IN A WEIGHTED AVERAGE ASSUMING MAXIMUM OVERLAP.

!     HERE WE TAKE THE DOUBLE RESOLUTION FLUXES AND EXTRACT JUST
!     THOSE THAT CORRESPOND TO THE STANDARD (VISIBLE) PRESSURE
!     LEVELS. THESE VALUES SHOULD BE SUPERIOR TO THOSE COMPUTED
!     WITHOUT DOUBLING

      DO L = NSOLP+1, NTOTAL
          K     =  1
          DO        J     =  1,NLAYER
            FNET(L,J)      =  DIRECTU(L,k)-DIREC(L,k)
            DIRECTU(L,J)   =  DIRECTU(L,K)
            DIREC(L,J)     =  DIREC(L,K)
            OPD(L,J)       =  OPD(L,K)
            TAUL(L,J)      =  TAUL(L,K)
            W0(L,J)        =  W0(L,K)
            G0(L,J)        =  G0(L,K)
            ug0(L,J)       =  ug0(L,K)
            uOPD(L,J)      =  uOPD(L,K)
            uTAUL(L,j)     =  uTAUL(L,K)
            uW0(L,J)       =  uW0(L,k)
            TMIU(L,J)      =  TMIU(L,k)
            TMID(L,J)      =  TMID(L,k)
            K     =  K+2
          ENDDO
      END DO



!     ATTENTION! THE FOLLOWING IS A MODEL-SPECIFIC MODIFICATION:
!     HERE WE PRESCRIBE THE BOTTOM BOUNDARY CONDITION NET FLUX IN THE IR.
!     BE AWARE: IT ALSO AFFECTS THE UPWARD FLUX AT THE BASE IN NEWFLUX.

      DO L = NSOLP+1,NTOTAL
          FNET(L,NLAYER)=FBASEFLUX
      END DO

!     SINCE WE ARE PLACING A CONSTRAINT ON THE NET FLUX AT THE BOTTOM OF
!     THE MODEL (BY DEFINING  NET FLUX AT THE BASE TO EQUAL FBASEFLUX)
!     TO BE SELF CONSISTENT, WE CAN REDIFINE THE UPWARD FLUX FROM THE
!     SURFACE TO BE NETFLUX MINUS DOWNWARD FLUX

!     HERE WE DERIVE THE UPWARD FLUX FROM THE NET FLUX, SELF CONSISTENT
!     WITH BOTTOM BOUNDARY CONDITION

!     MALSKY CHECK THAT THIS IS THE RIGHT BOUNDARY
!     MAYBE IT SHOULD BE NDBL
      DO L =  NSOLP+1,NTOTAL
          DIRECTU(L,NLAYER) = FBASEFLUX+DIREC(L,NLAYER)
      END DO

      DO J = 1, NLAYER
           DO L = 1, NTOTAL
               IF (L .LE. NSOLP) THEN
                   FNET(L,J) = FNET(L,J) * Beta_V(L)
               ELSE
                   FNET(L,J)    = FNET(L,J)    * Beta_IR(L - NSOLP)
                   DIREC(L,J)   = DIREC(L,J)   * Beta_IR(L - NSOLP)
                   DIRECTU(L,J) = DIRECTU(L,J) * Beta_IR(L - NSOLP)
               END IF
           END DO
      END DO


!     CALCULATE INFRAFRED AND SOLAR HEATING RATES (DEG/DAY),
      DO 500 J      =  1,NVERT
          HEATS(J)   =  0.0
          HEATI(J)   =  0.0

          TERM1      =  FDEGDAY/(DPG(J+1)*G)

          IF(incident_starlight_fraction.ge. 0) THEN
              DO 480 L     =  1,NSOLP
                  HEATS(J)   =  HEATS(J)+(FNET(L,J+1)-FNET(L,J)) * TERM1
 480          CONTINUE
          ENDIF

          IF (IR .NE. 0) THEN
              DO L    =  NSOLP+1,NTOTAL
                  HEATI(J)   =  HEATI(J)+(FNET(L,J+1)-FNET(L,J))*TERM1
              END DO
          ENDIF

          HEAT(J) = HEATS(J) + HEATI(J)

!         Load heating rates [deg_K/s] into interface common block
          heats_aerad(j) =  heats(j)/scday
          heati_aerad(j) =  heati(j)/scday

500   CONTINUE

      write(*,*) HEAT(:)


!     Load layer averages of droplet heating rates into interface common block
!     Calculate some diagnostic quantities (formerly done in radout.f) and
!     load them into the interface common block.  None of these presently
!     influence any microphysical processes -- hence, the following code
!     only needs to be executed before the aerosol model writes its output.
!     Not all of the calculated quantities are presently being
!     loaded into the interface common block.
!     Load optical depths into interface common block
!     <tsLu> and <total_downwelling> are total upwelling and downwelling solar
!     fluxes at top-of-atmosphere

      tsLu = 0.
      total_downwelling = 0.

!     <fupbs>, <fdownbs>, and <fnetbs> are total upwelling, downwelling,
!     and net solar fluxes at grid boundaries

      do 507 j = 1, nlayer
          fupbs(j)   = 0.
          fdownbs(j) = 0.
          fnetbs(j)  = 0.
          fdownbs2(j)= 0.
 507  continue

!     <fsLu> and <fsLd> are upwelling, downwelling, and net
!     solar fluxes at top-of-atmosphere (spectrally-resolved)
!     <alb_toa> is albedo at top-of-atmosphere (spectrally-resolved)

      do 509 i = 1, nsoL
          fsLu(i)    = 0.0
          fsLd(i)    = 0.0
          alb_toa(i) = 0.0
509   continue
!
!     <alb_tomi> and <alb_toai> are total solar albedos at top-of-model
!     and top-of-atmosphere
!
      alb_tomi = 0.
      alb_toai = 0.
!
!     CALCULATE SOLAR ABSORBED BY GROUND, SOLNET, AND UPWARD AND
!     DOWNWARD LONGWAVE FLUXES AT SURFACE
!
      SOLNET  = 0.0
      IF (solar_calculation_indexer .gE. 0) THEN
          DO 510 L       =  1,NSOLP
              SOLNET  = SOLNET - FNET(L,NLAYER)
              fp      = ck1(L,1) * eL2(L,1) - ck2(L,1) * em2(L,1) + cp(L,1)
              fsLu(L) = fsLu(L) + fp

              do 510 j = 1, NLAYER
                  fp  =  ck1(L,j) * eL1(L,j) + ck2(L,j) * em1(L,j) + cpb(L,j)
                  fm  =  ck1(L,j) * eL2(L,j) + ck2(L,j) * em2(L,j) + cmb(L,j)

                  fupbs(j)    = fupbs(j)    + fp * Beta_V(L)
                  fdownbs2(j) = fdownbs2(J) + fm * Beta_V(L)
                  fnetbs(j)   = fnetbs(j)   + fnet(L,j)

                  if (L .eq. nsolp) then
                      fdownbs(J) = (fupbs(j) - fnetbs(j))
                  endif

510      CONTINUE


          do  i = 1, nsoL
              fsLd(i) = psol_aerad*incident_starlight_fraction
              alb_toa(i) = fsLu(i)/fsLd(i)
              tsLu = tsLu + fsLu(i)
              total_downwelling = total_downwelling + fsLd(i)
          END DO

          alb_tomi = fupbs(1)/fdownbs(1)
          alb_toai = tsLu/total_downwelling
!
!         Load fluxes into interface common block
!
          do j = 1, nlayer
              fsl_up_aerad(j) = fupbs(j)
              fsl_dn_aerad(j) = fdownbs(j)
          enddo
      ENDIF


!     <tiru> is total upwelling infrared flux at top-of-atmosphere;
!     <fupbi>, <fdownbi>, and <fnetbi> are total upwelling, downwelling,
!     and net
!     infrared fluxes at grid boundaries

      tiru = 0.0

      do 606 j = 1, nlayer
          fupbi(j)    =  0.0
          fdownbi(j)  =  0.0
          fnetbi(j)   =  0.0
606   continue

!     <firu> is upwelling infrared flux at top-of-atmosphere
!     (spectrally-resolved)

      do 609 i = 1, nir
          firu(i) = 0.
 609  continue


      IF (IR .NE. 0) THEN
          DO 520 L        =  NSOLP+1,NTOTAL
             firu(L-nsol ) = firu( L-nsol ) + directu(L,1)

             do 520 j = 1, nlayer
                 fupbi(j)   = fupbi(j)   + (directu(L,j))
                 fdownbi(j) = fdownbi(j) + (direc(L,j))
                 fnetbi(j)  = fnetbi(j)  + (directu(L,j) - direc(L,j))
 520      CONTINUE


          do 529 i = 1, nir
              tiru = tiru + firu(i)
 529      continue

!         Load fluxes into interface common block

          do j = 1, nlayer
              jflip=nlayer+1-j
              fir_up_aerad(j)  = fupbi(jflip)
              fir_dn_aerad(j)  = fdownbi(jflip)
              fir_net_aerad(j) = fnetbi(jflip)
              fsl_up_aerad(j)  = fupbs(jflip)
              fsl_dn_aerad(j)  = fdownbs(jflip)
              fsl_net_aerad(j) = fnetbs(jflip)
          enddo
      ENDIF


C     RFLUXES  Array to hold fluxes at top and bottom of atmosphere
C     1st index - flux 1=SW, 2=LW
C     2nd index - Direction 1=DN, 2=UP
C     3rd index - Where 1=TOP, 2=SURFACE

      if (NSOLP .gt. 1) then
          RFLUXES_aerad(1,1,1) = fsl_dn_aerad(NLAYER) * Beta_V(1) +
     &                           fsl_dn_aerad(NLAYER) * Beta_V(2) +
     &                           fsl_dn_aerad(NLAYER) * Beta_V(3)

          RFLUXES_aerad(1,1,2) = fsl_dn_aerad(1)/(1.0-ALBSW) * Beta_V(1) +
     &                           fsl_dn_aerad(1)/(1.0-ALBSW) * Beta_V(2) +
     &                           fsl_dn_aerad(1)/(1.0-ALBSW) * Beta_V(3)

          RFLUXES_aerad(1,2,1) = fsl_up_aerad(NLAYER) * Beta_V(1) +
     &                           fsl_up_aerad(NLAYER) * Beta_V(2) +
     &                           fsl_up_aerad(NLAYER) * Beta_V(3)

          RFLUXES_aerad(1,2,2) = RFLUXES_aerad(1,1,2)*ALBSW * Beta_V(1) +
     &                           RFLUXES_aerad(1,1,2)*ALBSW * Beta_V(2) +
     &                           RFLUXES_aerad(1,1,2)*ALBSW * Beta_V(3)

          RFLUXES_aerad(2,1,1) = fir_dn_aerad(NLAYER) * Beta_IR(1) +
     &                           fir_dn_aerad(NLAYER) * Beta_IR(2)

          RFLUXES_aerad(2,1,2) = fir_dn_aerad(1) * Beta_IR(1) +
     &                           fir_dn_aerad(1) * Beta_IR(2)


          RFLUXES_aerad(2,2,1) = fir_up_aerad(NLAYER) * Beta_IR(1) +
     &                           fir_up_aerad(NLAYER) * Beta_IR(2)

          RFLUXES_aerad(2,2,2) = fir_up_aerad(1) * Beta_IR(1) +
     &                           fir_up_aerad(1) * Beta_IR(2)
      else
          RFLUXES_aerad(1,1,1)=fsl_dn_aerad(NLAYER)   ! SW down top
          RFLUXES_aerad(1,1,2)=fsl_dn_aerad(1)/(1.0-ALBSW)   ! SW down bottom
          RFLUXES_aerad(1,2,1)=fsl_up_aerad(NLAYER)  ! SW up top
          RFLUXES_aerad(1,2,2)=RFLUXES_aerad(1,1,2)*ALBSW   ! SW up bottom

          RFLUXES_aerad(2,1,1)=fir_dn_aerad(NLAYER)   ! LW down top
          RFLUXES_aerad(2,1,2)=fir_dn_aerad(1)       ! LW down bottom
          RFLUXES_aerad(2,2,1)=fir_up_aerad(NLAYER)       ! LW up top
          RFLUXES_aerad(2,2,2)=fir_up_aerad(1)   ! LW up bottom
      end if

      return
      END

