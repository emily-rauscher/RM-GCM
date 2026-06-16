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
     &  qrad,alb_tomi,alb_toai, num_layers, SLOPE, Y1, Y2, Y4, Y8, A1, A2, A3, A4, A5, A7, Y5,
     &  PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom, kount, ITSPD)

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
      use corrkmodule, only : TS_CORRK, PS_CORRK, TS_LOG_CORRK, 
     &           PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC, INT_SPEC, TAURAY_PER_DPG,
     &           OPAC_CORRK, PLANCK_INTS, PLANCK_TS, NWNO
      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, kount, itspd
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL), SOL(NTOTAL)
      REAL RAYPERBAR(NTOTAL),WEIGHT(NTOTAL)
      REAL GOL(NTOTAL,2*NL+2), WOL(NTOTAL,2*NL+2), WAVE(5+1), TT(NL+1), Y3(NTOTAL,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NTOTAL), PTEMPT(NTOTAL), G0(NTOTAL,2*NL+2), OPD( NTOTAL,2*NL+2), PTEMP(NTOTAL,2*NL+2)
      REAL uG0(NTOTAL,2*NL+2), uTAUL(NTOTAL,2*NL+2), W0(NTOTAL,2*NL+2), uW0(NTOTAL,2*NL+2), uopd(NTOTAL,2*NL+2),  U1S( NTOTAL)
      REAL U1I(NTOTAL), TOON_AK(NTOTAL,2*NL+2), B1(NTOTAL,2*NL+2), B2(  NTOTAL,2*NL+2), EE1( NTOTAL,2*NL+2), EM1(NTOTAL,2*NL+2)
      REAL EM2(NTOTAL,2*NL+2), EL1( NTOTAL,2*NL+2), EL2(NTOTAL,2*NL+2), GAMI(NTOTAL,2*NL+2), AF(NTOTAL,4*NL+4)
      REAL BF(NTOTAL,4*NL+4), EF(NTOTAL,4*NL+4), SFCS(NTOTAL), B3(NTOTAL,2*NL+2), CK1(NTOTAL,2*NL+2), CK2(NTOTAL,2*NL+2)
      REAL CP(NTOTAL,2*NL+2), CPB(NTOTAL,2*NL+2), CM(NTOTAL,2*NL+2), CMB(NTOTAL,2*NL+2), DIRECT(NTOTAL,2*NL+2), EE3(NTOTAL,2*NL+2)
      REAL EL3(NTOTAL,2*NL+2), FNET(NTOTAL,2*NL+2), TMI(NTOTAL,2*NL+2), AS(NTOTAL,4*NL+4), DF(NTOTAL,4*NL+4)
      REAL DS(NTOTAL,4*NL+4), XK(NTOTAL,4*NL+4), DIREC(NTOTAL,2*NL+2), DIRECTU(NTOTAL,2*NL+2), DINTENT(NTOTAL,3,2*NL+2)
      REAL UINTENT(NTOTAL,3,2*NL+2), TMID(NTOTAL,2*NL+2), TMIU(NTOTAL,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai

      REAL ttsub(2*NL+2)

      real wavea(NTOTAL),albedoa(NTOTAL),t(NZ)
      real maxopd(NTOTAL)
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL) :: Beta_V
      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL,TAUGAS,TAUAER
      real incident_starlight_fraction
      integer solar_calculation_indexer

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

      REAL, DIMENSION(NTOTAL,3,2*NL+2) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(NTOTAL,2*NL+2)   :: A1, A2, A3, A4, A5, A7, Y5

      REAL PI0_TEMP(NTOTAL, NL+1, 13)
      REAL G0_TEMP(NTOTAL, NL+1, 13)
      REAL tauaer_temp(NTOTAL, NL+1, 13)
      INTEGER j1
      real denom
      ! Corr-K common block:
    !   COMMON/CORRKGAS/OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC,
    !  &      INT_SPEC, TAURAY_PER_DPG
    !   REAL :: OPAC_CORRK(73, 20, NWNO, 8)
    !   REAL :: TS_CORRK(73), PS_CORRK(20), TS_LOG_CORRK(73), PS_LOG_CORRK(20), WGTS_CORRK(8) 
    !   REAL :: WNO_EDGES(NWNO+1), WNO_CTRS(NWNO), STEL_SPEC(NWNO), INT_SPEC(NWNO)
    !   REAL :: TAURAY_PER_DPG(73, 20, NWNO)

      integer L, J, K, chan_idx, stel_idx

      ! ADDING THESE
      real, dimension(NL+1) :: Tl, pl, pe, dpe
      real, dimension(NL+1) :: lTl, lpl, lpe, Te
      logical bezier_interpolation
      bezier_interpolation = .FALSE.

!     Reset flag for computation of solar fluxes
      if (incident_starlight_fraction .lt. 1e-10) then
          ISL = 0
      else
          ISL = 1
      endif

      IF (bezier_interpolation .eq. .TRUE.) THEN
          do J = 1, NL+1
             pe(J) = p_pass(J)
          end do

          DO J = 1, NL
              dpe(J) = pe(J+1) - pe(J)
              pl(J) = dpe(J) / log(pe(J+1)/pe(J))
          END DO

          if (tt(1) .ge. 100.0) then
              DO J = 1, NLAYER
                  Tl(J) = tt(J)
              END DO
          else
              DO J = 1, NL + 1
                  Tl(J) = t(J)
              END DO
          end if

          pl(NLAYER)  = 10.0 ** (LOG10(pl(NLAYER-1))  + (LOG10(pl(NLAYER-1))  - LOG10(pl(NLAYER-2))))
          Tl(NLAYER)  = Tl(NLAYER-1) + ABS(Tl(NLAYER-1) - Tl(NLAYER-2)) / 2.0

          lTl(:) = log10(Tl(:))
          lpl(:) = log10(pl(:))
          lpe(:) = log10(pe(:))

           do i = 2, NLAYER - 1
              !call bezier_interp(lpl(i-1:i+1), lTl(i-1:i+1), 3, lpe(i), Te(i))
              !Te(i) = 10.0 ** (Te(i))

              call bezier_interp(lpl(i-1:i+1), lTl(i-1:i+1), 3, lpe(i), TT(i))
              TT(i) = 10.0 ** (TT(i))
           end do

!           Te(1)  = 10.0 ** (log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/Te(2)))
!           Te(NLAYER) = 10.0 ** (log10(Tl(NLAYER-1)) + (log10(pe(NLAYER)/pe(NLAYER-1))/log10(pl(NLAYER-1)/pe(NLAYER-1)))
!     &                  * log10(Tl(NLAYER-1)/Te(NLAYER-1)))

           TT(1)  = 10.0 ** (log10(Tl(1)) + (log10(pe(1)/pe(2))/log10(pl(1)/pe(2))) * log10(Tl(1)/TT(2)))
           TT(NLAYER) = 10.0 ** (log10(Tl(NLAYER-1)) + (log10(pe(NLAYER)/pe(NLAYER-1))/log10(pl(NLAYER-1)/pe(NLAYER-1)))
     &                  * log10(Tl(NLAYER-1)/TT(NLAYER-1)))
      END IF

      !     create a T array for the double resolution IR by combinging the
      !     layer center and edge temperatures
      K  =  1
      DO J  = 1, NDBL-1,2
         L  =  J
         TTsub(L) = TT(K)
         L  =  L+1
         TTsub(L) = T(K)
         K  =  K+1
      END DO

!     Solar zenith angle
      u0 = incident_starlight_fraction

!     SURFACE REFLECTIVITY AND EMISSIVITY
!     Hack: use spectrally dependent surface albedo
      DO 20 L =  1,NSOL
         RSFX(L) = ALBSW
         EMIS(L) =  1.0 - RSFX(L)
 20   CONTINUE

!...Hack: specify EMIS based on RSFX rather than visa versa
      DO 30 L =  NSOL+1,NTOTAL
         EMIS(L) =  EMISIR
         RSFX(L) = 1.0 - EMIS(L)

         if( wave(nprob(L)).gt.wavea(NTOTAL) ) then
             rsfx(L) = albedoa(NTOTAL)
         endif

         EMIS(L) = 1.0 - RSFX(L)
 30   CONTINUE

!     CALCULATE THE OPTICAL PROPERTIES
    !   IF (AEROSOLCOMP.EQ. 'standard') THEN
      IF ((opacity_method .EQ. 'dogray') .or. (opacity_method .EQ. 'picket')) THEN
        ! write(*,*) "Using dogray or picket fence opacity method"
        CALL OPPRMULTI(TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG,
     &                   LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &                   EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
     &                   SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
     &                   GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
     &                   WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
     &                   uG0,uTAUL,W0,uW0,uopd,U1S,
     &                   U1I,TOON_AK,B1,B2,EE1,EM1,
     &                   EM2,EL1,EL2,GAMI,AF,
     &                   BF,EF,SFCS,B3,CK1,CK2,
     &                   CP,CPB,CM,CMB,DIRECT,EE3,
     &                   EL3,FNET,TMI,AS,DF,
     &                   DS,XK,DIREC,DIRECTU,DINTENT,
     &                   UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
     &                   tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
     &                   fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &                   qrad,alb_tomi,alb_toai, p_pass,
     &                   PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom,kount,itspd)
      ELSE
        ! write(*,*), "TAUGAS:", TAUGAS(71,1:NLAYER)
        CALL OPPRMULTI_CORRK(TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG,
     &                   LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &                   EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
     &                   SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
     &                   GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
     &                   WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
     &                   uG0,uTAUL,W0,uW0,uopd,U1S,
     &                   U1I,TOON_AK,B1,B2,EE1,EM1,
     &                   EM2,EL1,EL2,GAMI,AF,
     &                   BF,EF,SFCS,B3,CK1,CK2,
     &                   CP,CPB,CM,CMB,DIRECT,EE3,
     &                   EL3,FNET,TMI,AS,DF,
     &                   DS,XK,DIREC,DIRECTU,DINTENT,
     &                   UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
     &                   tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
     &                   fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &                   qrad,alb_tomi,alb_toai, p_pass,
     &                   PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom,kount,itspd)
        ! write(*,*) ""
        ! write(*,*) ""
        ! write(*,*) ""
        ! write(*,*) ""
        ! write(*,*) ""
        ! write(*,*) "TAUL:", TAUL(71,1:NLAYER)
        ! write(*,*) "TAUAER:", TAUAER(71,1:NLAYER)
        ! write(*,*) ""
        ! write(*,*) "TAUGAS:", TAUGAS
        ! write(*,*) "TAUAER:", TAUAER
        ! stop
      ENDIF
    !       write(*,*) TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG,
    !  &                   LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
    !  &                   EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
    !  &                   SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
    !  &                   GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
    !  &                   WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
    !  &                   uG0,uTAUL,W0,uW0,uopd,U1S,
    !  &                   U1I,TOON_AK,B1,B2,EE1,EM1,
    !  &                   EM2,EL1,EL2,GAMI,AF,
    !  &                   BF,EF,SFCS,B3,CK1,CK2,
    !  &                   CP,CPB,CM,CMB,DIRECT,EE3,
    !  &                   EL3,FNET,TMI,AS,DF,
    !  &                   DS,XK,DIREC,DIRECTU,DINTENT,
    !  &                   UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
    !  &                   tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
    !  &                   fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
    !  &                   qrad,alb_tomi,alb_toai, p_pass,
    !  &                   PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom,kount
        !   STOP
    !   ELSE IF (AEROSOLCOMP.EQ. 'none') THEN
    !       ! TK no aerosols, just gas
    !       DO J = 1,NDBL
    !         DO L = 1,NSOL
    !             TAUL(L,J) = TAUGAS(L,J) + TAURAY(L,J)
    !             W0(L,J) = TAURAY(L,J)/TAUL(L,J)
    !             G0(L,J) = 0.0 ! Rayleigh g=0
    !             OPD(L,J) = TAUGAS(L,J) + TAURAY(L,J) ! I think this is a test thing Mike put in that was never removed
    !         END DO
    !       END DO
    !       DO I = 1,NGAUSS
    !         DO L = NSOL+1,NTOTAL
    !             Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
    !         END DO
    !       END DO
    !       write(*,*) TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG,
    !  &                   LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
    !  &                   EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
    !  &                   SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
    !  &                   GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
    !  &                   WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
    !  &                   uG0,uTAUL,W0,uW0,uopd,U1S,
    !  &                   U1I,TOON_AK,B1,B2,EE1,EM1,
    !  &                   EM2,EL1,EL2,GAMI,AF,
    !  &                   BF,EF,SFCS,B3,CK1,CK2,
    !  &                   CP,CPB,CM,CMB,DIRECT,EE3,
    !  &                   EL3,FNET,TMI,AS,DF,
    !  &                   DS,XK,DIREC,DIRECTU,DINTENT,
    !  &                   UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
    !  &                   tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
    !  &                   fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
    !  &                   qrad,alb_tomi,alb_toai, p_pass,
    !  &                   PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom,kount
    !       STOP
        !   write(*,*) "ERROR! Thomas hasn't coded a workaround for corr-k clouds yet!"
        !   STOP
    !   ELSE
    !       write(*,*) 'ERROR! Something is off about your version. Check with Isaac'
    !       STOP
    !   ENDIF

      SLOPE(:,:) = 0.0
      DS(:,:)    = 0.0
      DF(:,:)    = 0.0
      
      IF ((opacity_method .EQ. 'dogray') .or. (opacity_method .EQ. 'picket')) THEN
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
      ELSE IF (opacity_method .EQ. 'correk') then
        CALL OPPR1_CORRK(TAUL, SLOPE, t,
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
      ENDIF

!     IF NO INFRARED SCATTERING THEN SET INDEX TO NUMBER OF SOLAR INTERVALS
      IF(IRS .EQ. 0) THEN
          LLA  =  NSOL
          write(*,*) "Something funny is going on, why no IR scattering?"
          stop
      ENDIF

      B1 = 0.
      B2 = 0.
      EL1 = 0.
      EL2 = 0.
      EM1 = 0.
      EM2 = 0.
      ck1 = 0.
      ck2 = 0.
      cpb = 0.
      cp  = 0.

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
    !   write(*,*) "before RADD: ", FNET(3,:NLAYER+1)
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
    !   write(*,*) "after RADD: ", FNET(3,:NLAYER+1)

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
     &  qrad,alb_tomi,alb_toai, num_layers, Y1, Y2, Y4, Y8, A1, A2, A3, A4, A5, A7, Y5)
      ENDIF


!     CLOUD FRACTION
!     NOW, IF WE ARE INCLUDING AEROSOLS, AND WE WOULD LIKE A  CLOUD
!     FRACTION LESS THAN UNITY, THEN RECOMPUTE THESE FLUXES FOR A
!     CLEAR SKY AND COMBINE IN A WEIGHTED AVERAGE ASSUMING MAXIMUM OVERLAP.

!     HERE WE TAKE THE DOUBLE RESOLUTION FLUXES AND EXTRACT JUST
!     THOSE THAT CORRESPOND TO THE STANDARD (VISIBLE) PRESSURE
!     LEVELS. THESE VALUES SHOULD BE SUPERIOR TO THOSE COMPUTED
!     WITHOUT DOUBLING

      DO L = NSOL+1, NTOTAL ! IR un-doubling?
          K =  1
          DO J =  1,NLAYER
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
            K = K+2
          ENDDO
      END DO
    !   write(*,*) "DIREC:" , DIREC(NSOL,:NLAYER+1)
    !   write(*,*) "DIRECTU:" , DIRECTU(NSOL,:NLAYER+1)

    !   WRITE(*,*) "FNET1: ", FNET(1,:NLAYER+1)
    !   WRITE(*,*) "FNET8: ", FNET(8,:NLAYER+1)
    !   WRITE(*,*) "FNET9: ", FNET(9,:NLAYER+1)

!     MALSKY CHECK THAT THIS IS THE RIGHT BOUNDARY
    !   write(*,*) "MU:" , U0
    !   DO L = 1, NSOL
    !     write(*,*) "FNET1: ", L, FNET(L,:NLAYER+1)
    !     ! write(*,*) "DIREC: ", L, DIREC(L,:NLAYER+1)
    !     write(*,*) "TAUL1: ", L, TAUL(L,:NLAYER+1)
    !     write(*,*) "OPD: ", L, OPD(L,:NLAYER+1)
    !   END DO
      IF ((opacity_method .EQ. 'dogray') .or. (opacity_method .EQ. 'picket')) THEN
!     ATTENTION! THE FOLLOWING IS A MODEL-SPECIFIC MODIFICATION:
!     HERE WE PRESCRIBE THE BOTTOM BOUNDARY CONDITION NET FLUX IN THE IR.
!     BE AWARE: IT ALSO AFFECTS THE UPWARD FLUX AT THE BASE IN NEWFLUX.
        DO L = NSOL+1,NTOTAL
          FNET(L,NLAYER)=FBASEFLUX
        END DO

!     SINCE WE ARE PLACING A CONSTRAINT ON THE NET FLUX AT THE BOTTOM OF
!     THE MODEL (BY DEFINING  NET FLUX AT THE BASE TO EQUAL FBASEFLUX)
!     TO BE SELF CONSISTENT, WE CAN REDIFINE THE UPWARD FLUX FROM THE
!     SURFACE TO BE NETFLUX MINUS DOWNWARD FLUX

!     HERE WE DERIVE THE UPWARD FLUX FROM THE NET FLUX, SELF CONSISTENT
!     WITH BOTTOM BOUNDARY CONDITION
        
        DO L = NSOL+1,NTOTAL
          DIRECTU(L,NLAYER) = FBASEFLUX+DIREC(L,NLAYER)
        END DO
        ! Prepare to combine channels
        DO J = 1, NLAYER
            DO L = 1, NTOTAL
                IF (L .LE. NSOL) THEN
                    FNET(L,J) = FNET(L,J) * Beta_V(L)
                ELSE
                    FNET(L,J)    = FNET(L,J)    * Beta_IR(L - NSOL)
                    DIREC(L,J)   = DIREC(L,J)   * Beta_IR(L - NSOL)
                    DIRECTU(L,J) = DIRECTU(L,J) * Beta_IR(L - NSOL)
                END IF
            END DO
       END DO
      ELSE IF (opacity_method .EQ. 'correk') THEN
        DO L = NSOL+1,NTOTAL
          chan_idx = MODULO(L-1,8)+1
        !   stel_idx = (L - chan_idx)/8 + 1
          stel_idx = MODULO((L-chan_idx)/8, NWNO) + 1
        !   write(*,*) "stel_idx: ", stel_idx
        !   write(*,*) "chan_idx: ", chan_idx
        !   write(*,*) "FNET: ", FNET(L,NLAYER)
          FNET(L,NLAYER)= FBASEFLUX * INT_SPEC(stel_idx)
        !   write(*,*) "FNET2: ", FNET(L,NLAYER)
        !   write(*,*) " "
          DIRECTU(L,NLAYER) = (FBASEFLUX * INT_SPEC(stel_idx)) + DIREC(L,NLAYER)
        END DO
        ! Prepare to combine channels
        DO J = 1, NLAYER
            DO L = 1, NTOTAL
                chan_idx = MODULO(L-1,8)+1
                stel_idx = MODULO((L-chan_idx)/8, NWNO) + 1
                ! stel_idx = (L - chan_idx)/8 + 1
                IF (L .LE. NSOL) THEN
                    FNET(L,J) = FNET(L,J) * WGTS_CORRK(chan_idx) ! * STEL_SPEC(stel_idx)
                ELSE
                    FNET(L,J)    = FNET(L,J)    * WGTS_CORRK(chan_idx)
                    DIREC(L,J)   = DIREC(L,J)   * WGTS_CORRK(chan_idx)
                    DIRECTU(L,J) = DIRECTU(L,J) * WGTS_CORRK(chan_idx)
                END IF
            END DO
       END DO
       ENDIF
    !    DO L = 1, NSOL
    !     write(*,*) "FNET2: ", L, FNET(L,:NLAYER+1)
    !     ! write(*,*) "DIREC: ", L, DIREC(L,:NLAYER+1)
    !     write(*,*) "TAUL2: ", L, TAUL(L,:NLAYER+1)
    !    END DO
      
    !   write(*,*) "FBASEFLUX: ", FBASEFLUX
    !   write(*,*) "intspec: ", INT_SPEC


!     CALCULATE INFRAFRED AND SOLAR HEATING RATES (DEG/DAY),
      DO 500 J      =  1,NVERT
          HEATS(J)   =  0.0
          HEATI(J)   =  0.0

          TERM1      =  FDEGDAY/(DPG(J+1)*G)

          IF(incident_starlight_fraction.ge. 0) THEN
              DO 480 L     =  MAX(solar_calculation_indexer,MINWNOSTEL*8+1),NSOL
                  HEATS(J)   =  HEATS(J)+(FNET(L,J+1)-FNET(L,J)) * TERM1
 480          CONTINUE
          ENDIF

          IF (IR .NE. 0) THEN
              DO L    =  NSOL+1,NTOTAL
                  HEATI(J)   =  HEATI(J)+(FNET(L,J+1)-FNET(L,J))*TERM1
              END DO
          ENDIF

          HEAT(J) = HEATS(J) + HEATI(J)

!         Load heating rates [deg_K/s] into interface common block
          heats_aerad(j) =  heats(j)/scday
          heati_aerad(j) =  heati(j)/scday
 


500   CONTINUE

    !   write(*,*) "HEATS: ", heats_aerad
    !   write(*,*) "HEATI: ", heati_aerad

!     Here we Calculate (4 * pi * mean_intensity) for the IR.
      IF (IR .NE. 0) THEN
        DO J = 1, NVERT
          DO L = NSOL+1, NTOTAL
            TMI(L,J) = TMIU(L,J)+TMID(L,J)
          end do
        end do
      ENDIF


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

      do 509 i = MAX(solar_calculation_indexer,MINWNOSTEL*8), nsoL
          fsLu(i)    = 0.0
          fsLd(i)    = 0.0
509   continue
!
!     <alb_tomi> and <alb_toai> are total solar albedos at top-of-model
!     and top-of-atmosphere
!
      alb_tomi = 0.
      alb_toai = 0.
      total_downwelling = 0.
      fupbs(1) = 0.0
      fdownbs(1) = 0.0




!     CALCULATE SOLAR ABSORBED BY GROUND, SOLNET, AND UPWARD AND
!     DOWNWARD LONGWAVE FLUXES AT SURFACE

      
      SOLNET   = 0.0
      IF (ISL .GT. 0) THEN
          if ((opacity_method .eq. 'picket') .or. (opacity_method .eq. 'dogray')) then
            DO L       =  MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
                SOLNET  = SOLNET - FNET(L,NLAYER)
                fp      = (ck1(L,1) * eL2(L,1) - ck2(L,1) * em2(L,1) + cp(L,1)) * Beta_V(L)

                fsLu(L) = fsLu(L) + fp

                do j = 1, NLAYER
                    fp  =  ck1(L,j) * eL1(L,j) + ck2(L,j) * em1(L,j) + cpb(L,j)
                    fm  =  ck1(L,j) * eL2(L,j) + ck2(L,j) * em2(L,j) + cmb(L,j)

                    fupbs(j)    = fupbs(j)    + fp * Beta_V(L)
                    fdownbs2(j) = fdownbs2(J) + fm * Beta_V(L)
                    fnetbs(j)   = fnetbs(j)   + fnet(L,j)

                    if (L .eq. NSOL) then
                        fdownbs(J) = (fupbs(j) - fnetbs(j))
                    endif
                end do
            END DO
          else if (opacity_method .eq. 'correk') then
            DO L       =  MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
                chan_idx = MODULO(L-1,8)+1
                ! stel_idx = (L - chan_idx)/8 + 1
                stel_idx = MODULO((L-chan_idx)/8, NWNO) + 1
                SOLNET  = SOLNET - FNET(L,NLAYER)
                ! This line shouldn't menanigful until reflection is added, but should have either a stel_spec of a wgts_corrk or both
                fp      = (ck1(L,1) * eL2(L,1) - ck2(L,1) * em2(L,1) + cp(L,1)) * WGTS_CORRK(chan_idx)!* STEL_SPEC(stel_idx)

                fsLu(L) = fsLu(L) + fp

                do j = 1, NLAYER
                    fp  =  ck1(L,j) * eL1(L,j) + ck2(L,j) * em1(L,j) + cpb(L,j)
                    fm  =  ck1(L,j) * eL2(L,j) + ck2(L,j) * em2(L,j) + cmb(L,j)

                    fupbs(j)    = fupbs(j)    + fp  * WGTS_CORRK(chan_idx) !* STEL_SPEC(stel_idx) !
                    fdownbs2(j) = fdownbs2(J) + fm  * WGTS_CORRK(chan_idx) !* STEL_SPEC(stel_idx) !
                    fnetbs(j)   = fnetbs(j)   + fnet(L,j)

                    if (L .eq. NSOL) then
                        fdownbs(J) = (fupbs(j) - fnetbs(j))
                    endif
                end do
            END DO
          end if
        !   write(*,*) "incident_starlight_fraction:", incident_starlight_fraction
        !   write(*,*) "fupbs:", fupbs
        !   write(*,*) "FDOWNBS:", fdownbs
        !   write(*,*) "FDOWNBS2:", fdownbs2
    ! 510       CONTINUE
          if ((opacity_method .eq. 'picket') .or. (opacity_method .eq. 'dogray')) then
            do  i = MAX(solar_calculation_indexer,MINWNOSTEL*8), nsoL              
                fsLd(i) = psol_aerad * incident_starlight_fraction * (Beta_V(i))
                alb_toa(i) = fsLu(i)/fsLd(i)
                tsLu = tsLu + fsLu(i)
                total_downwelling = total_downwelling + fsLd(i)
            END DO
          else if (opacity_method .eq. 'correk') then
            !psol_aerad is just solc_in here
            do  i = MAX(solar_calculation_indexer,MINWNOSTEL*8), nsoL
                chan_idx = MODULO(i-1,8)+1
                ! stel_idx = (i - chan_idx)/8 + 1
                stel_idx = MODULO((i-chan_idx)/8, NWNO) + 1
                fsLd(i) = psol_aerad * incident_starlight_fraction * WGTS_CORRK(chan_idx) * STEL_SPEC(stel_idx) !
                alb_toa(i) = fsLu(i)/fsLd(i)
                tsLu = tsLu + fsLu(i)
                total_downwelling = total_downwelling + fsLd(i)
            END DO
          end if





          alb_tomi = fupbs(1)/fdownbs(1)
          alb_toai = tsLu/total_downwelling

!         Load fluxes into interface common block

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

      do i = 1, nir
          firu(i) = 0.
      END DO


      IF (IR .NE. 0) THEN
          DO L        =  NSOL+1,NTOTAL
             firu(L-nsol ) = firu( L-nsol ) + directu(L,1)

             do j = 1, nlayer
                 fupbi(j)   = fupbi(j)   + (directu(L,j))
                 fdownbi(j) = fdownbi(j) + (direc(L,j))
                 fnetbi(j)  = fnetbi(j)  + (directu(L,j) - direc(L,j))
             END DO
          END DO

          do i = 1, nir
              tiru = tiru + firu(i)
          END DO

!         Load fluxes into interface common block

          do j = 1, nlayer
              fir_up_aerad(j)  = fupbi(nlayer+1-j)
              fir_dn_aerad(j)  = fdownbi(nlayer+1-j)
              fir_net_aerad(j) = fnetbi(nlayer+1-j)
              fsl_up_aerad(j)  = fupbs(nlayer+1-j)
              fsl_dn_aerad(j)  = fdownbs(nlayer+1-j)
              fsl_net_aerad(j) = fnetbs(nlayer+1-j)
          enddo
      ENDIF

C     RFLUXES  Array to hold fluxes at top and bottom of atmosphere
C     1st index - flux 1=SW, 2=LW
C     2nd index - Direction 1=DN, 2=UP
C     3rd index - Where 1=TOP, 2=SURFACE
      if (opacity_method .eq. 'correk') then
        do L = MAX(solar_calculation_indexer,MINWNOSTEL*8), NSOL
          chan_idx = MODULO(L-1,8)+1
          RFLUXES_aerad(1,1,1) = RFLUXES_aerad(1,1,1) + fsl_dn_aerad(NLAYER) * WGTS_CORRK(chan_idx)
          RFLUXES_aerad(1,1,2) = RFLUXES_aerad(1,1,2) + fsl_dn_aerad(1)/(1.0-ALBSW) * WGTS_CORRK(chan_idx)
          RFLUXES_aerad(1,2,1) = RFLUXES_aerad(1,2,1) + fsl_up_aerad(NLAYER) * WGTS_CORRK(chan_idx)
          RFLUXES_aerad(1,2,2) = RFLUXES_aerad(1,2,2) + RFLUXES_aerad(1,1,2)*ALBSW * WGTS_CORRK(chan_idx)
        end do
        do L=NSOL+1, NTOTAL
          chan_idx = MODULO(L-1,8)+1
          RFLUXES_aerad(2,1,1) = RFLUXES_aerad(2,1,1) + fir_dn_aerad(NLAYER) * WGTS_CORRK(chan_idx)
          RFLUXES_aerad(2,1,2) = RFLUXES_aerad(2,1,2) + fir_dn_aerad(1) * WGTS_CORRK(chan_idx)
          RFLUXES_aerad(2,2,1) = RFLUXES_aerad(2,2,1) + fir_up_aerad(NLAYER) * WGTS_CORRK(chan_idx)
          RFLUXES_aerad(2,2,2) = RFLUXES_aerad(2,2,2) + fir_up_aerad(1) * WGTS_CORRK(chan_idx)
        end do
          
      else if (NSOL .gt. 1) then ! picket-fence
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
      else ! double gray?
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

      subroutine bezier_interp(xi, yi, ni, x, y)
        implicit none

        integer, intent(in) :: ni
        real, dimension(ni), intent(in) :: xi, yi
        real x, y
        real xc, dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

        dx = xi(2) - xi(1)
        dx1 = xi(3) - xi(2)
        dy = yi(2) - yi(1)
        dy1 = yi(3) - yi(2)

        if (x > xi(1) .and. x < xi(2)) then
          w = dx1/(dx + dx1)
          yc = yi(2) - dx / 2.0 * (w*dy/dx + (1.0 - w)*dy1/dx1)
          t = (x - xi(1))/dx
          y = (1.0 - t)**2 * yi(1) + 2.0 *t*(1.0 - t)*yc + t**2*yi(2)
        else
          w = dx/(dx + dx1)
          yc = yi(2) + dx1 / 2.0 * (w*dy1/dx1 + (1.0 - w)*dy/dx)
          t = (x - xi(2))/(dx1)
          y = (1.0 - t)**2 * yi(2) + 2.0 *t*(1.0 - t)*yc + t**2*yi(3)
        end if
      end subroutine bezier_interp
