      subroutine radsub(iffirst,pr,p_pass,t_pass,qh2o_pass,
     &                  radheat,htlw,htsw,rfluxes,alat,alon,KOUNT,ITSPD,Beta_IR,Beta_V,
     &                  incident_starlight_fraction,TAURAY, TAUL, TAUGAS,TAUAER, solar_calculation_indexer,dpg,
     &                  ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &                  fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &                  pbar, dpgsub, pbarsub)
!     &   NPROB, SOL, RAYPERBAR, WEIGHT, GOL, WOL, TAUCONST,
!     &   WAVE, TT, Y3, PTEMPG, PTEMPT, G0, OPD,
!     &   PTEMP, uG0, uTAUL, W0, uW0,
!     &   uopd, U1S, U1I, TOON_AK, B1, B2, EE1,
!     &   EM1, EM2, EL1, EL2, GAMI, AF,
!     &   BF, EF, SFCS, B3, CK1, CK2,
!     &   CP, CPB, CM, CMB, DIRECT,
!     &   FNET, EE3,EL3, TMI, AS,
!     &   DF, DS, XK, DIREC, DIRECTU,
!     &   DINTENT, UINTENT, TMID, TMIU,
!     &   firu, fird, fsLu, fsLd, fsLn, alb_toa, fupbs, fdownbs,
!     &   fnetbs, fdownbs2, fupbi, fdownbi, fnetbi)


!     iffirst is just the indicator for numbering and runs the setup
!     deltaz--the layer thickness in meters
!     p_pass--the layer boundary pressures in pascal (NL+1)
!     t_pass--mid-layer temperatures in K and one bottom boundary temp
!     both p_ and t_ pass begin at the top and go down.

      include 'rcommons.h'
      PARAMETER(PI2=2.0*3.14159265359)
      integer iffirst
      real t_pass(NZ)
      real dpg(nl+1), pbar(nl+1)
      real dpgsub(2*nl+2), pbarsub(2*nl+2)
      real p_pass(NZ)
      real qh2o_pass(NZ)
      real radheat(NZ)
      real heats_aerad_tot(NZ), heati_aerad_tot(NZ), radheat_tot(NZ)
      real wave_pass(1)
      real cheats(NZ), cheati(NZ), cheat(NZ)
      real htlw(NZ), htsw(NZ)
      real PSOL,PSOL_aerad
      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY, TAUL, TAUGAS,TAUAER

      integer ifsetup
      real ibinm
      real rfluxes_aerad(2,2,2)
      real heati_aerad(NL+1)
      real heats_aerad(NL+1)
      real fsl_up_aerad(NL+1)
      real fsl_dn_aerad(NL+1)
      real fir_up_aerad(NL+1)
      real fir_dn_aerad(NL+1)
      real fir_net_aerad(NL+1)
      real fsl_net_aerad(NL+1)

      integer itime, ntime, solar_calculation_indexer

      ! Malsky add
      REAL AMU0, SOLC, DDAY, FORCE1DDAYS, DFAC, temporary_local_variable, ALON, ALAT, incident_starlight_fraction

      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL) :: Beta_V

      real rfluxes(2,2,2)
      REAL SSLON,SSLAT  ! ER:
      REAL DLENGTH  ! ER: half-length of solar day
      real PI2

!      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
!      REAL EPSILON, SOLNET, EMISIR, TPI, SQ3, SBK, AM, AVG, ALOS, SCDAY, RGAS
!      REAL U0, FDEGDAY, WOT, GOT, tslu, total_downwelling, alb_tot,tiru, alb_tomi, alb_toai
!      REAL HEATI(NL+1), HEATS(NL+1), HEAT(NL+1), GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5)
!      REAL NPROB(5), SOL(5), RAYPERBAR(5), WEIGHT(5), GOL(5,2*NL+2), WOL(5,2*NL+2),TAUCONST(5)
!      REAL WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2)
!      REAL PTEMP(5,2*NL+2), uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2)
!      REAL uopd(5,2*NL+2), U1S(5), U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(5,2*NL+2), EE1(5,2*NL+2)
!      REAL EM1(5,2*NL+2), EM2(5,2*NL+2), EL1(5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
!      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
!!      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2)
!      REAL FNET(5,2*NL+2), EE3(5,2*NL+2),EL3(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4)
!      REAL DF(5,4*NL+4), DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2)
!!      REAL DINTENT(5,3,2*NL+2), UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2)
 !     REAL firu(2), fird(2), fsLu(3), fsLd(3), fsLn(3), alb_toa(3), fupbs(NL+1), fdownbs(NL+1)
 !     REAL fnetbs(NL+1), fdownbs2(NL+1), fupbi(NL+1), fdownbi(NL+1), fnetbi(NL+1)

 582  FORMAT(I4,5(F12.3))

      ! Malsky what does this do???
      ibinm = ibinmin
      ifsetup = 0

      if( iffirst.eq. 1 ) THEN
          ifsetup = 1
      END IF

!     @ Keep an Eye on this, Mike
      if_diurnal = 0

      heats_aerad_tot = 0.
      heati_aerad_tot = 0.

!     @ The following lines of code are taken from cnikos and may require adjustment
C     ER modif for non-synchronous orbit

      IF (PORB.NE.0) THEN
         SSLON=(1./PORB-1.)*KOUNT*360./ITSPD
         SSLON=MOD(SSLON,360.)
      ELSE
         SSLON=0.  ! substellar longitude
      ENDIF
C ER modif for non-zero obliquity
      IF (OBLIQ.EQ.0) THEN
         SSLAT=0.  ! substellar latitude
         DLENGTH=PI/2.
      ELSE
         SSLAT=ASIN(SIN(OBLIQ*PI/180.)
     +        *SIN(PI2*KOUNT/ITSPD/PORB))*180./PI
         IF (SSLAT.GT.0) THEN
            IF (ALAT.GT.90.-SSLAT) THEN
               DLENGTH=PI
            ELSEIF (ALAT.LT.-90.+SSLAT) THEN
               DLENGTH=0.
            ELSE
               DLENGTH=ACOS(-1.*TAN(ALAT/360.*PI2)*TAN(SSLAT/360.*PI2))
            ENDIF
         ELSEIF (ALAT.LT.-90.-SSLAT) THEN
            DLENGTH=PI
         ELSEIF (ALAT.GT.90+SSLAT) THEN
            DLENGTH=0.
         ELSE
            DLENGTH=ACOS(-1.*TAN(ALAT/360.*PI2)*TAN(SSLAT/360.*PI2))
         ENDIF
      ENDIF

C Setup SW code                                                           
      IF (LBIN) THEN
        call BinaryFlux(SOLC,KOUNT,ITSPD)
        SOLC=SOLC*(1.0-TOAALB)
      ELSE
        SOLC=SOLC_IN*(1.0-TOAALB)
      ENDIF


C     globally averaged solar constant, vertical rays
      AMU0=1.0
      PSOL=SOLC/4.
      IF(.NOT.L1DZENITH) THEN
         DDAY=FORCE1DDAYS
         IF(DAY.GT.DDAY) THEN
            DFAC=MIN(1.0,(DAY - DDAY)/DDAY)
            IF(.NOT.LDIUR) THEN
               AMU0=(1.0-DFAC)*AMU0
     &              +DFAC*MAX(0.0,SIN(ALAT/360.*PI2)*SIN(SSLAT/360.*PI2)
     &                           +COS(ALAT/360.*PI2)*COS(SSLAT/360.*PI2)
     &                           *COS((ALON-SSLON)/360.*PI2))
               PSOL=(1.0-DFAC)*PSOL + DFAC*SOLC
            ELSE
               PSOL=(1.0-DFAC)*PSOL+DFAC*SOLC/PI*
     &              (SIN(ALAT/360.*PI2)*SIN(SSLAT/360.*PI2)*DLENGTH
     &              +COS(ALAT/360.*PI2)*COS(SSLAT/360.*PI2)*SIN(DLENGTH))
            ENDIF
         ENDIF
      ENDIF


      if ((AMU0.gt.0) .and. (AMU0.lt.1e-6)) THEN
          AMU0 = 0.0
      endif

      incident_starlight_fraction = MAX(0.0, AMU0)

      if (incident_starlight_fraction .lt. 1e-5) THEN
          solar_calculation_indexer = NSOLP + 1
      ELSE
          solar_calculation_indexer = 1
      END IF

      PSOL_aerad=PSOL
      ntime = 1

      if( if_diurnal.eq.1 ) ntime = 24

      do itime = 1, ntime
          t_pass(NLAYER) = t_pass(NLAYER-1)
          call setuprad_simple(Beta_V, Beta_IR, t_pass, pr, p_pass, incident_starlight_fraction,
     &                          TAURAY,TAUL,TAUGAS,TAUAER,
     &                          solar_calculation_indexer, DPG,
     &             ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &             fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &             pbar, dpgsub, pbarsub)
!     &   NPROB, SOL, RAYPERBAR, WEIGHT, GOL, WOL, TAUCONST,
!     &   WAVE, TT, Y3, PTEMPG, PTEMPT, G0, OPD,
!     &   PTEMP, uG0, uTAUL, W0, uW0,
!     &   uopd, U1S, U1I, TOON_AK, B1, B2, EE1,
!     &   EM1, EM2, EL1, EL2, GAMI, AF,
!     &   BF, EF, SFCS, B3, CK1, CK2,
!     &   CP, CPB, CM, CMB, DIRECT,
!     &   FNET, EE3,EL3, TMI, AS,
!     &   DF, DS, XK, DIREC, DIRECTU,
!     &   DINTENT, UINTENT, TMID, TMIU,
!     &   firu, fird, fsLu, fsLd, fsLn, alb_toa, fupbs, fdownbs,
!     &   fnetbs, fdownbs2, fupbi, fdownbi, fnetbi)



          call radtran(Beta_V, Beta_IR, incident_starlight_fraction,TAURAY,TAUL,TAUGAS,TAUAER,
     &                 solar_calculation_indexer, DPG, pr, t_pass, p_pass,
     &             ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &             fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &             pbar, dpgsub, pbarsub)
!     &   NPROB, SOL, RAYPERBAR, WEIGHT, GOL, WOL, TAUCONST,
!     &   WAVE, TT, Y3, PTEMPG, PTEMPT, G0, OPD,
!     &   PTEMP, uG0, uTAUL, W0, uW0,
!     &   uopd, U1S, U1I, TOON_AK, B1, B2, EE1,
!     &   EM1, EM2, EL1, EL2, GAMI, AF,
!     &   BF, EF, SFCS, B3, CK1, CK2,
!     &   CP, CPB, CM, CMB, DIRECT,
!     &   FNET, EE3,EL3, TMI, AS,
!     &   DF, DS, XK, DIREC, DIRECTU,
!     &   DINTENT, UINTENT, TMID, TMIU,
!     &   firu, fird, fsLu, fsLd, fsLn, alb_toa, fupbs, fdownbs,
!     &   fnetbs, fdownbs2, fupbi, fdownbi, fnetbi)


          cheats = 0.
          cheati = 0.
          cheat = 0.
          do iz = 1,NZ
              jz = NZ + 1 - iz
              radheat(iz) = heats_aerad(jz) + heati_aerad(jz)
              heats_aerad_tot(iz) = heats_aerad_tot(iz) + heats_aerad(jz) * SCDAY - cheats(iz)
              heati_aerad_tot(iz) = heati_aerad_tot(iz) + heati_aerad(jz) * SCDAY - cheati(iz)
              radheat_tot(iz) = radheat_tot(iz)+heats_aerad(jz) * SCDAY-cheats(iz)+ heati_aerad(jz)*SCDAY - cheati(iz)
          enddo
      enddo

      if( if_diurnal.eq. 1 ) then
        write(*,*)'if_diurnal ==1'
        write(*,*)'ntime',ntime
        heats_aerad_tot = heats_aerad_tot / ntime
        heati_aerad_tot = heati_aerad_tot / ntime
        radheat_tot = radheat_tot / ntime
      endif

      htlw    = heati_aerad_tot
      htsw    = heats_aerad_tot
      rfluxes = rfluxes_aerad

      return
      end

