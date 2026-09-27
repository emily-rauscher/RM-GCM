      subroutine radsub(iffirst,pr,p_pass,t, radheat,htlw,htsw,alat1,alon,KOUNT,ITSPD,Beta_IR,Beta_V,
     &                  incident_starlight_fraction,TAURAY, TAUL, TAUGAS,TAUAER, solar_calculation_indexer,dpg,
     &                  ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &                  fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &                  pbar, dpgsub, pbarsub,
     &                  LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, EMISIR,
     &                  EPSILON, HEATI, HEATS, HEAT, SOLNET,
     &                  TPI, SQ3, SBK,AM, AVG, ALOS,
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
     &  heats_aerad_tot, heati_aerad_tot, radheat_tot, cheati, cheats,
     &  dpe, Pl, Tl, pe,
     &  k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
     &  Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val, tau_IRe, tau_Ve,
     &  PI0_TEMP, G0_TEMP, tauaer_temp,j1,denom, fluxes, k_IRl, k_Vl, tau_ray_temp)


!     iffirst is just the indicator for numbering and runs the setup
!     deltaz--the layer thickness in meters
!     p_pass--the layer boundary pressures in pascal (NL+1)
!     both p_ and t_ pass begin at the top and go down.

      use corrkmodule, only: NWNO, MINWNOSTEL
      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,kount
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NBATCH), RSFX(NBATCH),NPROB(NBATCH), SOL(NBATCH)
      REAL RAYPERBAR(NBATCH),WEIGHT(NBATCH)
      REAL GOL(NBATCH,2*NL+2), WOL(NBATCH,2*NL+2), WAVE(NTOTAL+1), TT(NL+1), Y3(NBATCH,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NBATCH), PTEMPT(NBATCH), G0(NBATCH,2*NL+2), OPD(NBATCH,2*NL+2), PTEMP(NBATCH,2*NL+2)
      REAL uG0(NBATCH,2*NL+2), uTAUL(NBATCH,2*NL+2), W0(NBATCH,2*NL+2), uW0(NBATCH,2*NL+2), uopd(NBATCH,2*NL+2),  U1S(NBATCH)
      REAL U1I(NBATCH), TOON_AK(NBATCH,2*NL+2), B1(NBATCH,2*NL+2), B2(  5,2*NL+2), EE1(NBATCH,2*NL+2), EM1(NBATCH,2*NL+2)
      REAL EM2(NBATCH,2*NL+2), EL1(NBATCH,2*NL+2), EL2(NBATCH,2*NL+2), GAMI(NBATCH,2*NL+2), AF(NBATCH,4*NL+4)
      REAL BF(NBATCH,4*NL+4), EF(NBATCH,4*NL+4), SFCS(NBATCH), B3(NBATCH,2*NL+2), CK1(NBATCH,2*NL+2), CK2(NBATCH,2*NL+2)
      REAL CP(NBATCH,2*NL+2), CPB(NBATCH,2*NL+2), CM(NBATCH,2*NL+2), CMB(NBATCH,2*NL+2), DIRECT(NBATCH,2*NL+2), EE3(NBATCH,2*NL+2)
      REAL EL3(NBATCH,2*NL+2), FNET(NBATCH,2*NL+2), TMI(NBATCH,2*NL+2), AS(NBATCH,4*NL+4), DF(NBATCH,4*NL+4)
      REAL DS(NBATCH,4*NL+4), XK(NBATCH,4*NL+4), DIREC(NBATCH,2*NL+2), DIRECTU(NBATCH,2*NL+2), DINTENT(NBATCH,3,2*NL+2)
      REAL UINTENT(NBATCH,3,2*NL+2), TMID(NBATCH,2*NL+2), TMIU(NBATCH,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai, SLOPE(NBATCH,2*NL+2)

      REAL, DIMENSION(NBATCH,3,2*NL+2) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(NBATCH,2*NL+2)   :: A1, A2, A3, A4, A5, A7, Y5

      real, dimension(NKGAUSS, NL+1) :: k_IRl, tau_ray_temp
      real, dimension(NKGAUSS, NL+1) :: k_Vl

      REAL tau_IRe(NKGAUSS,NL+1), tau_Ve(NKGAUSS,NL+1)
      real, dimension(NL+1) :: dpe, Pl, Tl, pe
      real :: k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met
      real :: Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL)  :: Beta_V

      REAL PI0_TEMP(NBATCH, NL+1, 13)
      REAL G0_TEMP(NBATCH, NL+1, 13)
      REAL tauaer_temp(NBATCH, NL+1, 13)
      INTEGER j1
      REAL DENOM

      PARAMETER(PI2=2.0*3.14159265359)
      integer iffirst

      REAL PR(NL+1),T(NL+1), p_pass(NL+1)

      real dpg(nl+1), pbar(nl+1)
      real dpgsub(2*nl+2), pbarsub(2*nl+2)
      real radheat(NZ)
      real heats_aerad_tot(NZ), heati_aerad_tot(NZ), radheat_tot(NZ)
      real wave_pass(1)
      real cheats(NZ), cheati(NZ)
      real htlw(NZ), htsw(NZ)
      real PSOL,PSOL_aerad
      real, dimension(NBATCH,2*NL+2) :: TAURAY, TAUL, TAUGAS,TAUAER

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

      integer solar_calculation_indexer
      integer, AUTOMATIC :: itime, ntime
      integer, AUTOMATIC :: iband, band_solar_calc_idx, nbands
      real, AUTOMATIC :: tiru_acc, tslu_acc, total_downwelling_acc
      real, AUTOMATIC :: fir_up_acc(NL+1), fir_dn_acc(NL+1), fir_net_acc(NL+1)
      real, AUTOMATIC :: fsl_up_acc(NL+1), fsl_dn_acc(NL+1), fsl_net_acc(NL+1)

      ! Malsky add
      REAL AMU0, SOLC, DDAY, FORCE1DDAYS, DFAC, temporary_local_variable, ALON, ALAT1, incident_starlight_fraction

      real fluxes(2,2,2)
      REAL SSLON,SSLAT  ! ER:
      REAL DLENGTH  ! ER: half-length of solar day
      real PI2
 582  FORMAT(I4,5(F12.3))

      ! ibinmin was always 0 (BSS-init static); ibinm was removed from args, so assign 0 directly
      ibinm = 0
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
            IF (alat1.GT.90.-SSLAT) THEN
               DLENGTH=PI
            ELSEIF (alat1.LT.-90.+SSLAT) THEN
               DLENGTH=0.
            ELSE
               DLENGTH=ACOS(-1.*TAN(alat1/360.*PI2)*TAN(SSLAT/360.*PI2))
            ENDIF
         ELSEIF (alat1.LT.-90.-SSLAT) THEN
            DLENGTH=PI
         ELSEIF (alat1.GT.90+SSLAT) THEN
            DLENGTH=0.
         ELSE
            DLENGTH=ACOS(-1.*TAN(alat1/360.*PI2)*TAN(SSLAT/360.*PI2))
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
      ! multiply by sqrt3 here, then divide by it later either explicitly if L1DZENITH=F or implicitly 
      ! if L1DZENITH=T by setting the incident starlight fraction to 1/sqrt(3). This is my hacky way of making L1DZenith=T do a planet-average profile
      ! It's weird and complicated because the cosine of the zenith angle and the fraction of flux that is incident on a given column are identical in 3-D,
      ! but not if you want to take a 1-D average (see, e.g., Guillot+2010, Parmentier+Guillot 2014).
      PSOL=SOLC/4. * SQRT(3.0) 
      IF(.NOT.L1DZENITH) THEN
         DDAY=FORCE1DDAYS
         IF(DAY.GT.DDAY) THEN
            IF(DDAY.GT.0.0) THEN
               DFAC=MIN(1.0,(DAY - DDAY)/DDAY)
            ELSE
               DFAC=1.0
            ENDIF
            IF(.NOT.LDIUR) THEN
               AMU0=(1.0-DFAC)*AMU0
     &              +DFAC*MAX(0.0,SIN(alat1/360.*PI2)*SIN(SSLAT/360.*PI2)
     &                           +COS(alat1/360.*PI2)*COS(SSLAT/360.*PI2)
     &                           *COS((ALON-SSLON)/360.*PI2))
               PSOL=(1.0-DFAC)*PSOL/SQRT(3.0) + DFAC*SOLC
            ELSE
               PSOL=(1.0-DFAC)*PSOL/SQRT(3.0)+DFAC*SOLC/PI*
     &              (SIN(alat1/360.*PI2)*SIN(SSLAT/360.*PI2)*DLENGTH
     &              +COS(alat1/360.*PI2)*COS(SSLAT/360.*PI2)*SIN(DLENGTH))
            ENDIF
         ENDIF
        if ((AMU0.gt.0) .and. (AMU0.lt.1e-6)) THEN
          AMU0 = 0.0
        endif
        incident_starlight_fraction = MAX(0.0, AMU0)

      ELSE
         incident_starlight_fraction = 1.0 / SQRT(3.0)
         ! AMU0 = 1.0 / SQRT(3.0)
      ENDIF


      ! if ((AMU0.gt.0) .and. (AMU0.lt.1e-6)) THEN
      !     AMU0 = 0.0
      ! endif

      ! incident_starlight_fraction = MAX(0.0, AMU0)

      if (incident_starlight_fraction .lt. 1e-10) THEN
          solar_calculation_indexer = NSOL + 1
      ELSE
          solar_calculation_indexer = 1
      END IF

      PSOL_aerad=PSOL
      ntime = 1

      if( if_diurnal.eq.1 ) ntime = 24

      rfluxes_aerad = 0.

      do itime = 1, ntime
          t(NLAYER) = t(NLAYER-1)

          cheats = 0.
          cheati = 0.
          radheat = 0.
          radheat_tot = 0.
          tiru_acc = 0.
          tslu_acc = 0.
          total_downwelling_acc = 0.
          fir_up_acc = 0.
          fir_dn_acc = 0.
          fir_net_acc = 0.
          fsl_up_acc = 0.
          fsl_dn_acc = 0.
          fsl_net_acc = 0.

          IF (opacity_method .EQ. 'correk') THEN
              nbands = NWNO
          ELSE
              nbands = 1
          END IF

          DO iband = 1, nbands
              if (opacity_method .NE. 'correk') then
                  ! picket/dogray have no wavenumber-band/g-point structure to
                  ! batch over (unlike correk's NWNO bands x NKGAUSS g-points),
                  ! so just pass the indexer through as-is for the single pass.
                  band_solar_calc_idx = solar_calculation_indexer
              else if (solar_calculation_indexer .gt. NKGAUSS) then
                  band_solar_calc_idx = NKGAUSS + 1
              else if (iband .lt. MINWNOSTEL) then
                  band_solar_calc_idx = NKGAUSS + 1
              else
                  band_solar_calc_idx = 1
              end if

              call setuprad_simple(Beta_V((iband-1)*NKGAUSS+1), Beta_IR((iband-1)*NKGAUSS+1),
     &  t, pr, p_pass, incident_starlight_fraction,
     &  TAURAY,TAUL,TAUGAS,TAUAER,
     &  band_solar_calc_idx, DPG,
     &  ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &  fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &  pbar, dpgsub, pbarsub,
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
     &  tiru,firu((iband-1)*NKGAUSS+1),fird((iband-1)*NKGAUSS+1),
     &  fsLu((iband-1)*NKGAUSS+1),fsLd((iband-1)*NKGAUSS+1),
     &  fsLn((iband-1)*NKGAUSS+1),alb_toa((iband-1)*NKGAUSS+1),fupbs,
     &  fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &  qrad,alb_tomi,alb_toai, num_layers, SLOPE,
     &  dpe, Pl, Tl, pe,
     &  k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
     &  Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val,
     &  tau_IRe, tau_Ve, k_IRl, k_Vl, tau_ray_temp, iband)

              call radtran(Beta_V((iband-1)*NKGAUSS+1),
     &                 Beta_IR((iband-1)*NKGAUSS+1), incident_starlight_fraction,
     &                 TAURAY,TAUL,TAUGAS,TAUAER,
     &                 band_solar_calc_idx, DPG, pr, t, p_pass,
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
     &  tiru,firu((iband-1)*NKGAUSS+1),fird((iband-1)*NKGAUSS+1),
     &  fsLu((iband-1)*NKGAUSS+1),fsLd((iband-1)*NKGAUSS+1),
     &  fsLn((iband-1)*NKGAUSS+1),alb_toa((iband-1)*NKGAUSS+1),fupbs,
     &  fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &  qrad,alb_tomi,alb_toai, num_layers, SLOPE,
     &  Y1, Y2, Y4, Y8, A1, A2, A3, A4, A5, A7, Y5,
     &  PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom,kount,itspd,iband)

              tiru_acc = tiru_acc + tiru
              tslu_acc = tslu_acc + tslu
              total_downwelling_acc = total_downwelling_acc + total_downwelling

              fir_up_acc  = fir_up_acc  + fir_up_aerad
              fir_dn_acc  = fir_dn_acc  + fir_dn_aerad
              fir_net_acc = fir_net_acc + fir_net_aerad
              fsl_up_acc  = fsl_up_acc  + fsl_up_aerad
              fsl_dn_acc  = fsl_dn_acc  + fsl_dn_aerad
              fsl_net_acc = fsl_net_acc + fsl_net_aerad

              do iz = 1,NZ
                  jz = NZ + 1 - iz
                  radheat(iz) = radheat(iz) + heats_aerad(jz) + heati_aerad(jz)
                  heats_aerad_tot(iz) = heats_aerad_tot(iz) + heats_aerad(jz)*SCDAY
                  heati_aerad_tot(iz) = heati_aerad_tot(iz) + heati_aerad(jz)*SCDAY
                  radheat_tot(iz) = radheat_tot(iz) + heats_aerad(jz)*SCDAY
     &                                              + heati_aerad(jz)*SCDAY
              end do
          END DO ! iband

          tiru = tiru_acc
          tslu = tslu_acc
          total_downwelling = total_downwelling_acc
          fir_up_aerad  = fir_up_acc
          fir_dn_aerad  = fir_dn_acc
          fir_net_aerad = fir_net_acc
          fsl_up_aerad  = fsl_up_acc
          fsl_dn_aerad  = fsl_dn_acc
          fsl_net_aerad = fsl_net_acc
          if (total_downwelling_acc .gt. 0.) then
              alb_toai = tslu_acc / total_downwelling_acc
          end if
          if (fsl_dn_acc(NL) .gt. 0.) then
              alb_tomi = fsl_up_acc(NL) / fsl_dn_acc(NL)
          end if
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
      fluxes  = rfluxes_aerad

      !write(*,*) ALAT1, ',', ALON, ',', TT(1), ','
      !if (ALAT1 .lt. -19 .and. ALAT1 .gt. -20 .and. ALON .eq. 315) THEN
      !    write(*,*) 'the model has done every point'
      !    stop
      !END IF

      return
      end

