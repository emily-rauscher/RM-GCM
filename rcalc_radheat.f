!***********************************************************************  
!*                         SUBROUTINE CALC_RADHEAT                     *  
!*********************************************************************** 
      SUBROUTINE CALC_RADHEAT(pr,t,p_pass,alat1,alon,htlw,htsw,
     $             DOY,cf,ic,rfluxes,swalb,kount,itspd, incident_starlight_fraction, TAURAY, TAUL, TAUGAS, TAUAER,
     &             solar_calculation_indexer, dpg,
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
     &  qrad,alb_tomi,alb_toai, num_layers, SLOPE, Y1, Y2, Y4, Y8, A1, A2, A3, A4, A5, A7, Y5,
     &  heats_aerad_tot, heati_aerad_tot, radheat_tot, radheat, cheati, cheats,
     &  dpe, Pl, Tl, pe,
     &  k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
     &  Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val, tau_IRe, tau_Ve,
     &  PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom)

!      use physical_constants

      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5),NPROB(5), SOL(5),RAYPERBAR(5),WEIGHT(5)
      REAL GOL(5,2*NL+2), WOL(5,2*NL+2), WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2), PTEMP(5,2*NL+2)
      REAL uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2), uopd(5,2*NL+2),  U1S( 5)
      REAL U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(5,2*NL+2), EE1( 5,2*NL+2), EM1(5,2*NL+2)
      REAL EM2(5,2*NL+2), EL1( 5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2), EE3(5,2*NL+2)
      REAL EL3(5,2*NL+2), FNET(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4), DF(5,4*NL+4)
      REAL DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2), DINTENT(5,3,2*NL+2)
      REAL UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(2),fird(2),fsLu(3), fsLd(3),fsLn(3),alb_toa(3), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai, SLOPE(5,2*NL+2)
      real heats_aerad_tot(NL+1), heati_aerad_tot(NL+1), radheat_tot(NL+1), cheati(NL+1), cheats(NL+1), radheat(NL+1)

      REAL, DIMENSION(5,3,2*NL+2) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(5,2*NL+2)   :: A1, A2, A3, A4, A5, A7, Y5


      REAL tau_IRe(2,NL+1), tau_Ve(3,NL+1)
      real, dimension(NL+1) :: dpe, Pl, Tl, pe
      real :: k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met
      real :: Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val

      REAL PI0_TEMP(5, NL+1, 13)
      REAL G0_TEMP(5, NL+1, 13)
      REAL tauaer_temp(5, NL+1, 13)
      INTEGER j1
      REAL denom




      real, parameter :: BK = 1.38054e-16
      real, parameter :: L      = 2.5e10
       
      REAL PR(NL+1),T(NL+1),Cpd,p_pass(NL+1)
      real dpg(NLAYER), pbar(NLAYER)
      real dpgsub(NDBL), pbarsub(NDBL)
      real, dimension(NL+1) :: z, htsw,htlw
      real, dimension(2,2,2) :: rfluxes
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL) :: Beta_V
      real, dimension(NIR+NSOL,NDBL) :: TAURAY, TAUL, TAUGAS, TAUAER

      real alat1, alon, incident_starlight_fraction

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

      tgrnd=t(NZ)
      rfluxes=fluxes
      iffirst = 1
      R_AIR = GASCON*10000.
      Cpd = GASCON/AKAP !1.004e+7
      PREF = P0
      GRAV = GA * 100.
      eps= Rd/Rv
      RdCp = Rd/Cpd


      call radsub(iffirst,pr,p_pass,t,radheat,htlw,htsw,rfluxes,alat1,alon,KOUNT,ITSPD,Beta_IR,Beta_V,
     &            incident_starlight_fraction, TAURAY, TAUL, TAUGAS, TAUAER,solar_calculation_indexer, DPG,
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
     &  qrad,alb_tomi,alb_toai, num_layers, SLOPE, Y1, Y2, Y4, Y8, A1, A2, A3, A4, A5, A7, Y5,
     &  heats_aerad_tot, heati_aerad_tot, radheat_tot, cheati, cheats,
     &  dpe, Pl, Tl, pe,
     &  k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
     &  Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val, tau_IRe, tau_Ve,
     &  PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom)

        write(*,*) rfluxes(2,2,2)




      iffirst = 0

      IF ((LFLUXDIAG).AND.(KOUNTP-KOUTP.LT.NTSTEP_IN)) THEN !-1)) THEN
         WRITE(63,*) 'LATITUDE, LONGITUDE:',ALAT1,ALON
         WRITE(63,2010)'1)P(Bars)',
     $        'FLUXES(Wm-2): 2)LW UP','3)LW DOWN',         
     $        '4)LW NET',
     $        '5)SW UP','6)SW DOWN',
     $        '7)SW NET'
 2010    FORMAT(1X,A9,1X,A21,4X,A9,4X,A8,4x,A7,4x,A9,4x,A8)                               

C     ER Modif: output pressures in bar instead of mbar
         DO ILAY=NLAYER,1,-1                                                  
          WRITE(63,2013),p_pass(NLAYER-ILAY+1)*1e-5,FIR_UP_AERAD(ILAY),
     $    FIR_DN_AERAD(ILAY),FIR_NET_AERAD(ILAY),
     $    FSL_UP_AERAD(ILAY),FSL_DN_AERAD(ILAY),FSL_NET_AERAD(ILAY)                         
     $
 2013       FORMAT(2X,F12.6,3X,E12.5,3X,E12.5,3X,E12.5,3X,E12.5         
     $             ,3X,E12.5,3X,E12.5)

      END DO
!
         WRITE(63,2023)'PRESSURE (Bars)','HEATING RATES: SW (K/DAY)'          
     $        ,'LW (K/DAY)','TOTAL (K/DAY)'                                                 
 2023    FORMAT(2X,A15,2X,A25,4x,A10,4x,A13)                                        
         DO LHT=NLAYER,1,-1                                                  
         WRITE(63,2020) PR(NLAYER-LHT+1)*1e-5,HTSW(LHT),HTLW(LHT),
     $          HTSW(LHT)+HTLW(LHT)
     $  
 2020       FORMAT(2X,F12.6,19X,E12.5,3X,E12.5,3x,E12.5)   
         END DO
         WRITE(63,*)

!  HERE WE WRITE TO FILE 62 ADDTIONAL RADIATIVE TRANSFER BY PRODUCTS
         WRITE(62,*)'LATITUDE, LONGITUDE:',ALAT1,ALON
         WRITE(62,*)'  Cosine of the incidencd angle, mu0:',U0
         write(62,*)' Top of atmosphere albedo: ',alb_toai
         write(62,*)'  Top of model albedo: ',alb_tomi
         write(62,*)'  SW Flux absorbed at bottom boundary (Wm-2): '
     &                 ,SOLNET
         WRITE(62,*)'  Upwelling LW Flux at top-of-atmopshere: '
     &                 ,tiru
         Write(62,*)'  Up- & Downward SW flux at TOA: '
         WRITE(62,*)   tslu,total_downwelling
         WRITE(62,2031)'1)P(Bars)',
     &        '2)TAULS (SW & LW)', 
     &        '3)CUMMULATIVE TAUS',
     &        '4)PI0s',
     &        '5)G0s',
     &        '6)DIRECT SW (d-scaled)',
     &        '7)4PI*|INTENSITIES|(W/M^2)'
 2031    FORMAT(3X,A9,6X,A17,9X,A18,12X,A6,11x,A5,7x,A20,3x,A26)            
         
          DO IL=1,NLAYER                                                  
          WRITE(62,2033),p_pass(IL)*1e-5,uTAUL(1,IL),
     $    uTAUL(2,IL),uOPD(1,IL),uOPD(2,IL),uW0(1,IL),uW0(2,IL),uG0(1,IL),
     $    uG0(2,IL),DIRECT(1,IL),TMI(1,IL),TMI(2,IL)                   
 2033       FORMAT(F12.6,2X,E12.5,1X,E12.5,2X,E12.5,1x,E12.5,3X,
     $             F7.4,1x,F7.4,2X,F7.4,1X,F7.4,2X,E12.5,2X,
     $             E12.5,1X,E12.5) 
          END DO
          WRITE(62,*) ''
      ENDIF
      end

