C**********************************************************
C             SUBROUTINE RADIATION
C**********************************************************
      SUBROUTINE RADIATION(TROPHT,IH)

C     RADIATION SCHEME DERIVED FROM PREVIOUS CMORC.F AND THE
C     TOON CODES (TOON ET AL 1989). THE SCHEME IS CURRENTLY DOUBLE GRAY
C     AND APPLIES THE TWO-STREAM APPROXIMATION WITH QUADRATURE IN THE
C     VISIBLE AND HEMISPHERIC MEAN IN THE INFRARED.

C     It passes the pressure of the full sigma levels and the surface
C     to the Radiation scheme temperatures from TG and TSTAR
C
C     Determines model resolution
      use omp_lib
      include 'params.i'

C     Sets basic constants, especially those needed for array dimensions

      PARAMETER(MH=2,PI=3.14159265359,PI2=2.0*PI
     +,NNP=NN+1,MGPP=MG+2,JGP=JG+1,JGG=JG*NHEM,JGGP=JGG+1,MJP=NWJ2+NWJ2
     +,NLM=NL-1,NLP=NL+1,NLPP=NL+2,NLA=NL+3,NLB=NL+4,NL2=NL*NL
     +,IDA=(MG+MG+MG)/2+1,IDB=NWJ2*NL,IDC=IDB+IDB,IDD=MGPP*NL
     +,IDE=NL2*NN,IDF=NCRAY*(MG+1),IDG=JG*NL,IDH=JG*MG
     +,IDI=NNP/2,IDJ=IDI*IDI,IDK=NL*IDI,IDL=MGPP/2,IDM=NNP/2,IDN=IDM*NL
     +,NWW=1+(MM-1)/MOCT)

      PARAMETER(IGA=NWJ2*NHEM,IGB=IDB*NHEM,IGC=MGPP*NHEM,IGD=IDD*NHEM
     +,IGG=IDG*NHEM,IGL=IDL*NHEM,IGM=IDM*NHEM,IGN=IDN*NHEM
     +,IGO=IGA+IGA,IGP=IGB+IGB,NFTWG=(5+NTRAC)*NL+3
     +,NFTGW=(6+3*NTRAC)*NL+2,NFTGD=(3+NTRAC)*NL,NLTR=NL*NTRAC)

C     Number of 2D (surface) output fields. This value is
C     Doubled due to averaged and instantaneous fields.
      PARAMETER (N2DFLD=21,NGRPAD=N2DFLD*2*IGC)

C     Basic planetary parameters for run plus information about
C     vertical grid structure

C     Note that RD and GASCON are identical and CPD is set from RD,AKAP.
      COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)
     +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)
     +              ,TKP(NL),C(NL2),SQH(NNP)
     +              ,MF,MFP,JZF,NF
     +              ,AKAP,GA,GASCON,RADEA,WW,PFAC,EZ,AIOCT
     +              ,RD,RV,CPD,CLATNT
     +              ,P0,LRSTRT,LSHORT,LTVEC,LSTRETCH
     +              ,LFLUX
     +              ,LBALAN,LRESTIJ
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT
     +              ,LNOISE,NFP

      COMPLEX EZ,AIOCT

      LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ
     +       ,LFLUX,LNOISE
     +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT

C     Legendre polynomials and information about gaussian latitudes

      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC

C     Array ordering in GRIDP must correspond to that in SPECTR.
C     Real arrays: multi-level arrays are 2-dimensional.
C     the variables have been renamed to coincide with
C     variable names in bgcm5 DGRMLT

      COMMON/GRIDP/ CHIG(IGC,NL),SFG(IGC,NL),UG(IGC,NL),VG(IGC,NL)
     :              ,TTVD(IGC,NL),QTVD(IGC,NL),TG(IGC,NL)
     :              ,TRAG(IGC,NL,NTRAC)
     :              ,PLG(IGC),TYBL(IGC),TXBL(IGC)
     :              ,SPG(IGC),VPG(IGC),TTRD(IGC,NL)
     :              ,TNLG(IGC,NL),TRANLG(IGC,NL,NTRAC),UNLG(IGC,NL)
     :              ,VNLG(IGC,NL),TTLR(IGC,NL),UTRAG(IGC,NL,NTRAC)
     :              ,TTCR(IGC,NL),VTRAG(IGC,NL,NTRAC)
     :              ,UTVD(IGC,NL),VTVD(IGC,NL)
     :         ,ASSBL(IGC),ASHBL(IGC),ASLBL(IGC),ARRCR(IGC),ARRLR(IGC)
     :         ,arflux(igc,6),asfld(igc,6),acld(igc,4)
     :         ,SSBL(IGC),SHBL(IGC),SLBL(IGC),RRCR(IGC),RRLR(IGC)
     :         ,rflux(igc,6),sfld(igc,6),cld(igc,4)

       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN

       LOGICAL LPLOTMAP

C     Constant arrays and variables associated with time and vertical
C     differencing. Also counters.

      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)
     +              ,BEGDOY,DOY
C
      COMMON/PHYS/  CCR,RCON,DTBUOY,TSLA,TSLB,TSLC,TSLD,CUT1,CUT2
     :              ,TSTAR(IGC,JG),QSTAR(IGC,JG),FRAD(JG,NHEM)
     :              ,TSTARO(IGC,JG),TDEEPO(IGC,JG),smstar(igc,jg)
     :              ,tdeep(igc,jg),hsnow(igc,jg),sqstar(igc,jg)
     :              ,SALB(IGC,JG),SBAL(IGC,JG),BLCD(IGC)
     :              ,SVEGE(IGC,JG),CD,DRAG,BLVAD,BLA,BLRH,BLVB(IGC)
     :              ,AKVV,AKTV,AKQV,ESCONA,ESCONB,EPSIQ,CTQ,CCC
     :              ,ctqi,sdsn,shcs,shcsp,shcsn,skse,sksn,slhf,sd1,sd2,sdw
     :              ,ssmc,sdsnd,sasnow,saice,shsstar,shsmax
     :              ,LOC,LNOICE,LOLDBL,LCOND,LNNSK
     :              ,NLCR,CURHM,AKTC,AKQC,CUBMT,CBADJT,CBADJP
     :              ,SKAP(NL),SK(NLM),FWS(NL),CLR(NL),FB(NLM)
     :              ,TTDC(NL),QTDC(NL),TTMC(NL),QTMC(NL),TC(NL),QC(NL)
     :              ,CTCR(NL,NHEM),CTLR(NL,NHEM)
     :              ,LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ
     :              ,LSL,NAVRD,NAVWT,DELT2C,SHCO,SHCI,ITSLL,ITSLO,NCUTOP

      LOGICAL LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ,LSL,LOC,LNOICE,LOLDBL,LCOND,LNNSK

      COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM(3), AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS(3),with_TiO_and_VO, picket_fence_optical_depths

      LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF,RAYSCAT,AEROSOLS, picket_fence_optical_depths

      CHARACTER(30) :: AEROSOLMODEL

      REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),TCON(NL+1)

      LOGICAL DELTASCALE,HAZES,PICKET_FENCE_CLOUDS

      COMMON/CLOUDY/ AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &               CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &               ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &               MAXTAU,MAXTAULOC,TCON,AERSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS

      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC
     +              ,RNTAPO
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG
     +              ,LSPO(NL),LGPO(NL)
     $              ,LSHIST,LMINIH

      LOGICAL LSHIST,LMINIH

      LOGICAL LSPO,LGPO

C     Setup moisture variables by equivilencing them to
C     Tracer No. 1

      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1)),(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))
C
      COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)
     :     ,SNET(IGC,JG),RRFLUX(IGC,JG,6)
     :     ,TTSW(IGC,NL),TTLW(IGC,NL)

      COMMON/GSG/GSG(IGC,JG)

      REAL htnet

      COMMON /RADHT/ HTNET(NHEM,JG,MG,NL)

      REAL TAVE(IGP)

      REAL PR(NL+1),T(NL+1),p_pass(nl+1),htlw(nl+1),htsw(nl+1)
      real dpg(nl+1), pbar(nl+1)
      real dpgsub(2*nl+2), pbarsub(2*nl+2)

      real, dimension(5,2*NL+2) :: TAURAY, TAUL, TAUGAS, TAUAER

      ! Malsky is adding these
      integer solar_calculation_indexer, num_layers, malsky_test

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

      real PRB2T(NL+1),adum

      integer ifirst                ! If =1, first time reading o3
                                    ! and h2o (2 months' worth).
      real amfrac                   ! fraction through month
      integer ichange               ! =1 when in process of month change
      integer ifirstcol             ! =1 first time through column
                                    ! calculation (open new file).
      real p0
      real ps                       ! sfc pressure (used in
                                    ! interpolation from climatology
                                    ! to model).
      integer im                    ! Pointer for array plg (for getting sfc pressure).

C     Array to hold fluxes at top and bottom of atmosphere
C     1st index - flux 1=SW, 2=LW
C     2nd index - Direction 1=DN, 2=UP
C     3rd index - Where 1=TOP, 2=SURFACE

      real fluxes(2,2,2)

      real incident_starlight_fraction

c     The following for parallel testing --MTR
      integer TID, NTHREADS

      double precision test_wctime

      save                          ! Want to keep things like dcompl.

      DATA IFIRST/1/
      data ifirstcol/1/

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NL+1), HEATS(NL+1), HEAT(NL+1), SOLNET
      REAL TPI, SQ3, SBK, AM, AVG, ALOS

      ! I just put a huge chunk of these in
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
      real heats_aerad_tot(NL+1), heati_aerad_tot(NL+1), radheat_tot(NL+1), cheati(NL+1), cheats(NL+1), radheat(NL+1)

      REAL, DIMENSION(5,3,2*NL+2) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(5,2*NL+2)   :: A1, A2, A3, A4, A5, A7, Y5

      real, dimension(2, NL+1) :: k_IRl
      real, dimension(3, NL+1) :: k_Vl

      ! For the new picket fence stuff
      REAL tau_IRe(2,NL+1), tau_Ve(3,NL+1)
      real, dimension(2)  :: Beta_IR
      real, dimension(3)  :: Beta_V

      real, dimension(NL+1) :: dpe, Pl, Tl, pe
      real :: k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met
      real :: Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val


      REAL PI0_TEMP(5, NL+1, 13)
      REAL G0_TEMP(5, NL+1, 13)
      REAL tauaer_temp(5, NL+1, 13)
      INTEGER j1
      real denom



      COMMON /CLOUD_PROPERTIES/ TCONDS, KE_OPPR, PI0_OPPR, G0_OPPR,
     &                              DENSITY, FMOLW,
     &                              CORFACT,
     &                              input_particle_size_array_in_meters,
     &                              input_temperature_array,
     &                              particle_size_vs_layer_array_in_meters,
     &                              input_pressure_array_cgs,
     &                              HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg,
     &                              HAZE_PlanckMean_tau_per_bar,HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg,
     &                              HAZE_wav_tau_per_bar,HAZE_wav_pi0, HAZE_wav_gg,
     &                              haze_pressure_array_pascals, HAZE_WAV_GRID, CLOUD_WAV_GRID, exp_92_lnsig2_pi



      DATA EPSILON / 1d-6  /
      DATA SBK    / 5.6697E-8    /
      DATA AVG    /6.02252E+23/
      DATA ALOS   / 2.68719E19   /

      DATA GANGLE  /0.2123405382, 0.5905331356,0.9114120405/
      DATA GRATIO  /0.4679139346, 0.3607615730, 0.1713244924/
      DATA GWEIGHT /0.0698269799, 0.2292411064,0.2009319137/

      DATA RGAS   / 8.31430E+07  /
      DATA SCDAY  / 86400.0      /

      num_layers = NL

      RHSCL=288.0*GASCON/GA

      ! factor to non-dimensionalise heating rates
      CHRF=86400.*WW*CT

c     Skipping part set up here.
c     Radiation scheme only called every nskip longitudes
c     nskip must divide exactly into mg for longitude.
c     ntstep is the number of timesteps to skip.

      ntstep=NTSTEP_IN
      nskip=NSKIP_IN

      IOFM=0

      DO 800 ihem=1,nhem
        IF (mod(kount,ntstep) .eq. 0) THEN
          ilast=0

          ! schedule(guided), default(none),

          ! Do all the parallel stuff here
          !$OMP PARALLEL DO private(test_wctime,
!     &    im,idocalc, incident_starlight_fraction, RAYSCAT, solar_calculation_indexer, qrad, alb_toai,
!     &    dpe, Pl, Tl, pe, k_IRl, k_Vl,
!     &    PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom,
!     &    k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
!     &    Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val, tau_IRe, tau_Ve, Beta_IR, Beta_V,
!     &    heats_aerad_tot, heati_aerad_tot, radheat_tot, radheat, cheati, cheats,
!     &    EF, SFCS,
!     &    imp,PR,T,
!     &    TT, Y1, Y2, Y4, Y8, A1, A2, A3, A4, A5, A7, Y5,
!     &    imm,alat1,cf,ic,SWALB,alon,htlw, fluxes, GA,
!     &    htsw,HTNETO,a,b,
!     &    PRB2T, AEROPROF, ALBSW, AEROSOLS, AEROSOLMODEL,  IH,
!     &    EMISIR,
!     &    HEATI, HEATS, HEAT,
!     &    SOLNET,
!     &    LLA, LLS,
!     &    TPI,
!     &    AM,
!     &    EMIS, RSFX,NPROB,
!     &    SOL,RAYPERBAR, RAYPERBARCONS, WEIGHT,
!     &    GOL,
!     &    WOL,
!     &    iblackbody_above,
!     &    WAVE,
!     &    TAUGAS, TAURAY, TAUAER, TAUAEROSOL, TAUL, TAUA, uTAUL,
!     &    Y3,
!     &    U0,  ISL, IR, IRS,
!     &    WOT, GOT,
!     &    PTEMPG, PTEMPT,
!     &    G0, OPD,
!     &    PTEMP,
!     &    uG0,
!     &    W0,uW0,
!     &    uopd,
!     &    U1S,U1I,
!     &    TOON_AK,
!     &    B1,
!     &    B2, EE1,
!     &    EM1,
!     &    EM2, EL1,
!     &    EL2,    GAMI,
!     &    AF, BF, EF
!     &    SFCS,
!     &    B3,   CK1,
!     &    CK2,   CP,
!     &    CPB,   CM,
!     &    CMB,   DIRECT,
!     &    EE3,   EL3,
!     &    FNET,   TMI,
!     &    AS,     DF,
!     &    DS,     XK,
!     &    DIREC, DIRECTU, SLOPE,
!     &    DINTENT,
!     &    UINTENT,
!     &    TMID,TMIU,
!     &    tslu,total_downwelling,alb_tot,tiru,firu,fird,fsLu,
!     &    fsLd,fsLn,alb_toa,
!     &    fupbs,fdownbs,fnetbs,fdownbs2,
!     &    fupbi,fdownbi,fnetbi,
!     &    alb_tomi,alb_toai,
!     &    TCONDS, KE_OPPR, PI0_OPPR, G0_OPPR,
!     &    DENSITY, FMOLW, MOLEF,
!     &    CORFACT,
!     &    input_particle_size_array_in_meters,
!     &    input_temperature_array,
!     &    particle_size_vs_layer_array_in_meters,
!     &    input_pressure_array_cgs,
!     &    ifsetup, ibinm,
!     &    rfluxes_aerad,
!     &    psol_aerad,
!     &    heati_aerad, heats_aerad,
!     &    fsl_up_aerad, fsl_dn_aerad,
!     &    fir_up_aerad, fir_dn_aerad,
!     &    fir_net_aerad,fsl_net_aerad,
!     &    p_pass, dpg, pbar, dpgsub, pbarsub,
!     &    firstprivate(ilast),
!     &    lastprivate(ilast))

          DO i=1,mg

            im=i+iofm
            idocalc=0

            IF ((i.eq.1).or.(i-ilast.ge.nskip)) then
              idocalc=1
            ELSE

            IF (.FALSE.) write(*,*) ! This shouldn't be here, but necessary for the parallel work

            IF (LNNSK) THEN
              imp=im+1

              IF (imp.gt.(mg+iofm)) imp=1+iofm
                imm=im-1

                IF (((gsg(im,jh).gt.0.).and.(gsg(imp,jh).eq.0.)).or.
     $             ((gsg(im,jh).eq.0.).and.(gsg(imp,jh).gt.0.)).or.
     $             ((gsg(im,jh).gt.0.).and.(gsg(imm,jh).eq.0.)).or.
     $             ((gsg(im,jh).eq.0.).and.(gsg(imm,jh).gt.0.))) THEN
                  idocalc=1
                ENDIF
              ENDIF
            ENDIF

            IF (idocalc.eq.1) then
              DO LD=1,NL    ! Start of loop over column.
                L=NL-LD+2  ! Reverse index (Morc goes bottom up).
                PR(LD)=SIGMA(LD)*PLG(im)*P0 ! Pressure
                PRB2T(L)=PR(LD)
                T(LD)=TG(im,ld)*CT ! Temperature
                AEROPROF(LD)=0.0
              ENDDO

              AEROPROF(NL+1)=0.0
              PRB2T(1)=PLG(im)*P0
              PR(NL+1)=PLG(im)*P0
              T(NL+1)=((FBASEFLUX+rrflux(IM,JH,1))/5.6704e-8)**0.25

              alat1=alat(JH)*REAL(-(ihem*2.)+3)
              IF ((LFLUXDIAG).AND.(KOUNTP-KOUTP.LT.NTSTEP_IN)) THEN
                IF(JH.EQ.1.AND.IHEM.EQ.1.AND.I.EQ.1) THEN
                  REWIND(63) !! Rewind file for fluxes in nikosrad
                  REWIND(62) ! rwnd file for ancillary RT results

                  IF (PORB.NE.0) THEN
                    SSLON=(1./PORB-1.)*KOUNT*360./ITSPD
                    SSLON=MOD(SSLON,360.)
                  ELSE
                    SSLON=0.
                  ENDIF

                  SSLAT=ASIN(SIN(OBLIQ*PI/180.)*SIN(PI2*KOUNT/ITSPD/PORB))*180./PI
                  WRITE(63,2021) DAY,SSLON,SSLAT
                  WRITE(62,2021) DAY,SSLON,SSLAT

 2021             FORMAT('DAY:',F7.2,', SUBSTELLAR LON,LAT:',2F7.2)

                  WRITE(63,*)
                  WRITE(62,*)''
                ENDIF
              ENDIF

              SWALB=ALBSW

!             PR AND T ARE THE TEMPERATURE AND PRESSURE AT THE SIGMA LEVELS
!             AND BOTTOM BOUNDARY, AS USED BY THE DYNAMICAL CODE.
!             TO COMPUTE HEATING RATES AT THESE CENTERS, WE NEED TO DEFINE
!             LAYER EDGES AT WHICH FLUXES ARE COMPUTED, p_pass.

              DO LD    = 1,NL-1
                p_pass(LD+1)=(pr(LD)+pr(LD+1))/2.
              ENDDO

              p_pass(NL+1)=PR(NL+1)
              p_pass(1)=pr(1)*0.5

              alon=REAL(i-1)/REAL(mg)*360.0

!             PR in pascals for layer boundaries (NL+1), T in Kelvin for layer
!             centers + one layer for the bottom boundary. The top is n=1, the
!             bottom is n=NL+1

!             Extract a single column from the array AER4LAT(NLEV,LON,HEM)
              IF((AEROSOLS).AND.(AEROSOLMODEL.NE.'Global')) THEN
                DO  LD=1,NL +1
                  AEROPROF(LD)=TAUAEROSOL(LD,i,ihem,ih)
                ENDDO
              ENDIF

              rfluxes_aerad = 0.
              fluxes        = 0.
              fsl_dn_aerad  = 0.
              fir_up_aerad  = 0.

               Y3  = 0.
               Y1  = 0.
               Y5  = 0.
               Y3  = 0.
               Y4  = 0.

               k_IRl = 0
               k_Vl  = 0

              call calc_radheat(pr,t,p_pass,alat1,alon,htlw,htsw,
     &                          DOY,cf,ic,fluxes,swalb,kount,itspd,
     &                          incident_starlight_fraction,TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, dpg,
     &           ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &           fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &           pbar, dpgsub, pbarsub,
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
     &  PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom, Beta_IR, Beta_V, k_IRl, k_Vl)

              pr=prb2t

              PNET(IM,JH)=fluxes(1,1,1)-fluxes(1,2,1)+fluxes(2,1,1)-fluxes(2,2,1)
              SNET(IM,JH)=fluxes(1,1,2)-fluxes(1,2,2)+fluxes(2,1,2)-fluxes(2,2,2)

              rrflux(im,jh,1)=fluxes(1,1,2)
              rrflux(im,jh,2)=fluxes(1,2,2)
              rrflux(im,jh,3)=fluxes(2,1,2)
              rrflux(im,jh,4)=fluxes(2,2,2)
              rrflux(im,jh,5)=fluxes(1,1,1)-fluxes(1,2,1)
              rrflux(im,jh,6)=fluxes(2,2,1)

c             bottom heating rate is zero in morecret
              DO l=nl,1,-1
                LD=NL+1-L
                IM=I+IOFM
                HTNETO=HTNET(IHem,JH,I,LD)
                htnet(ihem,jh,i,ld)=(htlw(l+1)+htsw(l+1))
                TTRD(IM,LD)=(HTNETO+HTNET(IHEM,JH,I,LD))/(CHRF*2.0)
              ENDDO


              ilast=i
!             end of conditional execution of morcrette code

            ENDIF
          enddo
          !$OMP END PARALLEL DO

          IF (nskip.ne.0) then
             write(*,*),'CANNOT SKIP LONGITUDES IN PARALLEL!! ABORT'
             write(*,*),'Please set nskip=0 in fort.7'
             STOP
          ELSE
              ilast=mg
          ENDIF

          IF (ilast.ne.mg) then
            DO j=ilast+1,mg
              a=REAL(j-ilast)/REAL(mg+1-ilast)
              b=1.-a
              im=j+iofm
              DO l=nl,1,-1
                ld=nl+1-l
                HTNETO=HTNET(IHEM,JH,J,LD)
                htnet(ihem,jh,j,ld)=a*htnet(ihem,jh,1,ld)+b*htnet(ihem,jh,ilast,ld)
                TTRD(IM,LD)=(HTNET(IHEM,JH,J,LD)+HTNETO)/(CHRF*2.)

                IF (l.eq.nl) then
                  pnet(im,jh)=a*pnet(1+iofm,jh)+b*pnet(ilast+iofm,jh)
                  snet(im,jh)=a*snet(1+iofm,jh)+b*snet(ilast+iofm,jh)

                  DO k=1,6
                    rrflux(im,jh,k)=a*rrflux(1+iofm,jh,k)+b*rrflux(ilast+iofm,jh,k)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDIF

        ! This else goes with that very first if statement near the top lol
        ! This took forever to fix
        ELSE                   ! ntstep requirement
          DO i=1,mg
            DO LD=1,NL
              im=i+IOFM
              TTRD(im,LD)=(htnet(ihem,jh,i,ld))/CHRF
            ENDDO
          ENDDO
        ENDIF
        IOFM=MGPP
 800  CONTINUE

      IF (LSHORT.AND.(KOUNT.eq.1)) then
        DO l=1,nl
          DO i=1,igc
            ttrd(i,l)=ttrd(i,l)*2.
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END
