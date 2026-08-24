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
      INTEGER NIR,NSOL,NTOTAL
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
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS(3),with_TiO_and_VO, opacity_method

      LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF,RAYSCAT,AEROSOLS
      CHARACTER(len=6) :: opacity_method

      CHARACTER(30) :: AEROSOLMODEL

      REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),TCON(NL+1)
      LOGICAL DELTASCALE,HAZES,PICKET_FENCE_CLOUDS,GRAYCLDV

      COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &               CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &               ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &               MAXTAU,MAXTAULOC,TCON,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS,
     &               GRAYCLDV,C_TO_O

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

      real, dimension(NBATCH,2*NL+2) :: TAURAY, TAUL, TAUGAS, TAUAER

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
      ! integer TID, NTHREADS
      ! double precision test_wctime

      ! Thomas adding parallel stuff:
      INTEGER :: thread_num, istart, iend, nthreads
      REAL :: tstart, tend
      DATA IFIRST/1/
      data ifirstcol/1/

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NL+1), HEATS(NL+1), HEAT(NL+1), SOLNET
      REAL TPI, SQ3, SBK, AM, AVG, ALOS

      ! I just put a huge chunk of these in
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL), SOL(NTOTAL)
      REAL RAYPERBAR(NTOTAL),WEIGHT(NTOTAL)
      REAL GOL(NBATCH,2*NL+2), WOL(NBATCH,2*NL+2), WAVE(NTOTAL+1), TT(NL+1), Y3(NBATCH,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NBATCH), PTEMPT(NBATCH), G0(NBATCH,2*NL+2), OPD(NBATCH,2*NL+2), PTEMP(NBATCH,2*NL+2)
      REAL uG0(NBATCH,2*NL+2), uTAUL(NBATCH,2*NL+2), W0(NBATCH,2*NL+2), uW0(NBATCH,2*NL+2), uopd(NBATCH,2*NL+2),  U1S(NBATCH)
      REAL U1I(NBATCH), TOON_AK(NBATCH,2*NL+2), B1(NBATCH,2*NL+2), B2( NBATCH,2*NL+2), EE1(NBATCH,2*NL+2), EM1(NBATCH,2*NL+2)
      REAL EM2(NBATCH,2*NL+2), EL1( NBATCH,2*NL+2), EL2(NBATCH,2*NL+2), GAMI(NBATCH,2*NL+2), AF(NBATCH,4*NL+4)
      REAL BF(NBATCH,4*NL+4), EF(NBATCH,4*NL+4), SFCS(NBATCH), B3(NBATCH,2*NL+2), CK1(NBATCH,2*NL+2), CK2(NBATCH,2*NL+2)
      REAL CP(NBATCH,2*NL+2), CPB(NBATCH,2*NL+2), CM(NBATCH,2*NL+2), CMB(NBATCH,2*NL+2), DIRECT(NBATCH,2*NL+2), EE3(NBATCH,2*NL+2)
      REAL EL3(NBATCH,2*NL+2), FNET(NBATCH,2*NL+2), TMI(NBATCH,2*NL+2), AS(NBATCH,4*NL+4), DF(NBATCH,4*NL+4)
      REAL DS(NBATCH,4*NL+4), XK(NBATCH,4*NL+4), DIREC(NBATCH,2*NL+2), DIRECTU(NBATCH,2*NL+2), DINTENT(NBATCH,3,2*NL+2)
      REAL UINTENT(NBATCH,3,2*NL+2), TMID(NBATCH,2*NL+2), TMIU(NBATCH,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai, SLOPE(NBATCH,2*NL+2)
      real heats_aerad_tot(NL+1), heati_aerad_tot(NL+1), radheat_tot(NL+1), cheati(NL+1), cheats(NL+1), radheat(NL+1)

      REAL, DIMENSION(NBATCH,3,2*NL+2) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(NBATCH,2*NL+2)   :: A1, A2, A3, A4, A5, A7, Y5

      real, dimension(NKGAUSS, NL+1) :: k_IRl, tau_ray_temp
      real, dimension(NKGAUSS, NL+1) :: k_Vl

      ! For the new picket fence stuff
      REAL tau_IRe(NKGAUSS,NL+1), tau_Ve(NKGAUSS,NL+1)
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL)  :: Beta_V

      real, dimension(NL+1) :: dpe, Pl, Tl, pe
      real :: k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met
      real :: Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val


      REAL PI0_TEMP(NBATCH, NL+1, 13), G0_TEMP(NBATCH, NL+1, 13)
      REAL tauaer_temp(NBATCH, NL+1, 13)
      INTEGER j1
      real denom
      REAL, dimension (500) :: HAZE_WAV_GRID
      REAL, dimension (100)  :: CLOUD_WAV_GRID
      INTEGER printt
      REAL :: exp_92_lnsig2_pi

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

      COMMON /RAD_ALLLATS/ TG_forrad(IGC,NL,JG), PLG_forrad(IGC,JG),
     &                     HTNET_old(NHEM,JG,MG,NL)

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
      CHRF=86400.*WW*CT
      ntstep=NTSTEP_IN

      IOFM=0
      DO 800 ihem=1,nhem
        IF (mod(kount,ntstep) .eq. 0) THEN
          DO i=1,mg
            IM=I+IOFM
            DO l=nl,1,-1
              LD=NL+1-L
              TTRD(IM,LD)=(HTNET_old(ihem,IH,i,LD)+HTNET(ihem,IH,i,LD))
     &                    /(CHRF*2.0)
            ENDDO
          ENDDO
        ELSE
          DO i=1,mg
            DO LD=1,NL
              im=i+IOFM
              TTRD(im,LD)=(htnet(ihem,IH,i,ld))/CHRF
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


C**********************************************************
C             SUBROUTINE RADIATION_ALLLATS
C**********************************************************
      SUBROUTINE RADIATION_ALLLATS()
C     Parallelized over all JG*MG columns simultaneously.
C     Computes radiative heating rates and stores in HTNET.
C     HTNET_old is saved before the parallel loop.
C     TTRD is computed by the stripped-down RADIATION subroutine.

      use omp_lib
      INTEGER NIR,NSOL,NTOTAL
      include 'params.i'

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

      PARAMETER (N2DFLD=21,NGRPAD=N2DFLD*2*IGC)

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

      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC

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

      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)
     +              ,BEGDOY,DOY

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
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS(3),with_TiO_and_VO, opacity_method

      LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF,RAYSCAT,AEROSOLS
      CHARACTER(len=6) :: opacity_method

      CHARACTER(30) :: AEROSOLMODEL

      REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),TCON(NL+1)
      LOGICAL DELTASCALE,HAZES,PICKET_FENCE_CLOUDS,GRAYCLDV

      COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &               CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &               ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &               MAXTAU,MAXTAULOC,TCON,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS,
     &               GRAYCLDV,C_TO_O

      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC
     +              ,RNTAPO
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG
     +              ,LSPO(NL),LGPO(NL)
     $              ,LSHIST,LMINIH

      LOGICAL LSHIST,LMINIH
      LOGICAL LSPO,LGPO

      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1)),(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))

      COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)
     :     ,SNET(IGC,JG),RRFLUX(IGC,JG,6)
     :     ,TTSW(IGC,NL),TTLW(IGC,NL)

      COMMON/GSG/GSG(IGC,JG)

      REAL htnet

      COMMON /RADHT/ HTNET(NHEM,JG,MG,NL)

      COMMON /RAD_ALLLATS/ TG_forrad(IGC,NL,JG), PLG_forrad(IGC,JG),
     &                     HTNET_old(NHEM,JG,MG,NL)

      REAL TAVE(IGP)

      REAL PR(NL+1),T(NL+1),p_pass(nl+1),htlw(nl+1),htsw(nl+1)
      real dpg(nl+1), pbar(nl+1)
      real dpgsub(2*nl+2), pbarsub(2*nl+2)

      real, dimension(NBATCH,2*NL+2) :: TAURAY, TAUL, TAUGAS, TAUAER

      integer solar_calculation_indexer, num_layers

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

      integer ifirst
      real amfrac
      integer ichange
      integer ifirstcol
      real p0
      real ps
      integer im

      real fluxes(2,2,2)

      real incident_starlight_fraction

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NL+1), HEATS(NL+1), HEAT(NL+1), SOLNET
      REAL TPI, SQ3, SBK, AM, AVG, ALOS

      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL), SOL(NTOTAL)
      REAL RAYPERBAR(NTOTAL),WEIGHT(NTOTAL)
      REAL GOL(NBATCH,2*NL+2), WOL(NBATCH,2*NL+2), WAVE(NTOTAL+1), TT(NL+1), Y3(NBATCH,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NBATCH), PTEMPT(NBATCH), G0(NBATCH,2*NL+2), OPD(NBATCH,2*NL+2), PTEMP(NBATCH,2*NL+2)
      REAL uG0(NBATCH,2*NL+2), uTAUL(NBATCH,2*NL+2), W0(NBATCH,2*NL+2), uW0(NBATCH,2*NL+2), uopd(NBATCH,2*NL+2),  U1S(NBATCH)
      REAL U1I(NBATCH), TOON_AK(NBATCH,2*NL+2), B1(NBATCH,2*NL+2), B2( NBATCH,2*NL+2), EE1(NBATCH,2*NL+2), EM1(NBATCH,2*NL+2)
      REAL EM2(NBATCH,2*NL+2), EL1( NBATCH,2*NL+2), EL2(NBATCH,2*NL+2), GAMI(NBATCH,2*NL+2), AF(NBATCH,4*NL+4)
      REAL BF(NBATCH,4*NL+4), EF(NBATCH,4*NL+4), SFCS(NBATCH), B3(NBATCH,2*NL+2), CK1(NBATCH,2*NL+2), CK2(NBATCH,2*NL+2)
      REAL CP(NBATCH,2*NL+2), CPB(NBATCH,2*NL+2), CM(NBATCH,2*NL+2), CMB(NBATCH,2*NL+2), DIRECT(NBATCH,2*NL+2), EE3(NBATCH,2*NL+2)
      REAL EL3(NBATCH,2*NL+2), FNET(NBATCH,2*NL+2), TMI(NBATCH,2*NL+2), AS(NBATCH,4*NL+4), DF(NBATCH,4*NL+4)
      REAL DS(NBATCH,4*NL+4), XK(NBATCH,4*NL+4), DIREC(NBATCH,2*NL+2), DIRECTU(NBATCH,2*NL+2), DINTENT(NBATCH,3,2*NL+2)
      REAL UINTENT(NBATCH,3,2*NL+2), TMID(NBATCH,2*NL+2), TMIU(NBATCH,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai, SLOPE(NBATCH,2*NL+2)
      real heats_aerad_tot(NL+1), heati_aerad_tot(NL+1), radheat_tot(NL+1), cheati(NL+1), cheats(NL+1), radheat(NL+1)

      REAL, DIMENSION(NBATCH,3,2*NL+2) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(NBATCH,2*NL+2)   :: A1, A2, A3, A4, A5, A7, Y5

      real, dimension(NKGAUSS, NL+1) :: k_IRl, tau_ray_temp
      real, dimension(NKGAUSS, NL+1) :: k_Vl

      REAL tau_IRe(NKGAUSS,NL+1), tau_Ve(NKGAUSS,NL+1)
      real, dimension(NIR)  :: Beta_IR
      real, dimension(NSOL)  :: Beta_V

      real, dimension(NL+1) :: dpe, Pl, Tl, pe
      real :: k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met
      real :: Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val

      REAL PI0_TEMP(NBATCH, NL+1, 13), G0_TEMP(NBATCH, NL+1, 13)
      REAL tauaer_temp(NBATCH, NL+1, 13)

C     Large workspace arrays: put in COMMON /RADWORK/ so they are in
C     static (non-stack) storage; THREADPRIVATE gives each OMP thread
C     its own copy without the stack-overflow hazard of PARAMETER-sized
C     local arrays under -recursive.
      COMMON /RADWORK/ GOL, WOL, G0, OPD, PTEMP, uG0, uTAUL, W0,
     &  uW0, uopd, TOON_AK, B1, B2, EE1, EM1, EM2, EL1, EL2, GAMI,
     &  AF, BF, EF, B3, CK1, CK2, CP, CPB, CM, CMB, DIRECT, EE3,
     &  EL3, FNET, TMI, AS, DF, DS, XK, DIREC, DIRECTU, TMID, TMIU,
     &  DINTENT, UINTENT, SLOPE, Y3, Y1, Y2, Y4, Y8, A1, A2, A3,
     &  A4, A5, A7, Y5, TAURAY, TAUL, TAUGAS, TAUAER,
     &  PI0_TEMP, G0_TEMP, tauaer_temp
!$OMP THREADPRIVATE(/RADWORK/)

      INTEGER j1
      real denom
      REAL, dimension (500) :: HAZE_WAV_GRID
      REAL, dimension (100)  :: CLOUD_WAV_GRID
      INTEGER printt
      REAL :: exp_92_lnsig2_pi

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

C     Local scalars
      INTEGER :: col, i, jh_priv, ihem, ihem_save, nskip, ntstep
      INTEGER :: iofm, ld, l, thread_num, nthreads
      REAL :: CHRF, alat1, alon, SWALB, tstart, tend
      REAL :: SSLON, SSLAT

      num_layers = NL
      CHRF=86400.*WW*CT
      ntstep=NTSTEP_IN
      nskip=NSKIP_IN

C     Save HTNET_old before any updates (serial, before OMP region)
      DO ihem_save=1,NHEM
        DO jh_priv=1,JG
          DO i=1,MG
            DO ld=1,NL
              HTNET_old(ihem_save,jh_priv,i,ld) =
     &          HTNET(ihem_save,jh_priv,i,ld)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO 900 ihem=1,nhem
        iofm=0
        IF (ihem.eq.2) iofm=MGPP

C       LFLUXDIAG block: must be serial, placed before OMP region
        IF ((LFLUXDIAG).AND.(KOUNTP-KOUTP.LT.NTSTEP_IN)) THEN
          IF (ihem.eq.1) THEN
            REWIND(63)
            REWIND(62)
            IF (PORB.NE.0) THEN
              SSLON=(1./PORB-1.)*KOUNT*360./ITSPD
              SSLON=MOD(SSLON,360.)
            ELSE
              SSLON=0.
            ENDIF
            SSLAT=ASIN(SIN(OBLIQ*PI/180.)*SIN(PI2*KOUNT/ITSPD/PORB))
     &            *180./PI
            WRITE(63,2021) DAY,SSLON,SSLAT
            WRITE(62,2021) DAY,SSLON,SSLAT
 2021       FORMAT('DAY:',F7.2,', SUBSTELLAR LON,LAT:',2F7.2)
            WRITE(63,*)
            WRITE(62,*)''
          ENDIF
        ENDIF

        IF (nskip.ne.0) then
          write(*,*) 'CANNOT SKIP LONGITUDES IN PARALLEL!! ABORT'
          write(*,*) 'Please set nskip=0 in fort.7'
          STOP
        ENDIF

!$OMP   PARALLEL default(none)
!$OMP&  private(col, i, jh_priv, ld, l, PR, PRB2T, T, AEROPROF,
!$OMP&  p_pass, alon, rfluxes_aerad, fluxes, fsl_dn_aerad, fir_up_aerad,
!$OMP&  k_irl, k_vl, htlw, htsw, thread_num,
!$OMP&  psol_aerad, cheati, pl, dpg, tin, k_lowp,
!$OMP&  u0, sfcs, pbar, cheats,
!$OMP&  lla, heat, heati, heati_aerad,
!$OMP&  fnetbi, alb_toai, emis,
!$OMP&  fird, radheat_tot, alb_toa, freedman_met,
!$OMP&  beta_v, denom, weight, fupbs, tslu, fsl_net_aerad, irs,
!$OMP&  fsld, tiru, jdble, dpgsub, fslu, fdegday, dpe, u1s,
!$OMP&  solnet, freedman_p, jn2, sol, u1i, swalb, am,
!$OMP&  jn, pin, fupbi, qrad, isl, k_hip,
!$OMP&  alat1, beta_ir, fir_dn_aerad, nprob,
!$OMP&  fdownbi, total_downwelling, tau_ire, j1,
!$OMP&  ptempg, iblackbody_above, emisir,
!$OMP&  heats_aerad_tot, fnetbs, sq3, pressure_val, tl, ifsetup,
!$OMP&  heati_aerad_tot, lls, solar_calculation_indexer,
!$OMP&  alb_tomi, fsl_up_aerad, firu, fir_net_aerad, radheat, cf,
!$OMP&  got, fsln, temperature_val, rsfx,
!$OMP&  incident_starlight_fraction, fdownbs, wot,
!$OMP&  tt, tau_ve, k_ir, tpi, wave,
!$OMP&  ibinm, heats_aerad, jdbledble, fdownbs2, tl10, pbarsub,
!$OMP&  heats, alb_tot, pe, pl10, freedman_t,
!$OMP&  ptempt, ic, ir, im,
!$OMP&  gauss_idx, wave_idx, stel_idx, chan_idx, J, K, T_idx, P_idx,
!$OMP&  temp_idx, index_num, lo_temp_flag, it1, kindex,
!$OMP&  iffirst, tgrnd, ibinmin, log_start, log_end, log_step,
!$OMP&  P_pass_sub, ir_abs_coefficient, wavea, ttsub, albedoa,
!$OMP&  tau_ray_temp)
!$OMP&  shared(ntstep, nskip, lnnsk, sigma, GSG, P0, CT, FBASEFLUX,
!$OMP&  rrflux, alat, lfluxdiag, kountp, koutp,
!$OMP&  ntstep_in, porb, kount, itspd, obliq, day, albsw, aerosols,
!$OMP&  aerosolmodel, tauaerosol, doy, epsilon,
!$OMP&  avg, alos, SCDAY, RGAS, GANGLE, GWEIGHT, GRATIO, RAYPERBAR,
!$OMP&  sbk, num_layers, CHRF,
!$OMP&  PNET, SNET, HTNET, ihem, iofm,
!$OMP&  TG_forrad, PLG_forrad, T0)

!$OMP   DO SCHEDULE(GUIDED)
        DO col=1,JG*MG
          jh_priv = (col-1)/MG + 1
          i       = MOD(col-1,MG) + 1
          im = i + iofm

          DO LD=1,NL
            L=NL-LD+2
            PR(LD)=SIGMA(LD)*EXP(PLG_forrad(im,jh_priv))*P0
            PRB2T(L)=PR(LD)
            T(LD)=(TG_forrad(im,ld,jh_priv)+T0(LD))*CT
            AEROPROF(LD)=0.0
          ENDDO

          AEROPROF(NL+1)=0.0
          PRB2T(1)=EXP(PLG_forrad(im,jh_priv))*P0
          PR(NL+1)=EXP(PLG_forrad(im,jh_priv))*P0
          T(NL+1)=((FBASEFLUX+rrflux(IM,jh_priv,1))/5.6704e-8)**0.25

          alat1=alat(jh_priv)*REAL(-(ihem*2.)+3)

          SWALB=ALBSW

          DO LD=1,NL-1
            p_pass(LD+1)=(pr(LD)+pr(LD+1))/2.
          ENDDO
          p_pass(NL+1)=PR(NL+1)
          p_pass(1)=pr(1)*0.5

          alon=REAL(i-1)/REAL(mg)*360.0

          IF((AEROSOLS).AND.(AEROSOLMODEL.NE.'Global')) THEN
            DO LD=1,NL+1
              AEROPROF(LD)=TAUAEROSOL(LD,i,ihem,jh_priv)
            ENDDO
          ENDIF

          rfluxes_aerad = 0.
          fluxes        = 0.
          fsl_dn_aerad  = 0.
          fir_up_aerad  = 0.
          Y3  = 0.
          Y1  = 0.
          Y5  = 0.
          Y4  = 0.
          k_IRl = 0
          k_Vl  = 0

          call calc_radheat(pr,t,p_pass,alat1,alon,htlw,htsw,
     &                      DOY,cf,ic,fluxes,swalb,kount,itspd,
     &                      incident_starlight_fraction,TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, dpg,
     &       ifsetup, ibinm, rfluxes_aerad, psol_aerad, heati_aerad, heats_aerad,
     &       fsl_up_aerad, fsl_dn_aerad, fir_up_aerad, fir_dn_aerad, fir_net_aerad, fsl_net_aerad,
     &       pbar, dpgsub, pbarsub,
     &       LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,
     &       EMISIR, EPSILON, HEATI, HEATS, HEAT, SOLNET, TPI, SQ3, SBK, AM, AVG, ALOS,
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
     &  PI0_TEMP, G0_TEMP, tauaer_temp, j1, denom, Beta_IR, Beta_V, k_IRl, k_Vl, tau_ray_temp)

          pr=prb2t

          PNET(IM,jh_priv)=fluxes(1,1,1)-fluxes(1,2,1)
     &                    +fluxes(2,1,1)-fluxes(2,2,1)
          SNET(IM,jh_priv)=fluxes(1,1,2)-fluxes(1,2,2)
     &                    +fluxes(2,1,2)-fluxes(2,2,2)

          rrflux(im,jh_priv,1)=fluxes(1,1,2)
          rrflux(im,jh_priv,2)=fluxes(1,2,2)
          rrflux(im,jh_priv,3)=fluxes(2,1,2)
          rrflux(im,jh_priv,4)=fluxes(2,2,2)
          rrflux(im,jh_priv,5)=fluxes(1,1,1)-fluxes(1,2,1)
          rrflux(im,jh_priv,6)=fluxes(2,2,1)

          DO l=nl,1,-1
            LD=NL+1-L
            htnet(ihem,jh_priv,i,ld)=(htlw(l+1)+htsw(l+1))
          ENDDO

        ENDDO
!$OMP   END DO
!$OMP   END PARALLEL

 900  CONTINUE

      RETURN
      END
