!
!  Declare and define symbolic constants
!   (Implicit typing in effect unless explicitly specified)
!
!
!  Include implicit declarations
!
      include 'rprecision.h'
!
!
!  Include symbolic constants shared between aerosol and radiation models
!
      include 'raerad.h'
!
!
! DEFINE THE DIMENSIONS USED BY THE RADIATION MODEL
!
! NVERT  = MAXIMUM NUMBER OF LAYERS;
! NLAYER = MAXIMUM NUMBER OF LAYER BOUNDARIES
! NDBL   = TWICE THE MAXIMUM NUMBER OF LAYER BOUNDARIES
! NRAD   = MAXIMUM NUMBER OF AEROSOL RADIUS BINS;
!
      PARAMETER ( NVERT = NL )
      PARAMETER ( NRAD = NBIN )
      PARAMETER ( NLAYER = NVERT+1)
      PARAMETER ( NDBL = 2*NLAYER )
      PARAMETER ( NRADVER = NRAD*NVERT )
      PARAMETER ( NRADLAY = NRAD*NLAYER )
      PARAMETER ( NQPL = 4*NLAYER)
!
! NTOTAL = TOTAL NUMBER OF PROBABILITY INTERVALS;
! NSOLP  = NUMBER OF SOLAR PROBABILITY INTERVALS;
! NIRP   = NUMBER OF INFRARED PROBABILITY INTERVALS;
!
      PARAMETER ( NSOLP = 3)
      PARAMETER ( NIRP  = 2)
      PARAMETER ( NTOTAL = NSOLP + NIRP )
!
! NGAUSS = TOTAL NUMBER OF GAUSS QUADRATURE POINTS;
!
      PARAMETER ( NGAUSS = 3)
!
! NCOUNT = USED TO CALCULATE PLANK FUNCTION.
!
      PARAMETER (NLOW =2000)
      PARAMETER (NHIGH=2001)
      PARAMETER ( NCOUNT = NHIGH - NLOW )
!
! CONSTANT PARAMETERS THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL
!
      character*(50) prtofil
      character*(50) radofil
      logical is_grp_ice

      COMMON /irad1/
     &  O3MIX(NLAYER), O3MIXP(6), O3C, VRAT,
     &  PTOP, PBOT, RMIN(NGROUP), R(NRAD,NGROUP),
     &  is_grp_ice(NGROUP)
!
! TIME-DEPENDENT VARIABLES THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL
!
      COMMON /irad2/
     &  U0EXT,
     &  ALBEDO_SFC, EMISIR
!
! OUTPUT VARIABLES, CALCULATED BY THE RADIATION MODEL, THAT MIGHT
! BE USED BY THE EXTERNAL MODEL
!
      COMMON /irad3/
     &  HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER),
     &  SOLNET, XIRDOWN, XIRUP
!
! INITIATED IN SETUPRAD FOR RADIATION CALCULATION
!
      COMMON /irad4/
     &   LLA, LLS, JDBLE,JDBLEDBLE,JN,JN2, EPSILON, EXPMAX,
     &   TPI, SQ3, SBK,
     &   AM, AVG, ALOS, G, PI, SCDAY, RGAS,
     &   GANGLE(NGAUSS), GWEIGHT(NGAUSS), SFLX(NSOL), WVLN(NSOL),
     &   GRATIO(NGAUSS),
     &   EMIS(NTOTAL), RSFX(NTOTAL),LTEMP(NTOTAL),NPROB(NTOTAL),
     &   SOL(NTOTAL),RAYPERBAR(NTOTAL),WEIGHT(NTOTAL),
     &   GCLD(  NTOTAL,NDBL),   GOL(NTOTAL,NDBL),
     &   TAURAY( NTOTAL,NDBL),TAUAER(NTOTAL,NDBL),
     &   WCLD(NTOTAL,NDBL),
     &   TAUCLD(NTOTAL,NDBL),WOL(NTOTAL,NDBL),
     &   TREAL(2,NWAVE), TTMAG(2,NWAVE),
     &   contnm(NIRP), nprobi(NWAVE,2),TAUCONST(NWAVE)
!
      COMMON/irad5/
     &        iblackbody_above,       t_above,
     &        ACO2(NTOTAL),           AH2O(NTOTAL),
     &        AO2(NTOTAL),            AO3(NTOTAL),
     &        PACO2(NTOTAL,NLAYER),   PAH2O(NTOTAL,NLAYER),
     &        PAO2(NTOTAL,NLAYER),    PAO3(NTOTAL,NLAYER),
     &        PLANK(NIR+1,NCOUNT),
     &        PSCO2(NTOTAL),          PSH2O(NTOTAL),
     &        PSO2(NTOTAL),           PSO3(NTOTAL),
     &        SOLFX(NSOL),            WAVE(NWAVE+1),
     &        TAUGAS(NTOTAL,NDBL),
     &        XSECTA(NRAD,NGROUP),    RUP(NRAD,NGROUP),
     &       QSCAT(NRAD,NGROUP,NWAVE),
     &       QBRQS(NRAD,NGROUP,NWAVE),
     &        RDQEXT(NRAD,NGROUP,NWAVE)
!
      COMMON/irad6/   CO2(NLAYER), RDH2O(NLAYER),   O2(NLAYER),
     &                 O3(NLAYER), CAER(NRAD,NLAYER,NGROUP),
     &                 PBAR(NLAYER),PLAYER(NLAYER),
     &                 DPG(NLAYER), TT(NLAYER), Y3(NTOTAL,NGAUSS,NDBL),
     &                 PRESSMID(NLAYER), DPGsub(NDBL) , PBARsub(NDBL),
     &                 TTsub(NDBL),
     &                 TGRND,  U0,  ISL, IR, IRS, FDEGDAY
!
! DEFINED IN 'OPPROP'
!
      COMMON /irad7/
     &   WOT, GOT,
     &   PTEMPG(NTOTAL),      PTEMPT(NTOTAL),
     &   G0(NTOTAL,NDBL),     OPD( NTOTAL,NDBL),
     &   PTEMP(NTOTAL,NDBL),  TAUL(NTOTAL,NDBL),
     &   TAUH2O(NTOTAL,NDBL), TAUS(NWAVE,NDBL),
     &   TAUA(NWAVE,NDBL),    G01(NWAVE,NDBL),
     &   uG0(NTOTAL,NDBL),    uTAUL(NTOTAL,NDBL),
     &   W0(NTOTAL,NDBL),     uW0(NTOTAL,NDBL),
     &   uopd(NTOTAL,NDBL)

!
! DEFINED IN 'TWOSTR'
!
      COMMON /irad8/
     &  U1S( NTOTAL),           U1I( NTOTAL),
     &   ACON(NTOTAL,NDBL),   TOON_AK(  NTOTAL,NDBL),
     &   BCON(NTOTAL,NDBL),   B1(  NTOTAL,NDBL),
     &   B2(  NTOTAL,NDBL),   EE1( NTOTAL,NDBL),
     &   EM1(NTOTAL,NDBL),
     &   EM2(NTOTAL,NDBL),    EL1( NTOTAL,NDBL),
     &   EL2(NTOTAL,NDBL),    GAMI(NTOTAL,NDBL),
     &   AF(NTOTAL,NQPL), BF(NTOTAL,NQPL), EF(NTOTAL,NQPL)
!
! DEFINED IN 'ADD'
!
      COMMON /irad9/
     &   SFCS(NTOTAL),
     &   B3(  NTOTAL,NDBL),   CK1(   NTOTAL,NDBL),
     &   CK2( NTOTAL,NDBL),   CP(    NTOTAL,NDBL),
     &   CPB( NTOTAL,NDBL),   CM(    NTOTAL,NDBL),
     &   CMB( NTOTAL,NDBL),   DIRECT(NTOTAL,NDBL),
     &   EE3( NTOTAL,NDBL),   EL3(   NTOTAL,NDBL),
     &   FNET(NTOTAL,NDBL),   TMI(   NTOTAL,NDBL),
     &   AS(  NTOTAL,NQPL),     DF(    NTOTAL,NQPL),
     &   DS(  NTOTAL,NQPL),     XK(    NTOTAL,NQPL)
!
! DEFINED IN 'NEWFLUX1'
!
      COMMON /irad10/
     &   WEIT(   NTOTAL),
     &   DIREC(  NTOTAL,NDBL), DIRECTU(NTOTAL,NDBL),
     &   SLOPE(  NTOTAL,NDBL),
     &   DINTENT(NTOTAL,NGAUSS,NDBL),
     &   UINTENT(NTOTAL,NGAUSS,NDBL),
     &   TMID(NTOTAL,NDBL),TMIU(NTOTAL,NDBL)
!
! defined in 'radtran'
!
      common /irad11/
     &   tslu,tsld,alb_tot,tiru,firu(NIR),fird(NIR),fsLu(NSOL),
     &   fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL),
     &   fupbs(NLAYER),fdownbs(NLAYER),fnetbs(NLAYER),fdownbs2(nlayer),
     &   fupbi(NLAYER),fdownbi(NLAYER),fnetbi(NLAYER),
     &   qrad(NLAYER),alb_tomi,alb_toai
!
! defined for dharma
!
!  Define values of flag used to specify calculation of solar zenith angle
!
      parameter( I_FIXED = 0 )
      parameter( I_DIURNAL = 1 )

      logical is_bin_micro, do_mie, use_qv, use_ice

      common /irad12/
     &  is_bin_micro, do_mie, use_qv, use_ice,
     &  myproc, ifix_u0, ifix_sfc_alb, ienconc(NGROUP),
     &  sfc_alb, sfc_wind, zsin, zcos, t_offset,
     &  theta_o, rLscp, rLiscp,
     &  drop_n, drop_sig, haze_n, haze_rg, haze_sig,
     &  ship_n, ig_soot, r_soot,
     &  pienv(NVERT), rhobar(NVERT), dz(NLAYER)
!
! ensure all rad local variables are stored statically
!
      save