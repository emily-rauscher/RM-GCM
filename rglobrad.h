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
!
      PARAMETER (NLOW =2000)
      PARAMETER (NHIGH=2001)
      PARAMETER ( NCOUNT = NHIGH - NLOW )
!
! CONSTANT PARAMETERS THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL









!
! TIME-DEPENDENT VARIABLES THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL
!
!
! TIME-DEPENDENT VARIABLES THAT MIGHT BE SPECIFIED BY AN EXTERNAL MODEL
!
      COMMON /irradiation_parameters/
     &  EMISIR,
     &  HEATI(NLAYER),
     &  HEATS(NLAYER), HEAT(NLAYER),
     &  SOLNET,
     &  LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, EPSILON,
     &  TPI, SQ3, SBK,
     &  AM, AVG, ALOS, G, PI, SCDAY, RGAS,
     &  GANGLE(NGAUSS), GWEIGHT(NGAUSS),
     &  GRATIO(NGAUSS),
     &  EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL),
     &  SOL(NTOTAL),RAYPERBAR(NTOTAL),WEIGHT(NTOTAL),
     &  GOL(NTOTAL,NDBL),
     &  WOL(NTOTAL,NDBL),
     &  TAUCONST(NWAVE),
     &  iblackbody_above,
     &  WAVE(NWAVE+1),
     &  TT(NLAYER), Y3(NTOTAL,NGAUSS,NDBL),
     &  U0,  ISL, IR, IRS, FDEGDAY,
     &  WOT, GOT,
     &  PTEMPG(NTOTAL),
     &  PTEMPT(NTOTAL),
     &  G0(NTOTAL,NDBL),
     &  OPD( NTOTAL,NDBL),
     &  PTEMP(NTOTAL,NDBL),
     &  uG0(NTOTAL,NDBL),
     &  uTAUL(NTOTAL,NDBL),
     &  W0(NTOTAL,NDBL),
     &  uW0(NTOTAL,NDBL),
     &  uopd(NTOTAL,NDBL),
     &  U1S( NTOTAL),
     &  U1I( NTOTAL),
     &  TOON_AK(  NTOTAL,NDBL),
     &  B1(  NTOTAL,NDBL),
     &  B2(  NTOTAL,NDBL),
     &  EE1( NTOTAL,NDBL),
     &  EM1(NTOTAL,NDBL),
     &  EM2(NTOTAL,NDBL),
     &  EL1( NTOTAL,NDBL),
     &  EL2(NTOTAL,NDBL),
     &  GAMI(NTOTAL,NDBL),
     &  AF(NTOTAL,NQPL),
     &  BF(NTOTAL,NQPL),
     &  EF(NTOTAL,NQPL),
     &  SFCS(NTOTAL),
     &  B3(NTOTAL,NDBL),
     &  CK1(NTOTAL,NDBL),
     &  CK2(NTOTAL,NDBL),
     &  CP(NTOTAL,NDBL),
     &  CPB(NTOTAL,NDBL),
     &  CM(NTOTAL,NDBL),
     &  CMB(NTOTAL,NDBL),
     &  DIRECT(NTOTAL,NDBL),
     &  EE3(NTOTAL,NDBL),
     &  EL3(NTOTAL,NDBL),
     &  FNET(NTOTAL,NDBL),
     &  TMI(NTOTAL,NDBL),
     &  AS(NTOTAL,NQPL),
     &  DF(NTOTAL,NQPL),
     &  DS(NTOTAL,NQPL),
     &  XK(NTOTAL,NQPL),
     &  DIREC(  NTOTAL,NDBL),
     &  DIRECTU(NTOTAL,NDBL),
     &  DINTENT(NTOTAL,NGAUSS,NDBL),
     &  UINTENT(NTOTAL,NGAUSS,NDBL),
     &  TMID(NTOTAL,NDBL),
     &  TMIU(NTOTAL,NDBL),
     &  tslu,total_downwelling,alb_tot,tiru,firu(NIR),fird(NIR),fsLu(NSOL),
     &  fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL),
     &  fupbs(NLAYER),fdownbs(NLAYER),fnetbs(NLAYER),fdownbs2(nlayer),
     &  fupbi(NLAYER),fdownbi(NLAYER),fnetbi(NLAYER),
     &  qrad(NLAYER),alb_tomi,alb_toai

      common /irradiation_variables_kept_in_global/
     & zsin,O2(NLAYER), O3(NLAYER), AH2O(NTOTAL), TAUA(NWAVE,NDBL), dz(NLAYER)

! ensure all rad local variables are stored statically
      save