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
!                                                                       
! NTOTAL = TOTAL NUMBER OF PROBABILITY INTERVALS;                       
! NSOLP  = NUMBER OF SOLAR PROBABILITY INTERVALS;                       
! NIRP   = NUMBER OF INFRARED PROBABILITY INTERVALS;                    
!                                                                       
!      PARAMETER ( NSOLP = 77 )
!      PARAMETER ( NIRP = 71 )
      PARAMETER ( NSOLP = 1)
      PARAMETER ( NIRP = 1)
      PARAMETER ( NTOTAL = NSOLP + NIRP )          
!                                                                       
! NGAUSS = TOTAL NUMBER OF GAUSS QUADRATURE POINTS;                     
!                                                                       
      PARAMETER ( NGAUSS = 3)                                           
!                                                                       
! NCOUNT = USED TO CALCULATE PLANK FUNCTION.                            
!                                                                       
!      PARAMETER ( NLOW = 12500 )
!      PARAMETER ( NHIGH = 32500 )
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
     &  ALBEDO_SFC, EMISIR !,     
!     &  P(NVERT), T(NVERT), Q(NVERT)
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
     &   LLA, LLS, JDBLE, JN, EPSILON, EXPMAX,                               
     &   TPI, SQ3, SBK,                                                       
     &   AM, AVG, ALOS, G, PI, SCDAY, RGAS,                                   
     &   GANGLE(NGAUSS), GWEIGHT(NGAUSS), SFLX(NSOL), WVLN(NSOL),             
     &   GRATIO(NGAUSS),     
     &   EMIS(NTOTAL), RSFX(NTOTAL),LTEMP(NTOTAL),NPROB(NTOTAL),              
     &   SOL(NTOTAL),RAYPERBAR(NTOTAL),WEIGHT(NTOTAL),                           
     &   GCLD(  NTOTAL,NLAYER),   GOL(NTOTAL,NLAYER),                         
     &   TAURAY( NTOTAL,NLAYER),TAUAER(NTOTAL,NLAYER),                         
     &   WCLD(NTOTAL,NLAYER),                                                 
     &   TAUCLD(NTOTAL,NLAYER),WOL(   NTOTAL,NLAYER),                         
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
     &        TAUGAS(NTOTAL,NLAYER),       
     &        XSECTA(NRAD,NGROUP),    RUP(NRAD,NGROUP),     
     &       QSCAT(NRAD,NGROUP,NWAVE),      
     &       QBRQS(NRAD,NGROUP,NWAVE),      
     &        RDQEXT(NRAD,NGROUP,NWAVE)                                         
!                                                                       
      COMMON/irad6/   CO2(NLAYER), RDH2O(NLAYER),   O2(NLAYER),              
     &                 O3(NLAYER), CAER(NRAD,NLAYER,NGROUP),      
     &                 PBAR(NLAYER),PLAYER(NLAYER),                                          
     &                 DPG(NLAYER), TT(NLAYER), Y3(NTOTAL,NGAUSS,NLAYER),     
     &                 TGRND,  U0,  ISL, IR, IRS, FDEGDAY,FLXLIMDIF       
!                                                                       
! DEFINED IN 'OPPROP'                                                   
!
!MTR      COMMON /irad7/                                                        
!MTR     &   WOT, GOT,                                                            
!MTR     &   PTEMPG(NTOTAL),        PTEMPT(NTOTAL),     
!MTR     &   G0(   NTOTAL,NLAYER),  OPD( NTOTAL,NLAYER),                          
!MTR     &   PTEMP(NTOTAL,NLAYER),  TAUL(NTOTAL,NLAYER),     
!MTR     &   TAUH2O(NTOTAL,NLAYER), TAUS(NWAVE,NLAYER),     
!MTR     &   TAUA(NWAVE,NLAYER),    G01(NWAVE,NLAYER),     
!MTR     &   uG0(   NTOTAL,NLAYER), uTAUL(NTOTAL,NLAYER),     
!MTR     &   W0(   NTOTAL,NLAYER),                                                
!MTR     &   uW0(   NTOTAL,NLAYER) , uopd(NTOTAL,NLAYER)
      COMMON /irad7/
     &   WOT, GOT,
     &   PTEMPG(NTOTAL),        PTEMPT(NTOTAL),  
! SHOULD THESE BE NIR INSTEAD OF NTOTAL
     &   G0(   NTOTAL,NLAYER),  OPD( NTOTAL,NLAYER),
     &   PTEMP(NTOTAL,NLAYER),  TAUL(NTOTAL,NLAYER),
     &   TAUH2O(NTOTAL,NLAYER), TAUS(NWAVE,NLAYER),
     &   TAUA(NWAVE,NLAYER),    G01(NWAVE,NLAYER),
     &   uG0(   NTOTAL,NLAYER), uTAUL(NTOTAL,NLAYER),
     &   W0(   NTOTAL,NLAYER),
     &   uW0(   NTOTAL,NLAYER) , uopd(NTOTAL,NLAYER)

!
! DEFINED IN 'TWOSTR'                                                   
!
      COMMON /irad8/                                                       
     &  U1S( NTOTAL),           U1I( NTOTAL),                                
     &   ACON(NTOTAL,NLAYER),   AK(  NTOTAL,NLAYER),                          
     &   BCON(NTOTAL,NLAYER),   B1(  NTOTAL,NLAYER),                          
     &   B2(  NTOTAL,NLAYER),   EE1( NTOTAL,NLAYER),                          
     &   EM1(NTOTAL,NLAYER),                                                  
     &   EM2(NTOTAL,NLAYER),    EL1( NTOTAL,NLAYER),                          
     &   EL2(NTOTAL,NLAYER),    GAMI(NTOTAL,NLAYER),                          
     &   AF(NTOTAL,NDBL), BF(NTOTAL,NDBL), EF(NTOTAL,NDBL)               
!
! DEFINED IN 'ADD'                                                      
!
      COMMON /irad9/                                                        
     &   SFCS(NTOTAL),                                                        
     &   B3(  NTOTAL,NLAYER),   CK1(   NTOTAL,NLAYER),                        
     &   CK2( NTOTAL,NLAYER),   CP(    NTOTAL,NLAYER),                        
     &   CPB( NTOTAL,NLAYER),   CM(    NTOTAL,NLAYER),                        
     &   CMB( NTOTAL,NLAYER),   DIRECT(NTOTAL,NLAYER),                       
     &   EE3( NTOTAL,NLAYER),   EL3(   NTOTAL,NLAYER),                        
     &   FNET(NTOTAL,NLAYER),   TMI(   NTOTAL,NLAYER),                        
     &   AS(  NTOTAL,NDBL),     DF(    NTOTAL,NDBL),                          
     &   DS(  NTOTAL,NDBL),     XK(    NTOTAL,NDBL)                      
!
! DEFINED IN 'NEWFLUX1'                                                 
!
      COMMON /irad10/                                                        
     &   WEIT(   NTOTAL),                                                    
     &   DIREC(  NTOTAL,NLAYER), DIRECTU(NTOTAL,NLAYER),                      
     &   SLOPE(  NTOTAL,NLAYER), 
     &   DINTENT(NTOTAL,NGAUSS,NLAYER),                                       
     &   UINTENT(NTOTAL,NGAUSS,NLAYER),                                       
     &   TMID(NTOTAL,NLAYER),TMIU(NTOTAL,NLAYER)
!
! defined in 'radtran'
!
      common /irad11/     
     &   tslu,tsld,alb_tot,tiru,firu(NIR),fird(NIR),fsLu(NSOL),     
     &   fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL),     
     &   fupbs(NLAYER),fdownbs(NLAYER),fnetbs(NLAYER),     
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

