c
c  Declare and define symbolic constants
c  This simply cleans up the code and prevents me from
c  introducing a bunch of common blocks

       include 'params.i'
       include 'rglobrad.h'
C                                                                       
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN,TAULIMIT

       LOGICAL LPLOTMAP

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM, AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE, RAYPERBARCONS
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS

    
       CHARACTER(30) :: AEROSOLMODEL
       CHARACTER(30) :: AEROSOLCOMP
       REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),MAXTAU,TCON(NL+1)
       REAL MTLX
       INTEGER MAXTAULOC,AERLAYERS
       LOGICAL DELTASCALE
       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &   MAXTAU,MAXTAULOC,TCON,AEROSOLCOMP,MTLX,MOLEF,AERLAYERS
      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC                        
     +              ,RNTAPO                                               
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR                          
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY                          
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG                 
     +              ,LSPO(NL),LGPO(NL)                                    
     $              ,LSHIST,LMINIH                                        
      LOGICAL LSHIST,LMINIH                                               
      LOGICAL LSPO,LGPO

! need this for the logical switch for if binary or not.
      COMMON/BINVAL/PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2,STRAD1,STRAD2,STTEMP1,STTEMP2,LBIN

      LOGICAL LBIN


      INTEGER MXGAS,MXLEV,MXBAND,MXCL                                     
                                                                          
      PARAMETER(MXGAS=8)     ! Maximum number of gases                    
                                                                          
      PARAMETER(MXLEV=NL)     ! Maximum number of levels in the            
                             ! atmosphere (not including the surface)                                                                            
      PARAMETER(MXBAND=9)    ! Maximum number of spectral bands (not      
                             ! including the whole spectrum, 0-3000cm-1)  
                                                                          
      PARAMETER(MXCL=3)      ! Maximum number of cloud types              
                                                                          
C-----------------------------------------------------------------------  
C                                                                         
C  Switches for long wave radiation scheme                                
C                                                                         
      COMMON/RADLW/VMRCO2,VMRCH4,VMRN2O,VMRHALO,GAS(MXGAS),LLBLM          
      LOGICAL LLBLM                                                       
      INTEGER GAS                                                         
      REAL VMRCO2,VMRCH4,VMRN2O,VMRHALO                                   

C     Sets basic constants, especially those needed for array dimensions  
C                                                                         
      PARAMETER(MH=2                        
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
C                                                                         
C     Note that RD and GASCON are identical and CPD is set from RD,AKAP.  
      COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)                
     +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)     
     +              ,TKP(NL),C(NL2),SQH(NNP)
     +              ,MF,MFP,JZF,NF
     +              ,AKAP,GA,GASCON,RADEA,WW,PFAC,EZ,AIOCT
     +              ,RD,RV,CPD,CLATNT                                     
     +              ,P0,LRSTRT,LSHORT,LTVEC,LSTRETCH                         
     +              ,LFLUX,R_AIR                                                
     +              ,LBALAN,LRESTIJ                                       
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
     +              ,LNOISE,NFP                                               
      COMPLEX EZ,AIOCT                                                    
      LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ                 
     +       ,LFLUX,LNOISE                                                
     +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                           
C    
 
!      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
!     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
!     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
!     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)
!     +              ,BEGDOY,DOY 
      PARAMETER (NRLEV=MXLEV+1)                                            
********************************************************************      
       save
