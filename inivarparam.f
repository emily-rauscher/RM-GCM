C**********************************************************               
C             SUBROUTINE INIVARPARAM                                          
C**********************************************************               
      subroutine INIVARPARAM                                                 
C-----------------------------------------------------------------------  
C     Subroutine to initialise various model parameters            
C-----------------------------------------------------------------------  
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'                                                                                    
C                                                                         
C                                                                         
C     Sets basic constants, especially those needed for array dimensions  
C                                                                         
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
C
C      
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
                                                           
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        

      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C

       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN, TAULIMIT
      
       LOGICAL LPLOTMAP

C      Settling-dye parameters, shared with DGRMLT and ICTRAC. Their
C      own common, so the ~15 routines declaring VARPARAM need no
C      change. Read here rather than in ICTRAC because INIVARPARAM runs
C      unconditionally, while ICTRAC only runs when .NOT.LRSTRT.
C      ADYE is indexed by tracer, so ADYE(1) is the unused water slot.
       COMMON/DYEPAR/ADYE(NTRAC),RHODYE,PDYEFIX,PDYEUPPER,TRELAXORB

C      Lagged dye field for DGRMLT's implicit settling sweep. Seeded
C      to a negative sentinel here because INIVARPARAM is the one
C      routine that runs on both cold starts and restarts; DGRMLT
C      swaps in the real field on first touch. Mixing ratios are
C      non-negative, so a negative value is unambiguous.
       COMMON/DYEPRV/TRAPRV(IGC,NL,JG,NTRAC-1)

       COMMON/MAG/ BFIELD,TDRAG_MIN,RAMPUP,LBDRAG
       LOGICAL LBDRAG

       COMMON/BINVAL/PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2,STRAD1,STRAD2,STTEMP1,STTEMP2,LBIN

       LOGICAL LBIN

       NAMELIST/INVARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN,
     & ADYE, RHODYE, PDYEFIX, PDYEUPPER, TRELAXORB

       NAMELIST/INMAG/ LBDRAG,BFIELD,TDRAG_MIN,RAMPUP

       NAMELIST/INBINVAL/LBIN,PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2, STRAD1,STRAD2,STTEMP1,STTEMP2

       LPLOTMAP=.TRUE.

C      Defaults before the namelist read, so a fort.7 predating these
C      keys still gives a defined configuration.
       DO 5 KK=1,NTRAC
          ADYE(KK)=1.0E-6
    5  CONTINUE
       RHODYE=3.0E3
       PDYEFIX=1.0E5
       PDYEUPPER=1.0E6
       TRELAXORB=0.1

       DO 9 KK=1,NTRAC-1
          DO 8 JJ=1,JG
             DO 7 LL=1,NL
                DO 6 II=1,IGC
                   TRAPRV(II,LL,JJ,KK)=-1.0
    6           CONTINUE
    7        CONTINUE
    8     CONTINUE
    9  CONTINUE

       READ (7,INVARPARAM)
       WRITE(2,INVARPARAM)             

       READ (7,INMAG)
       WRITE(2,INMAG)
       READ(7,INBINVAL)
       WRITE(2,INBINVAL)
       
       TAULIMIT=2.*LOG(1.e-2)/1.66/(10.**(-1.*OOM_IN/NL)-10.**(OOM_IN/NL))

       ! Isaac is guessing
       !TAULIMIT=100

      END                                                                 
