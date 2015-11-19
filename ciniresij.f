C*************************************************************            
C               SUBROUTINE INIRESIJ                                       
C*************************************************************            
      SUBROUTINE INIRESIJ                                                 
C                                                                         
C     Sets up restoration variables and arrays. Sets NAMELIST             
C     variables to their default settings, then reads NAMELIST            
C     Sets up Rayleigh friction coefficients.                             
C                                                                         
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'
C      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
C     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
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
C     Restoration temperature field and constants which determine it,     
C     also contains timescales                                            
C                                                                         
      COMMON/RESTIJ/TTRES(IGB)                                            
     + ,DTNS,DTEP,DTTRP,FAC(NL),DDAMP(NL),TFRC(NL),YRLEN,TRS(NL)          
     +  ,RESTTT(NL),REDTEP(NL)
     +  ,ALR,ZTROP,TGR                                                    
      COMPLEX TTRES                                                       
C                                                                         
C     Basic planetary parameters for run plus information about           
C     vertical grid structure                                             
C                                                                         
C     Note that RD and GASCON are identical and CPD is set from RD,AKAP.  
      COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)                
     +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)     
     +              ,TKP(NL),C(NL2),SQH(NNP)                              
     +              ,MF,MFP,JZF,NF,NFP                                    
     +              ,AKAP,GA,GASCON,RADEA,WW,P0,PFAC,EZ,AIOCT             
     +              ,RD,RV,CPD,CLATNT                                     
     +              ,LRSTRT,LSHORT,LTVEC,LSTRETCH                         
     +              ,LFLUX                                                
     +              ,LBALAN,LRESTIJ                                       
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
     +              ,LNOISE                                               
      COMPLEX EZ,AIOCT                                                    
      LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ                 
     +       ,LFLUX,LNOISE                                                
     +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                           
C                                                                         
C                                                                         
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BM1(IDE),AK(NNP),AQ(NL2),G(NL2),TAU(NL2)              
     +              ,KOUNT,KITS,KSTART,KTOTAL,KRUN,BEGDAY,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,CTRA(NTRAC),KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C                                                                         
      DIMENSION RESTIM(NL)                                                
C                                                                         
      NAMELIST/INPRSIJ/ TFRC,RESTIM,RESTTT,REDTEP,DTNS,DTEP,ALR,                        
     +                  DTTRP,ZTROP,TGR,YRLEN                             
C                                                                         
      TFRC(NL)=0.1                                                         
      DO 19 L=1,NL-1                                                      
         TFRC(L) = 0.1! Rayleigh friction time scale IN DAYS               
! NB zero here equates to infinity - no friction.                         
 19   CONTINUE                                                            
      DO 20 L=1,NL                                                        
         RESTIM(L) = 15.0e28                      
         RESTTT(L) = 0.
         REDTEP(L) = 0.
 20   CONTINUE                                                            
      DTNS=0.                                                             
      DTEP=60.                                                            
      ALR=6.5E-03                                                         
      DTTRP=2.                                                            
      ZTROP=12.0E03                                                       
      TGR=288.                                                            
      YRLEN=0.                                                            
C                                                                         
C      write(*,*) 'reading fort.7 in iniresij'
      READ(7,INPRSIJ)                                                    
      WRITE(2,INPRSIJ)                                                    
C                                                                         
C     Dimensionless coefficient for Newtonian cooling friction            
C     and timestep. A day is 2*pi in non dimensional                      
C     units using omega as the unit of frequency.                         
C                                                                         
      DO 22 L=1,NL                                                        
         IF (RESTIM(L).GT.0.0) THEN                                       
            DDAMP(L)=1.0/(PI2*RESTIM(L))                                  
         ELSE                                                             
            DDAMP(L)=0.0                                                  
         ENDIF                                                            
         IF (TFRC(L).GT.0.0) THEN                                         
            TFRC(L)=1.0/(PI2*TFRC(L))                                     
         ELSE                                                             
            TFRC(L)=0.0                                                   
         ENDIF                                                            
 22   CONTINUE                                                            
C                                                                         
C     Make temperatures dimensionless                                     
C                                                                         
      DTNS=DTNS/CT                                                        
      DTEP=DTEP/CT                                                        
      DTTRP=DTTRP/CT                                                      
C     ER modif
      DO 23 L=1,NL
         REDTEP(L)=REDTEP(L)/CT
 23   CONTINUE
C                                                                         
C     Loop to set array FAC - this controls temperature gradients         
C     as a function of SIGMA in TTRES. It is a sine wave from one         
C     at SIGMA = 1 to zero at STPS (SIGMA at the tropopause).             
C                                                                         
C     First find SIGMA at ZTROP                                           
C                                                                         
      TTROP = TGR - ZTROP*ALR                                             
      STPS = (TTROP/TGR)**(GA/(ALR*GASCON))                               
      DO 600 L=1,NL                                                       
         THING=SIN(0.5*PI*(SIGMA(L)-STPS)/(1.-STPS))                      
         IF (THING.LT.0.) THEN                                            
            FAC(L)=0.                                                     
         ELSE                                                             
            FAC(L)=THING                                                  
         ENDIF                                                            
600   CONTINUE                                                            
C                                                                         
      END                                                                 
