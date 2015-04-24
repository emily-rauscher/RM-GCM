C**********************************************************               
C     SUBROUTINE XSECT3                                            
C**********************************************************               
      SUBROUTINE XSECT3                                                   
C                                                                         
C     Outputs a file (fort.29) with KE spectrum
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
C                                                                         
C ER Modif
      PARAMETER(NWP=1+MM/MOCT,MOCTP=MOCT+1)
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
C     Array ordering in SPECTR must correspond to that in GRIDP.          
C                                                                         
      COMMON/SPECTR/Z(IGB),D(IGB),T(IGB),TRA(IGB,NTRAC),SP(IGA),GS(IGA)   
     :              ,SPA(IGA),VP(IGA),DTE(IGB),TT(IGB)                    
     :              ,TRAT(IGB,NTRAC),DT(IGB),ZT(IGB)                      
     :              ,ZMI(IGB),DMI(IGB),TMI(IGB)                           
     :              ,TRAMI(IGB,NTRAC),SPMI(IGA)                           
      COMPLEX Z,D,T,TRA,SP,GS,SPA,VP,DTE,TT,TRAT,DT,ZT                    
     :       ,ZMI,DMI,TMI,TRAMI,SPMI                                      
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
C                                                                         
C                                                                         
      REAL RKE(NNP,NLP,2)  ! (,,1): total, (,,2): zonal waveno.
      COMPLEX ZL(NNP,NWP),DL(NNP,NWP)                                                        !!SPEC.4
C
C      CHARACTER IXLAB(2)*22                                !!NSPECT.27
C      DATA IXLAB/'TOTAL WAVENUMBER (N)*.'                                  !!NSPECT.32
C     :          ,'ZONAL WAVENUMBER (M)*.'/                                 !!NSPECT.33
C                                                                          !!NSPECT.34
 6000 FORMAT(' PROCESSING KE M & N SPECTRA  AT DAY',F10.2,' KOUNT=',I10)   !!NSPECT.35
C 6010 FORMAT(' KINETIC ENERGY SPECTRUM AS FUNCTION OF ',A20)              !!NSPECT.37
 6010 FORMAT(3I4)
 6020 FORMAT(1X,I4,11(E11.4))                                           !!NSPECT.38
C                                                                         
C     Calculate KE spectra for current analysis time.                      !!NSPECT.45
C                                                                          !!NSPECT.46
      DO 10 II=1,2                                                         !!NSPECT.53
      DO 10 L=1,NLP                                                        !!NSPECT.54
      DO 10 N=2,NNP                                                       !!NSPECT.55
   10 RKE(N,L,II)=0.                                                       !!NSPECT.56
C                                                                          !!NSPECT.57
      MP=0                                                                 !!NSPECT.58
      MFPP=MM+1
      DO 20 M=1,MFPP,MOCT                                                  !!NSPECT.59
      MP=MP+1                                                              !!NSPECT.60
      DO 20 JP=M,NNP                                                      !!NSPECT.61
      ZL(JP,MP)=(0.,0.)                                                    !!NSPECT.62
      DL(JP,MP)=(0.,0.)                                                    !!NSPECT.63
   20 CONTINUE                                                             !!NSPECT.64
C                                                                          !!NSPECT.65
      IOLD=0                                                               !!NSPECT.66
      DO 100 L=1,NL                                                        !!NSPECT.67
C                                                                          !!NSPECT.68
      DO 30 IHEM=1,NHEM                                                    !!NSPECT.69
      MP=0                                                                 !!NSPECT.70
      DO 30 M=1,MFP,MOCT                                                   !!NSPECT.71
      MP=MP+1                                                              !!NSPECT.72
      DO 30 JP=M,NFP,MH                                                    !!NSPECT.73
      IOLD=IOLD+1                                                          !!NSPECT.74
      DL(JP+IHEM-1,MP)=D(IOLD)                                             !!NSPECT.75
      ZL(JP+2-IHEM,MP)=Z(IOLD)                                             !!NSPECT.76
   30 CONTINUE                                                             !!NSPECT.77
C                                                                          !!NSPECT.78
      DO 50 N=2,NNP                                                       !!NSPECT.79
      IF (MOCTP.GT.N) GOTO 50                                              !!NSPECT.80
      MP=1                                                                 !!NSPECT.81
      MMAX=MIN0(N,MFP)                                                     !!NSPECT.82
      DO 40 M=MOCTP,MMAX,MOCT                                              !!NSPECT.83
      MP=MP+1                                                              !!NSPECT.84
   40 RKE(N,L,1)=RKE(N,L,1)+(ZL(N,MP)*CONJG(ZL(N,MP))                      !!NSPECT.85
     :                      +DL(N,MP)*CONJG(DL(N,MP)))/SQ(N)               !!NSPECT.86
   50 CONTINUE                                                             !!NSPECT.87
C                                                                          !!NSPECT.88
      MN=1                                                                 !!NSPECT.89
      DO 60 MP=MOCTP,MFP,MOCT                                              !!NSPECT.90
      MN=MN+1                                                              !!NSPECT.91
      DO 60 N=MP,NNP                                                      !!NSPECT.92
   60 RKE(MN,L,2)=RKE(MN,L,2)+(ZL(N,MN)*CONJG(ZL(N,MN))                    !!NSPECT.93
     :                        +DL(N,MN)*CONJG(DL(N,MN)))/SQ(N)             !!NSPECT.94
C                                                                          !!NSPECT.95
  100 CONTINUE                                                             !!NSPECT.96

      REWIND 29
C                                                                          !!NSPECT.115
C     Loop over spectra.  Set wavenumber limits, dimensionalise,           !!NSPECT.116
C     calculate vertical integrals (make approximation surface             !!NSPECT.117
C     pressure equal to P0)
C     Print values
C                                                                          !!NSPECT.120
      WRITE(29,6010) NNP-1,MN-1,NLP

      DO 150 II=1,2                                                        !!NSPECT.121
C                                                                          !!NSPECT.122
      IF (II.EQ.1) THEN                                                    !!NSPECT.123
         MMAX=NNP                                                         !!NSPECT.124
         IF (2*(NN/2).EQ.NN.AND.2*(MOCT/2).EQ.MOCT) MMAX=NFP               !!NSPECT.125
         NPTS=MMAX-MOCT                                                    !!NSPECT.126
         MIN=MOCT                                                          !!NSPECT.127
         MAX=MMAX-1                                                        !!NSPECT.128
         INC=1                                                             !!NSPECT.129
         MINE=MIN+1                                                        !!NSPECT.130
         MAXE=MAX+1                                                        !!NSPECT.131
         INCE=INC                                                          !!NSPECT.132
      ELSE                                                                 !!NSPECT.133
         MMAX=MN                                                           !!NSPECT.134
         NPTS=MMAX-1                                                       !!NSPECT.135
         MIN=MOCT                                                          !!NSPECT.136
         MAX=MF                                                            !!NSPECT.137
         INC=MOCT                                                          !!NSPECT.138
         MINE=2                                                            !!NSPECT.139
         MAXE=MMAX                                                         !!NSPECT.140
         INCE=1                                                            !!NSPECT.141
      ENDIF                                                                !!NSPECT.142
C                                                                          !!NSPECT.143
      RKEFAC=CG*P0/(2.*GA)                                                 !!INITAL.371
C
      DO 110 L=1,NL                                                        !!NSPECT.144
      DO 110 N=2,MMAX                                                      !!NSPECT.145
  110 RKE(N,L,II)=RKE(N,L,II)*RKEFAC                                       !!NSPECT.146
C                                                                          !!NSPECT.147
      DO 120 L=1,NL                                                        !!NSPECT.148
      DO 120 N=2,MMAX                                                      !!NSPECT.149
  120 RKE(N,NLP,II)=RKE(N,NLP,II)+RKE(N,L,II)*DSIGMA(L)                    !!NSPECT.150
C                                                                          !!NSPECT.152
C      WRITE(29,6010) IXLAB(II)                                    !!NSPECT.153
C     Will print spectrum as total, then zonal waveno.
C     (last value is vertically integrated, L=nlp)
      DO 130 N=2,MMAX                                                      !!NSPECT.154
      NM=(N-1)*INC                                                         !!NSPECT.155
C      WRITE(29,6020) KOUNT,NM,(RKE(N,L,II),L=1,NLP)                      !!NSPECT.156
      WRITE(29,6020) NM,(RKE(N,L,II),L=1,NLP)                      !!NSPECT.156
  130 CONTINUE                                                             !!NSPECT.157
C                                                                          !!NSPECT.158
  150 CONTINUE                                                             !!NSPECT.166

         WRITE(29,'(///)')                                               !!NSPECT.48
         WRITE(29,6000) DAY,KOUNT                                        !!NSPECT.49

         END                                                                 
