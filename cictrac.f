C**********************************************************               
C             SUBROUTINE ICTRAC                                           
C**********************************************************               
      SUBROUTINE ICTRAC                                                   
C                                                                         
C     *********PMF Version 1.0 (8.4.97)*******************                
C     Subroutine which initialises the tracer fields.                     
C     Tracers can be started as any function of the usual grid            
C     point fields.                                                       
C     Mass mixing ratios converted to volume mixing ratios for            
C     flux form of advection.                                             
C                                    John Methven, 10.04.95               
C     Slighly modified to add more tracers and                            
C                 add array GWORK Piers Forster 22.01.97                  
C                                                                         
C     Changed again to combine with phys2 update                          
C     to include surface water and temperature                            
C     initilization.                                                      
C     Also lets water vapour profile be cacluated by setting              
C     surface RH     Piers Forster  18.02.97                              
C                                                                         
C      Flag LFLUX added for FLUX form of advection                        
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
      PARAMETER(NLTRG=NLTR*NHEM,NITRWG=7*NL*NHEM                          
     :         ,IDDAGFL=NITRWG*MGPP,IDPLG=3*IGC)                          
C                                                                         
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
     +              ,LFLUX                                                
     +              ,LBALAN,LRESTIJ                                       
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
     +              ,LNOISE,NFP                                               
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
C                                                                         
C     Constants and arrays needed for the fast Fourier transforms         
C                                                                         
      COMMON/COMFFT/NTGW,NRSTGW,NTWG,NRSTWG,NTGD,NRSTGD                   
     +              ,TRIG(IDA),WORK(IDF),IFAX(10)                         
C                                                                         
C                                                                         
C     Polynomial used to aid vectorization of Legendre transforms         
C                                                                         
      COMMON/POLYNO/POLY(NWJ2,2,4),CMPA(IGL)                              
      COMPLEX CMPA                                                        
C                                                                         
C                                                                         
C     Array ordering in GRIDP must correspond to that in SPECTR.          
C     Real arrays: multi-level arrays are 1-dimensional.                  
C                                                                         
      COMMON/GRIDP/ CHIG(IGD),SFG(IGD),UG(IGD),VG(IGD)                    
     *              ,ZG(IGD),DG(IGD),TG(IGD)                              
     +              ,TRAG(IGD,NTRAC)                                      
     *              ,PLG(IGC),PJG(IGC),PMG(IGC)                           
     *              ,SPG(IGC),VPG(IGC),EG(IGD)                            
     +              ,TNLG(IGD),TRANLG(IGD,NTRAC),FUG(IGD),FVG(IGD)        
     +              ,UTG(IGD),UTRAG(IGD,NTRAC)                            
     +              ,VTG(IGD),VTRAG(IGD,NTRAC),FVGT(IGD),FUGT(IGD)        
     $              ,GRPAD(NGRPAD)                                        
C                                                                         
C                                                                         
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
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
     : ,ctqi,sdsn,shcs,shcsp,shcsn,skse,sksn,slhf,sd1,sd2,sdw             
     :        ,ssmc,sdsnd,sasnow,saice,shsstar,shsmax
     :     ,LOC,LNOICE,LOLDBL,LCOND,LNNSK
     :              ,NLCR,CURHM,AKTC,AKQC,CUBMT,CBADJT,CBADJP
     :              ,SKAP(NL),SK(NLM),FWS(NL),CLR(NL),FB(NLM)             
     :              ,TTDC(NL),QTDC(NL),TTMC(NL),QTMC(NL),TC(NL),QC(NL)    
     :              ,CTCR(NL,NHEM),CTLR(NL,NHEM)
     :              ,LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ
     :              ,LSL,NAVRD,NAVWT,DELT2C,SHCO,SHCI,ITSLL,ITSLO,NCUTOP
                LOGICAL LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ,LSL,LOC                    
     :       ,LNOICE,LOLDBL,LCOND,LNNSK                                   
C                                                                         
      COMMON/COMGRM/DUDLSG(IGC,NL),DVDLSG(IGC,NL),DTDLSG(IGC,NL)          
C                                                                         
      LOGICAL LNSURF,LRH                                                  
      REAL RH(NL),TDUM(JG),QDUM(JG)                                       
      COMPLEX GWORK(IGL,(3+NTRAC)*NL+2)                                   
      REAL DAGFL(IDDAGFL),DAGPLG(IDPLG),TEMPTR(IGD),DAFTR(NLTRG*MGPP)     
      EQUIVALENCE (DAGFL(1),CHIG(1)),(DAGPLG(1),PLG(1))                   
     :           ,(DAFTR(1),TRAG(1,1))                                    
      CHARACTER LHEM(2)*25                                                
      DATA LHEM/'  N HEM  EQUATOR.....POLE'                               
     :         ,'  S HEM  EQUATOR.....POLE'/                              
C                                                                         
      NAMELIST/INQ/LRH,RH,LNSURF                                          
C                                                                         
 2010 FORMAT('0SURFACE TEMPERATURE IN DEGC',A25/20(1X,F8.2))              
 2020 FORMAT('0SURFACE SPECIFIC HUMIDITY IN G/KG',A25/20(1X,f8.4))        
 2030 FORMAT(3X,20(1X,F5.1))                                              
C                                                                         
      LRH=.FALSE.                                                         
      LNSURF=.FALSE.                                                      
      DO L=1,NL                                                           
        RH(L)=90.0*SIGMA(L)                                               
      ENDDO                                                               
C      write(*,*) 'now reading fort.7 in ictrac'
      READ(7,INQ)                                                         
      WRITE(2,INQ)                                                        
C                                                                         
      IF (.NOT.LRH.AND.NTRACO.LT.NTRAC) THEN                              
         WRITE(2,'(A25)') 'No water Vapour. Need to set LRH to TRUE'      
         CALL ABORT                                                       
      ENDIF               

C convert RH from % to fraction                                           
      DO L=1,NL                                                           
        RH(L)=RH(L)/100.0                                                 
      ENDDO                                                               
C     Surface characteristics now read in in INISURF, not here            
C                                                                         
C     Save old tracers in the array TRAT and reset TRA to zero.           
C                                                                         
      IF (NTRACO.GT.0) THEN                                               
         DO 10 KK=1,NTRACO                                                
            DO 20 I=1,IGB                                                 
               TRAT(I,KK)=TRA(I,KK)                                       
               TRA(I,KK)=0.0                                              
 20         CONTINUE                                                      
 10      CONTINUE                                                         
      ENDIF                                                               
C                                                                         
C     Loop over latitude                                                  
C                                                                         
      JL=1                                                                
      IF(JGL.EQ.1) REWIND 25                                              
      DO 400 IH=1,JG                                                      
         JH=IH                                                            
         IF(JGL.EQ.1) READ(25) ALP,DALP,DP,DQ                             
C                                                                         
C     Spectral to grid transforms for initial model variables.            
C     Calculate the horizontal derivatives of temperature as well         
C     as the usual derivatives done by LTI.                               
C                                                                         
         CALL LTI                                                         
         CALL LTIDT                                                       
         NTALL=(NITRWG-1)/NCRAY                                           
         NRSTALL=NITRWG-NCRAY*NTALL                                       
C                                                                         
C  Sets up surface temperature here                                       
C                                                                         
C                                                                         
         IF(NTALL.EQ.0) GOTO 40                                           
         DO 30 I=1,NTALL                                                  
 30      CALL FFT991(DAGFL(1+(I-1)*NCRAY*MGPP),WORK,TRIG,IFAX             
     :           ,1,MGPP,MG,NCRAY,1)                                      
 40      CALL FFT991(DAGFL(1+NTALL*NCRAY*MGPP),WORK,TRIG,IFAX             
     :        ,1,MGPP,MG,NRSTALL,1)                                       
         CALL FFT991(DAGPLG(1),WORK,TRIG,IFAX,1,MGPP,MG,3*NHEM,1)         

C                                                                         
C     Grid point calculations for surface pressure,                       
C     temperature T (=TO+T') and vertical derivatives of U,V and T.       
C                                                                         
         DO 50 I=1,IGC                                                    
            PLG(I)=EXP(PLG(I))                                            
 50      CONTINUE                                                         
         K=0                                                              
         DO 60 L=1,NL                                                     
            DO 65 I=1,IGC                                                 
               K=K+1                                                      
               TG(K)=TG(K)+T0(L)                                          
               IF (LNSURF.AND.L.EQ.NL) THEN                                        
                  TSTAR(I,IH)=TG(K)/SKAP(NL)                                          
               ENDIF                                                               
 65         CONTINUE                                                      
 60      CONTINUE                                                         
C                                                                         
         CALL DLSGCR(IGC,UG,RGG,DUDLSG,IGC,NL)                            
         CALL DLSGCR(IGC,VG,RGG,DVDLSG,IGC,NL)                            
         CALL DLSGCR(IGC,TG,RGG,DTDLSG,IGC,NL)                            
C                                                                         
C     Grid point calculations for new tracer fields at one                
C     latitude (IH).                                                      
C                                                                         
C initilise water vapour if LRH set  depends on RH                        
C       goes into tracer No. 1                                            
         IOFM=0                                                           
         DO 810 IHEM=1,NHEM                                               
               SPZ=SSUM(MG,PLG(1+IOFM),1)/REAL(MG)                        
            IF (LRH) THEN                                                 
               DO I=1,MG                                                  
                  J=I+IOFM                                                
                  ESCON=1./PLG(J)                                         
                  K=J                                                     
                  DO L=1,NL                                               
                   TRAG(K,1)=ESCON*RH(L)*PQSAT(TG(K))/SIGMA(L)            
                   K=K+IGC                                                
                  ENDDO                                                   
               ENDDO                                                      
            ENDIF                                                         
            IF (LNSURF) THEN                                                    
               ESCON=1./SPZ                                                        
               DO K=1,IGC                                                          
                  QSTAR(K,IH)=ESCON*PQSAT(TSTAR(K,IH))                                
CCCC      ESCON=1./SPZ                                                    
CCCC      QSTAR(J,JH)=ESCON*PQSAT(TSTAR(J,JH))                            
               ENDDO                                                               
             PRINT '(i2,x,10(F5.1,x))',ih,(ct*tstar(k,ih)-273.15,k=1,10)         
            ENDIF                                                               
            IOFM=MGPP                                                     
 810     CONTINUE                                                         
                                                                          
         DO 70 KK=NTRACO+1,NTRAC                                          
            DO 300 I=1,IGD                                                
               TEMPTR(I)=TRAG(I,KK)                                       
 300        CONTINUE                                                      
            IF(KOLOUR(KK).EQ.1) THEN                                      
               CALL PVCR(TEMPTR,1)                                        
            ENDIF                                                         
            IF(KOLOUR(KK).EQ.2) THEN                                      
               CALL PVCR(TEMPTR,2)                                        
            ENDIF                                                         
            DO 310 I=1,IGD                                                
               TRAG(I,KK)=TEMPTR(I)                                       
 310        CONTINUE                                                      
C                                                                         
C     Convert mass mixing ratios to volume mixing ratios.                 
C     If LFLUX set                                                        
            IF (LFLUX) THEN                                               
            DO 95 L=1,NL                                                  
               LL=(L-1)*IGC                                               
               DO 94 I=1,IGC                                              
                  TRAG(LL+I,KK)=PLG(I)*TRAG(LL+I,KK)                      
 94            CONTINUE                                                   
 95         CONTINUE                                                      
            ENDIF                                                         
 70      CONTINUE                                                         
C                                                                         
C     Transform all tracer fields to spectral space at each latitude.     
C     First do inverse FFT to give half-transforms and then call          
C     HANAL for direct Legendre transform.                                
C                                                                         
         NTTR=(NLTRG-1)/NCRAY                                             
         NRSTTR=NLTRG-NCRAY*NTTR                                          
         IF(NTTR.EQ.0) GO TO 110                                          
         DO 100 I=1,NTTR                                                  
 100        CALL FFT991(DAFTR(1+(I-1)*NCRAY*MGPP),WORK,TRIG,IFAX          
     :           ,1,MGPP,MG,NCRAY,-1)                                     
 110     CALL FFT991(DAFTR(1+ NTTR*NCRAY*MGPP),WORK,TRIG,IFAX             
     :        ,1,MGPP,MG,NRSTTR,-1)                                       
         CALL HANAL(TRAG,GWORK,TRA,NLTR,2)                                
C                                                                         
C     End of latitude loop. Note that if JGL=JG then JINC=1               
C          if JGL=1  then JINC=0                                          
C     Therefore JL=JH(=IH) if JGL=JG, or JL=1 if JGL=1.                   
C                                                                         
         JL=JL+JINC                                                       
 400  CONTINUE                                                            
C                                                                         
C     Return old tracers from the restart record into array TRA and       
C     set the t-1 record for the new tracers.                             
C                                                                         
      IF (NTRACO.GT.0) THEN                                               
         DO 410 KK=1,NTRACO                                               
            DO 420 I=1,IGB                                                
               TRA(I,KK)=TRAT(I,KK)                                       
               TRAT(I,KK)=0.0                                             
 420        CONTINUE                                                      
 410     CONTINUE                                                         
      ENDIF                                                               
      DO 430 KK=NTRACO+1,NTRAC                                            
         DO 440 I=1,IGB                                                   
            TRAMI(I,KK)=TRA(I,KK)                                         
 440     CONTINUE                                                         
 430  CONTINUE                                                            
C                                                                         
C    writes out surface temp and humdidity to channel 14                  
C     Continue to write out restart in old format....                     
      if (LNSURF) WRITE(14,*) ((tstar(i,j),j=1,jg),i=1,igc),              
     $     ((qstar(i,j),j=1,jg),i=1,igc)                                  
      WRITE(2,'(A40)') 'O DEG LONG X-SECTS SHOWN BELOW'                   
                                                                          
      DO 820 IHEM=1,NHEM                                                  
      DO 510 J=1,JG                                                       
      TDUM(J)=CT*TSTAR(IHEM,J)-273.15                                     
 510  QDUM(J)=CQ*QSTAR(IHEM,J)                                            
      IF(JG.LE.20) THEN                                                   
        DO 520 NCHAN=2,2                                                  
        WRITE(NCHAN,2010) LHEM(IHEM),(TDUM(J),J=JG,1,-1)                  
  520   WRITE(NCHAN,2020) LHEM(IHEM),(QDUM(J),J=JG,1,-1)                  
      ELSE                                                                
        DO 530 NCHAN=2,2                                                  
        WRITE(NCHAN,2010) LHEM(IHEM),(TDUM(J),J=JG,2,-2)                  
        WRITE(NCHAN,2030)            (TDUM(J),J=JG-1,1,-2)                
        WRITE(NCHAN,2020) LHEM(IHEM),(QDUM(J),J=JG,2,-2)                  
  530   WRITE(NCHAN,2030)            (QDUM(J),J=JG-1,1,-2)                
      ENDIF                                                               
  820 CONTINUE                                                            
      RETURN                                                              
      END                                                                 
