C**********************************************************               
C             SUBROUTINE XSECT                                            
C**********************************************************               
      SUBROUTINE XSECT(ISKIP)                                             
C                                                                         
C     This subroutine gives quick look sigma-latitude x-sections          
C     of [U], [T] and [V*T*].                                             
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
C                                                                         
C     Switches counters and constants controlling type and frequency of   
C     model output                                                        
C                                                                         
      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC                        
     +              ,RNTAPO                                               
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR                          
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY                          
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG                 
     +              ,LSPO(NL),LGPO(NL)                                    
     $              ,LSHIST,LMINIH                                        
      LOGICAL LSHIST,LMINIH                                               
      LOGICAL LSPO,LGPO                                                   
C                                                                         
      COMMON/GRIDSS/ZG1(IGD,JG),DG1(IGD,JG),UG1(IGD,JG),VG1(IGD,JG),      
     :              TG1(IGD,JG),SPG1(IGC,JG)                              
C                                                                         
C                                                                         
C     Restoration fields and timescale                                    
C                                                                         
      COMMON/RESTOR/ZRES(IGN),DRES(IGN),TRES(IGN),SPRES(IGM),DAMP         
C                                                                         
      INTEGER LU(IGG),LT(IGG),LVT(IGG)                                    
c                                                                         
      integer nmaxdims,nall                                               
      parameter(nmaxdims=4,nall=10)                                       
      integer ndim,nvar,natts(nall),nattsvar(nall),                       
     :        vdims(nall),vadims(nmaxdims,nall),                          
     :        ndims(nall)                                                 
      character*200 dimname(nall),varname(nall),                          
     :              attdimname(2,nmaxdims,nall),                          
     :              attvarname(2,nmaxdims,nall)                           
      integer nc1,iddim1(nall),idvar1(nall)                               
      common/netcdfvars/nc1,iddim1,idvar1                                 
c                                                                         
      REAL LUG(MG,JGG,NL),LTG(MG,JGG,NL)                                  
      LOGICAL LFIRST                                                      
      INTEGER INETCOUNT                                                   
      SAVE LFIRST,INETCOUNT                                               
      DATA LFIRST/ .TRUE. /                                               
      DATA INETCOUNT/0/                                                   
C                                                                         
C     Set ISLT, The default being to have 16 I5 integers per level.       
C                                                                         
      IF (ISKIP.EQ.0) THEN                                                
         RETURN                                                           
      ELSE IF (ISKIP.GT.0) THEN                                           
         ISLT = ISKIP                                                     
      ELSE                                                                
         ISLT = (NHEM*JG-1)/18 + 1                                        
      ENDIF                                                               
      IF (LFIRST) THEN                                                    
CC         call ini_netcdf                                                
         LFIRST=.FALSE.                                                   
      END IF                                                              
      INETCOUNT=INETCOUNT+1                                               
C                                                                         
C     Output is wanted. Read grid point fields from stream 24             
C                                                                         
C      REWIND 24                                                          
      RMG=1./REAL(MG)                                                     
      WRITE (2,104) DAY,DOY                                               
  104 FORMAT(/' CROSS SECTIONS FOR DAY, DOY',2F6.1)                       
      DO 10 J=1,JG                                                        
C                                                                         
C        Loop over latitudes to read gridpoint data, dimensionalise       
C        and calculate zonal means                                        
C                                                                         
         SEC=ALAT(J)/57.29578                                             
         SEC=10.*CV/COS(SEC)                                              
C         READ(24) ZG,DG,UG,VG,TG,SPG                                     
         DO I=1,IGD                                                       
            ZG(I)=ZG1(I,J)                                                
            DG(I)=DG1(I,J)                                                
            UG(I)=UG1(I,J)                                                
            VG(I)=VG1(I,J)                                                
            TG(I)=TG1(I,J)                                                
         END DO                                                           
         DO I=1,IGC                                                       
            SPG(I)=SPG1(I,J)                                              
         END DO                                                           
C                                                                         
C        Loop for [U] and [T] sections                                    
C                                                                         
         DO L=1,NL                                                        
         DO IHEM=1,NHEM                                                   
               IF (IHEM.EQ.1) THEN                                        
                  JJ=J                                                    
               ELSE                                                       
                  JJ=JGGP-J                                               
               ENDIF                                                      
            IOF=MGPP*(IHEM-1)+(L-1)*IGC                                   
            DO I=1,MG                                                     
               LUG(I,JJ,L)=0.1*SEC*UG(I+IOF)                              
            END DO                                                        
            DO I=1,MG                                                     
               LTG(I,JJ,L)=CT*(TG(I+IOF)+T0(L))                           
            END DO                                                        
        end do                                                            
        end do                                                            
         DO 11 L=1,NL                                                     
            DO 12 IHEM=1,NHEM                                             
               IP=(L-1)*IGC+(IHEM-1)*MGPP                                 
               UB=0.                                                      
               DO 30 I=1,MG                                               
                  UB=UB+UG(IP+I)                                          
   30          CONTINUE                                                   
               UB=UB*RMG                                                  
               TB=0.                                                      
               DO 40 I=1,MG                                               
                  TB=TB+TG(IP+I)                                          
   40          CONTINUE                                                   
               TB=TB*RMG+T0(L)                                            
C                                                                         
C              Check whether hemisphere or globe                          
C                                                                         
               IF (IHEM.EQ.1) THEN                                        
                  JJ=J                                                    
               ELSE                                                       
                  JJ=JGGP-J                                               
               ENDIF                                                      
               IX=(L-1)*JGG+JJ                                            
               LU(IX)=NINT(UB*SEC)                                        
               LT(IX)=NINT(TB*CT)                                         
   12       CONTINUE                                                      
   11    CONTINUE                                                         
C                                                                         
C        Loop for temperature flux section                                
C                                                                         
         DO 13 L=1,NL                                                     
            DO 14 IHEM=1,NHEM                                             
               IP=(L-1)*IGC+(IHEM-1)*MGPP                                 
               VT=0.                                                      
               DO 15 I=1,MG                                               
                  VT=VT+TG(IP+I)*VG(IP+I)                                 
   15          CONTINUE                                                   
               VT=RMG*VT                                                  
               TB=0.                                                      
               DO 50 I=1,MG                                               
                  TB=TB+TG(IP+I)                                          
   50          CONTINUE                                                   
               TB=TB*RMG                                                  
               VB=0.                                                      
               DO 60 I=1,MG                                               
                  VB=VB+VG(IP+I)                                          
   60          CONTINUE                                                   
               VB=VB*RMG                                                  
C                                                                         
C              Check whether hemisphere or globe                          
C                                                                         
               IF (IHEM.EQ.1) THEN                                        
                  JJ=J                                                    
               ELSE                                                       
                  JJ=JGGP-J                                               
               ENDIF                                                      
               IX=(L-1)*JGG+JJ                                            
               LVT(IX)=NINT((VT-VB*TB)*SEC*CT)                            
   14       CONTINUE                                                      
   13    CONTINUE                                                         
   10 CONTINUE                                                            
C                                                                         
CC      call writevar2(nc1,idvar1(1),lug,                                 
CC     :                     1,mg,1,jgg,1,nl,inetcount,inetcount)         
CC      call writevar2(nc1,idvar1(2),ltg,                                 
CC     :                     1,mg,1,jgg,1,nl,inetcount,inetcount)         
C     Printing sections - separate loop for each field.                   
C     First the zonal wind                                                
C                                                                         
      WRITE (2,100)                                                       
  100 FORMAT(' ZONAL WIND IN 0.1 m/s')                                    
      DO 20 L=1,NL                                                        
         J0=(L-1)*JGG                                                     
         WRITE (2,101) (LU(J0+J),J=1,JGG,ISLT)                            
   20 CONTINUE                                                            
101   FORMAT(18I4)                                                        
C                                                                         
C     Secondly the temperature                                            
C                                                                         
      WRITE(2,102)                                                        
  102 FORMAT(' TEMPERATURE IN K')                                         
      DO 21 L=1,NL                                                        
         J0=(L-1)*JGG                                                     
         WRITE(2,101) (LT(J0+J),J=1,JGG,ISLT)                             
   21 CONTINUE                                                            
C                                                                         
C     Thirdly temperature flux                                            
C                                                                         
      WRITE(2,103)                                                        
  103 FORMAT(' POLEWARD TEMPERATURE FLUX IN 0.1Km/s')                     
      DO 22 L=1,NL                                                        
         J0=(L-1)*JGG                                                     
         WRITE(2,101) (LVT(J0+J),J=1,JGG,ISLT)                            
   22 CONTINUE                                                            
C                                                                         
C      REWIND 24                                                          
      RETURN                                                              
      END                                                                 
                                                                          
