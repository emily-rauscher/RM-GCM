C**************************************************************           
C                   SUBROUTINE INISTR                                     
C**************************************************************           
      SUBROUTINE INISTR                                                   
C                                                                         
C     Reads in data for a start/restart run                               
C                                                                         
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'
C      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
C     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
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
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
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
C                                                                         
C     Restoration fields and timescale                                    
C                                                                         
      COMMON/RESTOR/ZRES(IGN),DRES(IGN),TRES(IGN),SPRES(IGM),DAMP         
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
      COMMON/STATS/GMSP0,GMSPMI,LMASCOR,LMASOLD,LMASPRT                   
      LOGICAL LMASCOR,LMASOLD,LMASPRT                                     
C                                                                         
C-----------------------------------------------------------------------  
C                                                                         
C     Constants and arrays needed for balancing                           
C                                                                         
      COMMON/BALAN/BFILT(NL),RGT0(NL),RG(NL2),TMEAN(NL)                   
     +            ,EP1(IGA),EP2(IGA),KBAL,MFTBAL,SRGT0,LTBAL              
      LOGICAL LTBAL                                                       
C                                                                         
       COMMON/GSG/GSG(IGC,JG)                                             
       REAL htnet                                                         
       COMMON /RADHT/ HTNET(NHEM,JG,MG,NL)                                
       REAL TAVE(IGP)                                                     
C                                                                         
 2000 FORMAT(' ***WARNING*** KITS < 1 FOR AN INITIAL RUN.'/)              
 2010 FORMAT(/' ***ABORT*** THE HISTORY RECORDS READ FROM CHANNEL '       
     + ,I3,/' ARE NOT IN CORRECT FORMAT ')                                
 2020 FORMAT(/' ***ABORT*** THE RUN NUMBER IN THE HISTORY RECORD ',F10.3  
     +       /' DIFFERS FROM RNTAPO IN THE NAMELIST ',F10.3)              
 2030 FORMAT(/' ***ABORT*** CANNOT FIND THE CORRECT HISTORY RECORD.'/     
     + ' LOOKING FOR DAY',F8.2/,' BUT THE NEAREST RECORD FOUND',          
     + ' IS FOR DAY',F8.2)                                                
 2040 FORMAT(/' HISTORY RECORD READ FROM CHANNEL ',I3,/                   
     + ' KOUNT  RMTAPE  DAY  DOY =',I8,3F12.3)                            
 2011 FORMAT(/' ***ABORT*** THE RESTART RECORDS READ FROM CHANNEL '       
     + ,I3,/' ARE NOT IN CORRECT FORMAT ')                                
 2021 FORMAT(/' ***ABORT*** THE RUN NUMBER IN THE RESTART RECORD ',F10.3  
     +       /' DIFFERS FROM RNTAPO IN THE NAMELIST ',F10.3)              
 2031 FORMAT(/' ***ABORT*** CANNOT FIND THE CORRECT RESTART RECORD.'/     
     + ' LOOKING FOR DAY',F8.2/,' BUT THE NEAREST RECORD FOUND',          
     + ' IS FOR DAY',F8.2)                                                
 2041 FORMAT(/' RESTART RECORD READ FROM CHANNEL ',I3,/                   
     + ' KOUNT  RMTAPE  DAY  DOY =',I8,3F12.3)                            
 2012 FORMAT(/' ***ABORT*** THE RESTORATION RECORDS READ FROM CHANNEL '   
     + ,I3,/' ARE NOT IN CORRECT FORMAT ')                                
 2022 FORMAT(/' ***ABORT*** THE RUN NUMBER IN THE RESTORATION RECORD '    
     +       ,F10.3/' IS NOT THE SAME AS RNTAPO IN THE NAMELIST ',F10.3)  
 2032 FORMAT(/' ***ABORT*** CANNOT FIND THE CORRECT RESTORATION RECORD.'  
     + /' LOOKING FOR DAY',F8.2/,' BUT THE NEAREST RECORD FOUND',         
     + ' IS FOR DAY',F8.2)                                                
 2042 FORMAT(' RESTORATION RECORD READ FROM CHANNEL ',I3,/                
     + ' KOUNT  RMTAPE  DAY  DOY =',I8,3F12.3)                            
 2050 FORMAT(/' SPECTRAL ARRAYS ARE SET TO ZERO ')                        
 2100 FORMAT(/' MASS INFORMATION READ AT RESTART:'                        
     :       /'    INITIAL (REFERENCE) MASS     = ',1PE20.12,' (PA)'      
     :       /'    TIME-LAGGED MASS             = ',1PE20.12,' (PA)')     
 2110 FORMAT(/' ***ABORT: MASS RESTART RECORD HAS WRONG LENGTH')          
 2130 FORMAT(/' ***ABORT: EOF READING MASS RESTART RECORD AT DAY ',F8.2)  
C                                                                         
C     Initialize spectral arrays to zero and overwrite as desired         
C                                                                         
      DO 1 I=1,IGA                                                        
         SP(I)=(0.0,0.0)                                                  
         SPMI(I)=(0.,0.)                                                  
         GS(I)=(0.0,0.0)                                                  
    1 CONTINUE                                                            
      DO 5 I=1,IGB                                                        
         Z(I)=(0.0,0.0)                                                   
         D(I)=(0.0,0.0)                                                   
         T(I)=(0.0,0.0)                                                   
         ZMI(I)=(0.0,0.0)                                                 
         DMI(I)=(0.0,0.0)                                                 
         TMI(I)=(0.0,0.0)                                                 
    5 CONTINUE                                                            
      DO 6 KK=1,NTRAC                                                     
         DO 6 I=1,IGB                                                     
            TRA(I,KK)=(0.0,0.0)                                           
            TRAMI(I,KK)=(0.0,0.0)                                         
    6 CONTINUE                                                            
      DO 8 l=1,nl                                                         
        DO 8 i=1,mg                                                       
          DO 8 j=1,jg                                                     
            DO 8 ihem=1,nhem                                              
              htnet(ihem,j,i,l)=0.                                        
    8 CONTINUE                                                            
C                                                                         
C     Initialise current and global reference mass to zero (unset).       
C                                                                         
      GMSP0=0.0                                                           
      GMSPMI=0.0                                                          
C                                                                         
      IF (.NOT.LRSTRT) THEN                                               
C                                                                         
C        Initial run                                                      
C                                                                         
         IF ( KITS .LT. 1) THEN                                           
            WRITE(2,2000)                                                 
         ENDIF                                                            
         DAY=0.0                                                          
         IF (KITS.EQ.0) THEN                                              
            KTOTAL=KRUN                                                   
            KOUTP=0                                                       
            KOUTE=0                                                       
            KOUTH=0                                                       
            KOUTR=0                                                       
         ELSE                                                             
            KTOTAL=KRUN+KITS-1                                            
            KOUTP=1-KITS                                                  
            KOUTE=1-KITS                                                  
            KOUTH=1-KITS                                                  
            KOUTR=1-KITS                                                  
         END IF                                                           
C                                                                         
C        Initialise restoration array                                     
C                                                                         
         IF (LRESTIJ) CALL SETZT                                          
C                                                                         
C        Initialise spectral arrays                                       
C                                                                         
         IF (LBALAN) THEN                                                 
           CALL INIBAL                                                    
           KOUNT=-KBAL                                                    
           KSTART=0                                                       
         ELSE                                                             
           IF (LRESTIJ) THEN                                              
              CALL INISP                                                  
           ELSE                                                           
              WRITE (2,2050)                                              
           ENDIF                                                          
           KOUNT=0                                                        
           KSTART=KOUNT                                                   
         ENDIF                                                            
      ELSE                                                                
C                                                                         
C        Code for restart and normal mode perturbation runs.              
C        assume spectral data is set up non-dimensionalised               
C        on a history (LSHORT) or restart (.NOT.LSHORT) record.
C
         DAYNEAR=-9999.5 ! allows desired day no. = 0                     
         IF (LSHORT) THEN
            ID=10
  180       IF (NTRACO.EQ.0) THEN                                         
              READ(ID,END=1000)RKOUNT,RM1TAPE,DAY,Z,D,T,SP,RM2TAPE        
            ELSE                                                          
              READ(ID,END=1000)RKOUNT,RM1TAPE,DAY,DOY,Z,D,T               
     +             ,((TRA(I,KK),I=1,IGB),KK=1,NTRACO),SP,RM2TAPE          
            ENDIF                                                         
            IF (ABS(RM1TAPE-RM2TAPE) .GT. 1.0E-03) THEN
               WRITE(2,2010) ID                                           
               CALL ABORT                                                 
            ENDIF                                                         
            IF (ABS(RM1TAPE-RNTAPO) .GT. 1.0E-03) THEN
               WRITE(2,2020) RM1TAPE,RNTAPO                               
               CALL ABORT                                                 
            ENDIF
C Extra read if restart record contains mass correction.                  
            IF (LMASOLD) THEN                                             
               READ(ID,END=1010) RK,RMT1,RDAY,RREC,GMSP0,GMSPMI,RMT2      
               IF (RMT2.NE.RMT1) THEN                                     
                  WRITE(2,2110)                                           
                  CALL ABORT                                              
               ENDIF                                                      
            ENDIF


            IF (ABS(DAY-BEGDAY) .GT. 1.0E-02) THEN                        
               IF (ABS(DAY-BEGDAY) .LT. ABS(DAYNEAR-BEGDAY)) THEN         
                  DAYNEAR = DAY                                           
               ENDIF                                                      
               GOTO 180                                                   
            ELSE                                                          
               GOTO 200                                                   
            ENDIF                                                         
 1000       WRITE(2,2030) BEGDAY,DAYNEAR                                  
            CALL ABORT                                                    
 1010       WRITE(2,2130) DAY                                             
            CALL ABORT


C                                                                         
  200       KOUNT=NINT(RKOUNT)                                            
            WRITE(2,2040) ID,KOUNT,RM1TAPE,BEGDAY,DOY                     
            KOUNT=0                                                       
            IF (KITS.EQ.0) THEN                                           
               KTOTAL=KRUN                                                
               KOUTP=0                                                    
               KOUTE=0                                                    
               KOUTH=0                                                    
               KOUTR=0                                                    
            ELSE                                                          
               KTOTAL=KRUN+KITS-1                                         
               KOUTP=1-KITS                                               
               KOUTE=1-KITS                                               
               KOUTH=1-KITS                                               
               KOUTR=1-KITS                                               
            ENDIF                                                         
            KSTART=KOUNT                                                  
            DO 183 I=1,IGA                                                
               SPMI(I)=SP(I)                                              
  183       CONTINUE                                                      
            DO 184 J=1,IGB                                                
               ZMI(J)=Z(J)                                                
               DMI(J)=D(J)                                                
               TMI(J)=T(J)                                                
  184       CONTINUE                                                      
            DO 185 KK=1,NTRAC                                             
               DO 185 J=1,IGB                                             
                  TRAMI(J,KK)=TRA(J,KK)                                   
  185       CONTINUE                                                      
         ELSE
            ID=11
  181       IF (NTRACO.EQ.0) THEN
C               READ(ID,END=1001)RKOUNT,RM1TAPE,DAY,Z,D,T,SP,RM2TAPE       
C     +         ,ZMI,DMI,TMI,SPMI,HTNET,RM3TAPE                            
               READ(ID,END=1001)RKOUNT,RM1TAPE,DAY,DOY,Z,D,T,TRA,SP       
     +              ,RM2TAPE,ZMI,DMI,TMI,TRAMI,SPMI,HTNET,RM3TAPE                            
            ELSE
               READ(ID,END=1001)RKOUNT,RM1TAPE,DAY,DOY,Z,D,T              
     +             ,((TRA(I,KK),I=1,IGB),KK=1,NTRACO),SP,RM2TAPE          
     +             ,ZMI,DMI,TMI                                           
     +             ,((TRAMI(I,KK),I=1,IGB),KK=1,NTRACO),SPMI,             
     +             HTNET,RM3TAPE                                          
            ENDIF
            IF (ABS(RM1TAPE-RM2TAPE) .GT. 1.0E-03 .OR.
     +          ABS(RM1TAPE-RM3TAPE) .GT. 1.0E-03) THEN
               WRITE(2,2011) ID                                           
               CALL ABORT                                                 
            ENDIF
            IF (ABS(RM1TAPE-RNTAPO) .GT. 1.0E-03) THEN
               WRITE(2,2021) RM1TAPE,RNTAPO                               
               CALL ABORT                                                 
            ENDIF
C Extra read if restart record contains mass correction.                  
            IF (LMASOLD) THEN                                             
               READ(ID,END=1011) RK,RMT1,RDAY,RREC,GMSP0,GMSPMI,RMT2      
               IF (RMT2.NE.RMT1) THEN                                     
                  WRITE(2,2110)                                           
                  CALL ABORT                                              
               ENDIF                                                      
            ENDIF                                                         
            IF (ABS(DAY-BEGDAY) .GT. 1.0E-02) THEN                        
               IF (ABS(DAY-BEGDAY) .LT. ABS(DAYNEAR-BEGDAY)) THEN         
                  DAYNEAR = DAY                                           
               ENDIF                                                      
               GOTO 181                                                   
            ELSE                                                          
               GOTO 201                                                   
            ENDIF                                                         
 1001       WRITE(2,2031) BEGDAY,DAYNEAR                                  
            CALL ABORT                                                    
 1011       WRITE(2,2130) DAY                                             
            CALL ABORT                                                    
C

  201       KOUNT=NINT(RKOUNT)                                            
            WRITE(2,2041) ID,KOUNT,RM1TAPE,BEGDAY,DOY                     
            KTOTAL=KOUNT+KRUN                                             
            KSTART=KOUNT                                                  
            KTEMP=KOUNT                                                   
            IF(KITS.GT.0) KTEMP=KOUNT+1-KITS                              
            KOUTP=KTEMP-KOUNTP*(KTEMP/KOUNTP)                             
            KOUTE=KTEMP-KOUNTE*(KTEMP/KOUNTE)                             
            KOUTH=KTEMP-KOUNTH*(KTEMP/KOUNTH)                             
            KOUTR=KTEMP-KOUNTR*(KTEMP/KOUNTR)                             
            IF ((KTEMP.GT.0).AND.(KOUTH.EQ.0)) THEN                       
               KOUTH=KOUNTH                                               
            ENDIF                                                         
         END IF
C
         IF (LMASOLD) THEN                                                
            WRITE(2,2100) GMSP0*P0,GMSPMI*P0                              
         ENDIF                                                            
C                                                                         
C        Read in restoration state from separate file                     
C                                                                         
         IF (LRESTIJ) THEN                                                
            ID=13                                                         
  182       READ(ID,END=1002)RKOUNT,RM1TAPE,DAY,DOY,TTRES,RM2TAPE         
            IF (ABS(RM1TAPE-RM2TAPE) .GT. 1.0E-03) THEN
              WRITE(2,2012) ID                                            
              CALL ABORT                                                  
            ENDIF                                                         
            IF (ABS(RM1TAPE-RNTAPO) .GT. 1.0E-03) THEN
               WRITE(2,2022) RM1TAPE,RNTAPO                               
               CALL ABORT                                                 
            ENDIF                                                         
            IF (ABS(DAY-BEGDAY) .GT. 1.0E-02) THEN                        
               IF (ABS(DAY-BEGDAY) .LT. ABS(DAYNEAR-BEGDAY)) THEN         
                 DAYNEAR = DAY                                            
               ENDIF                                                      
               GOTO 182                                                   
            ELSE                                                          
               GOTO 202                                                   
            ENDIF                                                         
 1002       WRITE(2,2032) BEGDAY,DAYNEAR                                  
            CALL ABORT                                                    
C                                                                         
  202       WRITE(2,2042) ID,NINT(RKOUNT),RM1TAPE,BEGDAY,DOY              
         ELSE                                                             
            IF(DAMP.GT.0.0) THEN                                          
               REWIND 13                                                  
               READ(13)ZRES,DRES,TRES,SPRES                               
            END IF                                                        
         ENDIF                                                            
C                                                                         
C  Initialise any new tracer fields and surface state.                    
C                                                                         
C         CALL ICTRAC                                                      
C                                                                         
      END IF                                                              
C                                                                         
C   READS in T21 orography  part                                          
      IF (LOROG) THEN                                                     
c read in spectral field                                                  
        open(50,file='orogdata/t21.50',status='old')                      
        read(50,*) GS                                                     
        close(50)                                                         
C read in gridpoint field                                                 
C skipping extra line on equator                                          
c tstar goes from 1 at pole to jg at equator, gsg made the same           
        open(50,file='orogdata/t21.59',status='old')                      
                                                                          
        read(50,*)((gsg(i,j),i=1,mg),j=1,jg),(adum,i=1,mg),               
     $      ((gsg(i+mgpp,j),i=1,mg),j=jg,1,-1)                            
        close(50)                                                         
c print out land mask                                                     
        write(2,*) 'land mask looks like:'                                
        do j=1,jg                                                         
          write (2,'(128I1)')                                             
     $         (nint((gsg(i,j)/(gsg(i,j)+.1))),i=1,mg)                    
        enddo                                                             
        do j=jg,1,-1                                                      
          write (2,'(128I1)')                                             
     $         (nint((gsg(i+mgpp,j)/(gsg(i+mgpp,j)+.1))),i=1,mg)          
        enddo                                                             
      ENDIF

                                                                          
C If BEGDOY=0 sets DOY from restart record - else                         
C sets it from namelist                                                   
      IF (BEGDOY.GT.0.9) DOY=BEGDOY                                       
      WRITE(2,*) 'STARTING FROM ',DOY,' DAY OF YEAR '                     
      if (lshort) then                                                    
C        print *,'resetting DAY to 0.0'                                    
        DAY=0.0                                                           
      else                                                                
C        print*,' begdoy day doy ',begdoy,day,doy                          
      end if                                                              
      END                                                                 
