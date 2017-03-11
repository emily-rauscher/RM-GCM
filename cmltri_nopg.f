
      PROGRAM IGCM3
      call IGCM3_SUB()
      END


      SUBROUTINE IGCM3_SUB
C!$ use omp_lib
C**********************************************************************   
C                    IGCM3_1                                              
C                                                                         
C FULL PHYSICS + MIXED LAYER OCEAN                                        
C                                                                         
C     DEPARTMENT OF METEOROLOGY      UNIVERSITY OF READING                
C                                                                         
C**********************************************************************   
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
      PARAMETER(IDDAF=NFTGW*IGC,IDDAG=NFTWG*IGC,IDDAD=NFTGD*IGC           
     *         ,IDDZ=(11*NL+4)*JGG)                                       
C                                                                         
C     Constants and arrays needed for balancing                           
C                                                                         
      COMMON/BALAN/BFILT(NL),RGT0(NL),RG(NL2),TMEAN(NL)                   
     +            ,EP1(IGA),EP2(IGA),KBAL,MFTBAL,SRGT0,LTBAL              
      LOGICAL LTBAL                                                       
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

                                                                      
C     Constants and arrays needed for the fast Fourier transforms         
C                                                                         
      COMMON/COMFFT/NTGW,NRSTGW,NTWG,NRSTWG,NTGD,NRSTGD                   
     +              ,TRIG(IDA),WORK(IDF),IFAX(10)                         
C                                                                         
      COMMON/COMGRM/DUDLSG(IGC,NL),DVDLSG(IGC,NL),DTDLSG(IGC,NL)          
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
C                                                                         
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN, TAULIMIT 
      
       LOGICAL LPLOTMAP

                                                                     
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
C                                                                         
C     Polynomial used to aid vectorization of Legendre transforms         
C                                                                         
      COMMON/POLYNO/POLY(NWJ2,2,4),CMPA(IGL)                              
      COMPLEX CMPA                                                        
C                                                                         
      COMMON/PTENDZ/UTVDZ(IGG),VTVDZ(IGG),TTVDZ(IGG),QTVDZ(IGG)           
     :              ,TTCRZ(IGG),QTCRZ(IGG),TTLRZ(IGG),QTLRZ(IGG)          
     :              ,TTRDZ(IGG),CTCRZ(IGG),CTLRZ(IGG)                     
     :              ,UTBLZ(JGG),VTBLZ(JGG),TTBLZ(JGG),QTBLZ(JGG)          
     :              ,AUTVDZ(IGG),AVTVDZ(IGG),ATTVDZ(IGG),AQTVDZ(IGG)      
     :              ,ATTCRZ(IGG),AQTCRZ(IGG),ATTLRZ(IGG),AQTLRZ(IGG)      
     :              ,ATTRDZ(IGG),ACTCRZ(IGG),ACTLRZ(IGG)                  
     :              ,AUTBLZ(JGG),AVTBLZ(JGG),ATTBLZ(JGG),AQTBLZ(JGG)      
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
C     Restoration fields and timescale                                    
C                                                                         
      COMMON/RESTOR/ZRES(IGN),DRES(IGN),TRES(IGN),SPRES(IGM),DAMP         
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
      COMMON/GRIDSS/ZG1(IGD,JG),DG1(IGD,JG),UG1(IGD,JG),VG1(IGD,JG),      
     :              TG1(IGD,JG),SPG1(IGC,JG)                              
C                                                                         
       COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     :               ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     :               ,TTSW(IGC,NL),TTLW(IGC,NL)                           
      COMMON/STATS/GMSP0,GMSPMI,LMASCOR,LMASOLD,LMASPRT                   
      LOGICAL LMASCOR,LMASOLD,LMASPRT                                     
C                                                                         
C                    
      COMMON/NEWFORC/TTRESN(IGB)
      COMPLEX TTRESN
      REAL TMPLAT
C     ER Modif to manage outputs
      INTEGER IDAYS(6),ITSOUT,IFTOUT,ISFOUT
      REAL FDAY
C
C-----------------------------------------------------------------------  
      REAL DAG(IDDAG),DAF(IDDAF),DAD(IDDAD)                               
     *    ,DDZ(IDDZ),ADDZ(IDDZ)                                           
      EQUIVALENCE (DAG(1),UG(1)),(DAF(1),SPG(1)),(DAD(1),TNLG(1))         
     : ,(DDZ(1),UTVDZ(1)),(ADDZ(1),AUTVDZ(1))                             
      INTEGER IODSIZE                                                     
      PARAMETER (IODSIZE=MAX(NGRPAD+IGD,IGP*(3+NTRAC)+IGO)+7)             
      REAL*4 RODATA(IODSIZE)                                              
      REAL R8DATA((3+NTRAC)*IGP)                                          
      REAL R8SP(IGO)                                                      
      EQUIVALENCE (R8DATA(1),Z(1))                                        
      REAL R8TT(IGP)                                                      
      EQUIVALENCE(R8TT(1),TT(1))                                          
      EQUIVALENCE(R8SP(1),SP(1))                                          
       REAL htnet                                                         
       COMMON /RADHT/ HTNET(NHEM,JG,MG,NL)                                
       REAL TAVE(IGP) 
       real POSLATS(JG)

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM, AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS

       CHARACTER(30) :: AEROSOLMODEL
       REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(nl+1)
       LOGICAL DELTASCALE
       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF
     
!       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
!     &   CLOUDTOP,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
!     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,AERO4LAT,AEROPROF,
!     &   TAUAEROSOL
!       CHARACTER(20) :: AEROSOLMODEL       
!       REAL AERO4LAT(NL+1,MG,2),AEROPROF(NL+1)
!       REAL TAUAEROSOL(nl+1,mg,2,jg) 
!       LOGICAL DELTASCALE

!       REAL TAUAEROSOL(nl+1,mg,2,jg)

      NAMELIST/COMMENT/THECOMMENT
      CHARACTER(70) :: THECOMMENT
C                                                                      
 2000 FORMAT(/' RESTART RECORD WRITTEN TO CHANNEL ',I3,/                  
     +        ' RKOUNT  RNTAPE  DAY  DOY  =',4F12.3)                      
 2010 FORMAT(/' HISTORY RECORD WRITTEN TO CHANNEL ',I3,/                  
     +        ' RKOUNT  RNTAPE  DAY  DOY  =',4F12.3)                      
 2020 FORMAT(/' RESTORATION RECORD WRITTEN TO CHANNEL ',I3,/              
     +        ' RKOUNT  RNTAPE  DAY  DOY  =',4F12.3)                      

C!$omp parallel 
C!$     print *,'thread',int(omp_get_thread_num()),'started'
C!$omp do
C      do i=1,16
C      write(*,*) 'Hello',i
C      enddo
C!$omp end do
C!$omp end parallel
C     FIRST THE COMMENT
      READ(7,COMMENT)
      WRITE(2,COMMENT)
      WRITE(*,*) THECOMMENT  


      tripped=0.0
      CALL INISET
      tripped=1.0                                                         
C     if not restart run, then use modified start to set TTRES, otherwise skip
      IF(.NOT.LRSTRT) THEN
C     &&&&&&&&&&&&&MODIFIED START &&&&&&&&&&&&&&
         CALL INITAL         
         DO I=1,IGA          
            VP(I) =0.0       
            SPA(I)=0.0       
         ENDDO               
         DO I=1,IGB          
            ZT(I)=0.0        
            DT(I)=0.0        
            TT(I)=0.0        
            DTE(I)=0.0       
         ENDDO               
 
!      HERE WE MAKE OR READ IN THE AEROSOLS     
!      CALL READ_AEROSOLS
!           IF (AEROSOLS) THEN
!         DO      L  = 1, NL
!          DO     I  = 1, MG
!           DO    JH = 1, JG
!            DO   IHEM = 1,NHEM
!         READ(71,*)  CLAT,CLAY,CHEM,CLON,TAUVAL
C          write(*,*) CLAT,CLAY,CHEM,CLON,TAUVAL
          IF(AEROSOLS) THEN
           DO JH=1,JG
           poslats(JH)=alat(JH)
           ENDDO
          CALL MAKECLOUDS(poslats,p0,sigma)
!          CALL MAKECLOUDS(poslats,mg,jg,nl,p0,sigma)
          ELSE
          PI0AERSW   =  0.0
          ASYMSW     =  0.0
          EXTFACTLW  =  0.0
          PI0AERLW   =  0.0
          ASYMLW     =  0.0
          ENDIF
C         The loop for the radiative transfer code is as follows:
!        DO_LATLOOP (below in cmltr_nopg.f)
!            DO_HEMLOOP (in rradiation.f aka formally cmorc)
!              {START PARALLEL}
!                 DO_LONLOOP (inrradiation.f)
!                      DO_VERTLOOP (actually matrix inversion)
!                          RADIATIVE TRANSFER
!                      ENDDO_VERTLOOP
!                 ENDDO_LONLOOP
!               {END PARALLEL}
!            ENDDO_HEMLOOP
!        ENDDO_LATLOOP
!        So for speed, innermost is leftmost (and should ideally be
!        smallest)  
!        THEREFORE, INDEX ORDER SHOULD BE: TAUAER(NL,MG,IHEM,JH)


C                 
C     Main loop over latitudes ( so says previous author, but it's NOT so for radiation  ~mtr)
C
         JL=1
C
         DO  IH=1,JG        
            JH=IH           
            IF(JGL.EQ.1) THEN
               REWIND(25)
               READ(25) ALP,DALP,RLP,RDLP 
            ENDIF
C                                                
C  Here we define arbitrary temperature profile in grid space.
            IDUM=-1
 
            DO  L=1,NL
               DO  IHEM=1,NHEM  
                  IP=(L-1)*IGC+(IHEM-1)*MGPP 
                  DO  I=1,MG                 
                     IF(IHEM.EQ.1) THEN
                        TMPLAT=PI/180.0*ALAT(IH)
                     ELSE
                        TMPLAT=PI/180.0*ALAT(JG*NHEM-(IH-1))
                     ENDIF

                     TNLG(IP+I)=(TRS(L)-T0(L)) !*RANF(IDUM)

C     For hot Jupiter:
                     IF (COS(1.0*PI2*(I-1)/MG).GE.0) THEN
                        TNLG(IP+I)=TNLG(IP+I) - TRS(L)
     &                       + ( (TRS(L))**4 
     &                       + ( (TRS(L)+REDTEP(L))**4 - (TRS(L))**4 )
     &                       *COS(1*TMPLAT)*COS(1.0*PI2*(I-1)/MG) 
     &                       )**0.25
                     ENDIF
C     (if nightside nothing needed from scheme because TRS=RESTTT=tnight,eq for all lat,long)

                  ENDDO
               ENDDO
            ENDDO
C        Go from grid point space to spectral space using   
C        direct Legendre and Fourier transforms           
C                                                         
            DO I=1,NTGW                
               CALL FFT991(DAF(1+(I-1)*NCRAY*MGPP),WORK,TRIG,IFAX,1,MGPP
     +              ,MG,NCRAY,-1)                                         
            ENDDO                                                            
            
            CALL FFT991(DAF(1+ NTGW*NCRAY*MGPP),WORK,TRIG,IFAX,1,MGPP,MG     
     +           ,NRSTGW,-1)                                           
            
            CALL LTD                                         
            JL=JL+JINC                                       
         ENDDO                                               
C                                                         
         TTRESN=TT

         REWIND(25)
         REWIND(7)
         REWIND(2)
         REWIND(60)
         CALL INISET                                                         
C        &&&&&&&&&&&&&&& END MODIFIED START &&&&&&&&&&&&&&& 
      ENDIF

      kflag=0                                                             
      CALL INITAL                                                         
C     ER modif for output management
      FDAY=KRUN/ITSPD
      IDAYS(1)=5
      IDAYS(2)=INT(FDAY/4.-1.)
      IDAYS(3)=INT(FDAY/2.-1.)
      IDAYS(4)=INT(FDAY*3./4.-1.)
      IDAYS(5)=INT(FDAY-6.)
      IDAYS(6)=INT(FDAY-2.)
      ITSOUT=2600
      IFTOUT=5000
      ISFOUT=6400

!@@@@@ HERE IS WHERE THE TIME STEP ITERATION LOOP BEGINS @@@@@@@@@@@@@@
!@@@@@                                                   @@@@@@@@@@@@@@
!

 1    CONTINUE                                                            
C                                                                         
C     Adiabatic part of timestep. Preset tendencies to zero.              
C                                                                         
      DO 31 I=1,IGA                                                       
         VP(I) =0.0                                                       
         SPA(I)=0.0                                                       
 31   CONTINUE                                                            
      DO 32 I=1,IGB                                                       
         ZT(I)=0.0                                                        
         DT(I)=0.0                                                        
         TT(I)=0.0                                                        
         DTE(I)=0.0                                                       
 32   CONTINUE                                                            
      DO 35 KK=1,NTRAC                                                    
         DO 35 I=1,IGB                                                    
            TRAT(I,KK)=0.0                                                
 35   CONTINUE                                                            
C                                                                         
      IF (KOUNT.EQ.0) THEN                                                
C                                                                         
C        Initialise new tracer fields.                                    
         IF (.NOT. LRSTRT) CALL ICTRAC                                    
C                                                                         
C        Add white noise perturbation                                     
C                                                                         
         IF (LNOISE.AND..NOT.LRSTRT) CALL NOISE                           
C                                                                         
C        Calculate restoration arrays                                     
C                                                                         
         IF (.NOT. LRESTIJ) CALL SETRES                                   
C                                                                         
      ENDIF                                                               
C                                                                         
      IF (JGL.EQ.1) REWIND 25                                             
      KKOUT=KOUNT*(KOUNTP-KOUTP)                                          
      JL=1                                                                
C                                                                         
C set surface pressure to 976mb                                           
C      IF (DAY.LT.0.1) THEN                                                
C        print *,' setting surface pressure to 976mb'                      
C        SP(1)=CMPLX((sqr2*(0.976-1.0)),0.0)                               
C      ENDIF                                                               

!!      print*,kount                                                        
C     Main loop over latitudes                                            
C                                                                         
C!$omp parallel default(shared) private(IH,JH,JL,WORK) 
C!$omp do
      DO 5 IH=1,JG                                                        
         JH=IH                                                            
         IF(JGL.EQ.1) READ(25) ALP,DALP,RLP,RDLP                          
C                                                                         
C        Go from spectral space to grid point space using                 
C        inverse Legendre and Fourier transforms                          
C                                                                         
         CALL LTI                                                         

CC!$omp parallel default(shared) private(I,WORK) 
CC!$omp do
         DO 10 I=1,NTWG                                                   
            CALL FFT991(DAG(1+(I-1)*NCRAY*MGPP),WORK,TRIG,IFAX,1,MGPP,MG  
     +                 ,NCRAY,1)                                          
 10      CONTINUE                                                         
CC!$omp end do 
CC!$omp end parallel

         CALL FFT991(DAG(1+ NTWG*NCRAY*MGPP),WORK,TRIG,IFAX,1,MGPP,MG     
     +              ,NRSTWG,1)                                            
C                                                                         
C        Calculate nonlinear terms                                        
C                                                                         
         CALL MGRMLT                                                      
C                                                                         
C        Save grid point fields for use in XSECT                          
C                                                                         
         DO I=1,IGD                                                       
            ZG1(I,IH)=ZG(I)                                               
            DG1(I,IH)=DG(I)                                               
            UG1(I,IH)=UG(I)                                               
            VG1(I,IH)=VG(I)                                               
            TG1(I,IH)=TG(I)                                               
         END DO                                                           
         DO I=1,IGC                                                       
            SPG1(I,IH)=SPG(I)                                             
         END DO                                                           
C         IF (KKOUT.EQ.0.AND.NLAT.GT.0) WRITE(24)ZG,DG,UG,VG,TG,SPG       
C                                                                         
C        Go from grid point space to spectral space using                 
C        direct Legendre and Fourier transforms                           
C                                                                         

CC!$omp parallel default(shared) private(I,WORK) 
CC!$omp do
         DO 20 I=1,NTGW                                                   
            CALL FFT991(DAF(1+(I-1)*NCRAY*MGPP),WORK,TRIG,IFAX,1,MGPP,MG  
     +                 ,NCRAY,-1)                                         
 20      CONTINUE                                                         
CC!$omp end do 
CC!$omp end parallel

         CALL FFT991(DAF(1+ NTGW*NCRAY*MGPP),WORK,TRIG,IFAX,1,MGPP,MG     
     +              ,NRSTGW,-1)                                           
         CALL LTD                                                         
         JL=JL+JINC   
!!         JL=1+JINC*IH
!!         write(*,*) IH, JH, JL, JINC
 5    CONTINUE                                                            
C!$omp end do 
C!$omp end parallel
C                                                                         
      IF (LBALAN) THEN                                                    
C        Balance spectral fields                                          
C                                                                         
         IF (KOUNT.LT.0) THEN                                             
            KOUNT=KOUNT+1                                                 
            IF (.NOT.LTBAL) CALL BALANC                                   
            IF (     LTBAL) CALL TBAL                                     
            CALL ENERGY                                                   
            GO TO 1                                                       
         ENDIF                                                            
      ENDIF                                                               
C                                                                         
C      IF (KKOUT.EQ.0.AND.NLAT.GT.0) REWIND 24                            
C                                                                         
C     First timestep - output history and diagnostics                     
C                                                                         
      IF (KOUNT.EQ.0) THEN                                                
c         REWIND 9                                                        
         RKOUNT=KOUNT                                                     
         IF (LMINIH) THEN                                                 
           RODATA(1)=RKOUNT                                               
           RODATA(2)=RNTAPE                                               
           RODATA(3)=DAY                                                  
           RODATA(4)=DOY                                                  
           DO J=1,IGP*(3+NTRAC)                                           
             RODATA(J+4)=R8DATA(J)                                        
           ENDDO                                                          
           DO J=1,IGO                                                     
             RODATA(J+4+IGP*(3+NTRAC))=R8SP(J)                            
           ENDDO                                                          
           RODATA(IGO+5+IGP*(3+NTRAC))=RNTAPE                             
           WRITE (9) (RODATA(J),J=1,5+IGO+IGP*(3+NTRAC))                  
C           write (*,*) IODSIZE,5+IGO+IGP*(3+NTRAC)                       
         ELSE                                                             
           WRITE(9) RKOUNT,RNTAPE,DAY,DOY,Z,D,T,TRA,SP,RNTAPE             
         ENDIF                                                            
         WRITE(2,2010)9,RKOUNT,RNTAPE,DAY,DOY                             
         IF (LRESTIJ) THEN                                                
           WRITE(13)RKOUNT,RNTAPE,DAY,DOY,TTRES,RNTAPE                    
           WRITE(2,2020)13,RKOUNT,RNTAPE,DAY,DOY                          
         ENDIF                                                            
         CALL XSECT(INLAT)
         CALL XSECT2
         CALL SPOP                                                        
         CALL ENERGY                                                      
      ELSE IF (KOUNT.EQ.KSTART) THEN                                      
         CALL ENERGY                                                      
      ENDIF                                                               
      IF (LRESTIJ) THEN                                                   
C                                                                         
C        Write a restoration record                                       
C                                                                         
         IF (KOUTH.EQ.KOUNTH.OR.KOUTR.EQ.KOUNTR) THEN                     
            RKOUNT=KOUNT                                                  
            WRITE(13)RKOUNT,RNTAPE,DAY,DOY,TTRES,RNTAPE                   
            WRITE(2,2020)13,RKOUNT,RNTAPE,DAY,DOY                         
         ENDIF                                                            
      ENDIF                                                               
C                                                                         
C     Write a restart record                                              
C                                                                         
      IF (KOUTR.EQ.KOUNTR) THEN                                           
         RKOUNT=KOUNT                                                     
         WRITE(11) RKOUNT,RNTAPE,DAY,DOY,Z,D,T,TRA,SP,RNTAPE              
     +     ,ZMI,DMI,TMI,TRAMI,SPMI,HTNET,RNTAPE                           
         IF (LMASCOR)THEN                                                 
           RREC=2.0                                                       
           WRITE(11) RKOUNT,RNTAPE,DAY,RREC,GMSP0,GMSPMI,RNTAPE           
         ENDIF                                                            
         WRITE(19) -999.999                                               
         WRITE(19) RKOUNT,RNTAPE,DAY,DOY,TSTAR,TDEEP,SMSTAR,QSTAR,        
     $        HSNOW,SQSTAR,SALB,SBAL,TSTARO,TDEEPO,SNET,RNTAPE            
         WRITE(2,2000)11,RKOUNT,RNTAPE,DAY,DOY                            
         KOUTR=0                                                          
      ENDIF                                                               
C                                                                         
C     Write a history record                                              
C                                                                         
      IF (KOUTH.EQ.KOUNTH) THEN                                           
         RKOUNT=KOUNT                                                     
         if (LMINIH) THEN                                                 
            RODATA(1)=RKOUNT                                                
            RODATA(2)=RNTAPE                                                
            RODATA(3)=DAY                                                   
            RODATA(4)=DOY                                                   
C          DO j=1,igb                                                     
C            RODATA(4+j)=z(j)                                             
C          ENDDO                                                          
            DO J=1,IGP*(3+NTRAC)                                            
               RODATA(J+4)=R8DATA(J)                                         
            ENDDO                                                           
            DO J=1,IGP                                                      
               RODATA(J+4+2*IGP)=TAVE(J)/REAL(KOUNTH)                         
               TAVE(J)=0.0                                                    
            ENDDO                                                           
            DO J=1,IGO                                                      
               RODATA(J+4+IGP*(3+NTRAC))=R8SP(J)                             
            ENDDO                                                           
            RODATA(IGO+5+IGP*(3+NTRAC))=RNTAPE                              
            WRITE (9) (RODATA(J),J=1,5+IGO+IGP*(3+NTRAC))                   
         ELSE                                                              
            WRITE(9) RKOUNT,RNTAPE,DAY,DOY,Z,D,T,TRA,SP,RNTAPE              
         ENDIF                                                             
         WRITE(2,2010)9,RKOUNT,RNTAPE,DAY,DOY                              
         KOUTH=0                                                           
      ENDIF                                                               
C                                                                         
C     Output diagnostics                                                  
C                                                                         
      IF (KOUTP.EQ.KOUNTP) THEN                                           
         CALL XSECT(INLAT)
         CALL XSECT2
         IF (INT(DAY).EQ.IDAYS(1)) THEN
            CALL FILECOPY(31,41,81)
         ENDIF
         IF (INT(DAY).EQ.IDAYS(2)) THEN
            CALL FILECOPY(32,42,82)
         ENDIF
         IF (INT(DAY).EQ.IDAYS(3)) THEN
            CALL FILECOPY(33,43,83)
         ENDIF
         IF (INT(DAY).EQ.IDAYS(4)) THEN
            CALL FILECOPY(34,44,84)
         ENDIF
         IF (INT(DAY).EQ.IDAYS(5)) THEN
            CALL FILECOPY(35,45,85)
         ENDIF
         IF (INT(DAY).EQ.IDAYS(6)) THEN
            CALL FILECOPY(36,46,86)
         ENDIF
C        ER Modif for KE spectrum output
         CALL XSECT3
CC      call plotfields(0)                                                
         kflag=1                                                             
         CALL SPOP                                                        
         KOUTP=0                                                          
      ENDIF                                                               
C     ER modif for 90 outputs during last orbit of planet
      IF (KOUNT.GT.(KTOTAL-ABS(PORB)*ITSPD)) THEN       ! Note: if porb=0, won't do this
         CALL FINALORB(KOUNT,KTOTAL,ABS(PORB*ITSPD),ITSOUT,IFTOUT,ISFOUT)
      ENDIF
      IF (KOUTE.EQ.KOUNTE) THEN                                           
         CALL ENERGY                                                      
         KOUTE=0                                                          
      ENDIF                                                               
C                                                                         
      IF (KOUNT.LT.KTOTAL) THEN                                           
         KOUTP=KOUTP+1                                                    
         KOUTE=KOUTE+1                                                    
         KOUTH=KOUTH+1                                                    
         KOUTR=KOUTR+1                                                    
         KOUNT=KOUNT+1                                                    
         IF(KOUNT.EQ.1.AND.KITS.GT.0) DAY=DAY+DELT/PI2                    
         DAY=DAY+DELT/PI2                                                 
         IF (.NOT. LPERPET) THEN                                          
           IF(KOUNT.EQ.1.AND.KITS.GT.0) DOY=DOY+DELT/PI2                  
           DOY=DOY+DELT/PI2                                               
           IF (DOY.GE.361) DOY = DOY-360.0                                
         ENDIF                                                            
C                                                                         
C Adiabatic part of timestep                                              
C                                                                         
         CALL TSTEP                                                       
C                                                                         
C Mass correction.                                                        
C                                                                         
         IF (LMASCOR) THEN                                                
            CALL MASCOR                                                   
         ELSE IF(LRESTIJ) THEN                                            
C Simple fix to maintain mass.                                            
            SP(1)=CMPLX(0.,0.)                                            
         ENDIF                                                            
C                                                                         
C Diabatic part of timestep. Preset tendencies to zero.                   
C                                                                         
         DO 33 I=1,IGB                                                    
            ZT(I)=0.0                                                     
            DT(I)=0.0                                                     
            TT(I)=0.0                                                     
   33    CONTINUE                                                         
         DO 36 KK=1,NTRAC                                                 
            DO 36 I=1,IGB                                                 
               TRAT(I,KK)=0.0                                             
   36    CONTINUE                                                         
C                                                                         
C        Set parameters for optional Newtonian cooling.                   
C                                                                         
         IF (LRESTIJ) CALL SETTEE                                         
C                                                                         
C        Preset accumulated diagnostics.                                  
C                                                                         
         IF (KOUTH.EQ.1) THEN                                             
            DO 232 I=1,IDDZ                                               
               ADDZ(I)=0.0                                                
  232       CONTINUE                                                      
         ENDIF                                                            
C                                                                         
C        Loop over latitude for spectral transforms                       
C        and calculation of diabatic tendencies.                          
C                                                                         
         IF (JGL.EQ.1) REWIND(25)                                         
C      REWIND NAVRD                                                       
C      REWIND NAVWT 

                                                      
!@@@@@@  HERE IS WHERE THE ITERATION OVER LATITUDE FOR  @@@@@@@@@@@@@@@@
!@@@@@@  THE RADIATIVE TRANSFER BEGINS. SHOULD BE PARALELLIZED @@@@@@@@@

         JL=1                                                             
         DO 260 IH=1,JG                                                   
            JH=IH                                                         
            IF(JGL.EQ.1) READ(25) ALP,DALP,RLP,RDLP                       
 
C     FOR RADIATIVE TRANSFER
C       IF the aerosol keyword is set true, then here's where we extract
C       a row from the aerosol array, representing vertical profiles of
C       aersols for each longitude for the given latitude being currently indexed.
C       aer4thislat=aerosolarry(alllays,allons,IH)

C         The loop for the radiative transfer code is as follows:
!        DO_LATLOOP (here in cmltr_nopg.f)
!            DO_HEMLOOP (in rradiation.f aka formally cmorc)
!              {START PARALLEL}
!                 DO_LONLOOP (inrradiation.f)
!                      DO_VERTLOOP (actually matrix inversion)
!                          RADIATIVE TRANSFER
!                      ENDDO_VERTLOOP
!                 ENDDO_LONLOOP
!               {END PARALLEL}
!            ENDDO_HEMLOOP
!        ENDDO_LATLOOP
!
!       IF(AEROSOLS) THEN
!       AERO4LAT=TAUAEROSOL(:,:,:,IH) 
!       ENDIF
!         WRITE(*,*) 'TEST',TAUAEROSOL(:,1,1,16)
!         WRITE(*,*) 'TEST',TAUAEROSOL(:,1,10,16)
!         WRITE(*,*) 'TEST',TAUAEROSOL(:,1,20,16)
!         WRITE(*,*) 'TEST',TAUAEROSOL(:,1,30,16)
!         WRITE(*,*) 'TEST',TAUAEROSOL(:,1,40,16)
!         WRITE(*,*) 'TEST',TAUAEROSOL(:,1,50,16)
!         WRITE(*,*) 'TEST',TAUAEROSOL(:,1,60,16)
!         write(*,*) tauaerosol
C                                                                         
C        Go from spectral space to grid point space using                 
C        inverse Legendre and Fourier transforms                          
C                                                                         
            CALL LTI                                                      

CC!$omp parallel default(shared) private(I,WORK) 
CC!$omp do
            DO 210 I=1,NTWG                                               
               CALL FFT991(DAG(1+(I-1)*NCRAY*MGPP),WORK,TRIG,IFAX,        
     +                 1,MGPP,MG,NCRAY,1)                                 
 210        CONTINUE                                                      
CC!$omp end do 
CC!$omp end parallel

            CALL FFT991(DAG(1+ NTWG*NCRAY*MGPP),WORK,TRIG,IFAX,           
     +                  1,MGPP,MG,NRSTWG,1)                               
C                                                                         
C        Calculate diabatic terms                                         
C                                                                         
            CALL DGRMLT(IH)                                                   
C                                                                         
C        Write accumulated diagnostics to history file.                   
            if (kflag.eq.1.and.nlat.gt.0) write(24)grpad                        
            RKOUNT=KOUNT                                                  
            RIH=IH                                                        
            IF (KOUTH.EQ.KOUNTH) THEN                                     
               RODATA(1)=RKOUNT                                            
               RODATA(2)=RNTAPE                                            
               RODATA(3)=DAY                                               
               RODATA(4)=DOY                                               
               RODATA(5)=RIH                                               
               IF (LSHIST) THEN                                            
                  DO J=1,NGRPAD/2                                           
                     RODATA(J+5)=GRPAD(J)                                    
                  ENDDO                                                     
                  RODATA(NGRPAD/2+6)=RNTAPE                                 
                  DO J=1,IGD                                                
                     RODATA(NGRPAD/2+6+J)=TNLG(J)                            
                  ENDDO                                                     
                  RODATA(NGRPAD/2+7+IGD)=RNTAPE                             
                  WRITE (9) (RODATA(J),J=1,NGRPAD/2+IGD+7)                  
               ELSE                                                        
                  DO J=1,NGRPAD                                             
                     RODATA(J+5)=GRPAD(J)                                    
                  ENDDO                                                     
                  RODATA(NGRPAD+6)=RNTAPE                                   
                  DO J=1,IGD                                                
                     RODATA(NGRPAD+6+J)=TNLG(J)                              
C  TNLG (gridded heating) is not required in history file. When new       
C  version of flux programme is created TNLG should be removed            
C  from RODATA.                                                           
                  ENDDO                                                     
                  RODATA(NGRPAD+7+IGD)=RNTAPE                               
                  WRITE (9) (RODATA(J),J=1,NGRPAD+IGD+7)                    
               ENDIF                                                       
            ENDIF                                                         
C                                                                         
C        Go from grid point space to spectral space using                 
C        direct Legendre and Fourier transforms                           
C                                                                         
CC!$omp parallel default(shared) private(I,WORK) 
CC!$omp do
            DO 220 I=1,NTGD                                               
               CALL FFT991(DAD(1+(I-1)*NCRAY*MGPP),WORK,TRIG,IFAX,        
     +                     1,MGPP,MG,NCRAY,-1)                            
 220        CONTINUE                                                      
CC!$omp end do 
CC!$omp end parallel
            CALL FFT991(DAD(1+ NTGD*NCRAY*MGPP),WORK,TRIG,IFAX,           
     +                  1,MGPP,MG,NRSTGD,-1)                              
C                                                                         
            CALL LTDDIA                                                   
            JL=JL+JINC                                                    
 260     CONTINUE                                                         
C                                                                         
C Write zonally averaged diagnostics and spectral heating                 
C to history file.                                                        
C                                                                         
      if (kflag.eq.1) then                                                
         if (nlat.gt.0) rewind 24                                         
CC         call plotfields2(1)                                            
CC      call netout2                                                      
         kflag=0                                                          
      end if                                                              
         IF(KOUTH.GE.1.AND.KOUTH.LE.KOUNTH) THEN                          
           DO 270 I=1,IDDZ                                                
  270      ADDZ(I)=ADDZ(I)+DDZ(I)                                         
         ENDIF                                                            
         IF(KOUTH.EQ.KOUNTH) THEN                                         
           RKP=1.0/REAL(KOUNTH)                                           
           DO 280 I=1,IDDZ                                                
  280      ADDZ(I)=ADDZ(I)*RKP                                            
           RKOUNT=KOUNT                                                   
           RIH=JGP                                                        
           IF (LMINIH) THEN                                               
             RODATA(1)=RKOUNT                                             
             RODATA(2)=RNTAPE                                             
             RODATA(3)=DAY                                                
             RODATA(4)=DOY                                                
             RODATA(5)=RIH                                                
             DO j=1,iddz                                                  
               RODATA(J+5)=ADDZ(J)                                        
               RODATA(J+IDDZ+5)=DDZ(J)                                    
             ENDDO                                                        
             RODATA(6+IDDZ*2)=RNTAPE                                      
             WRITE (9) (RODATA(J),J=1,IDDZ*2+6)                           
             RIH=RIH+1                                                    
             RODATA(5)=RIH                                                
             DO J=1,IGP                                                   
               RODATA(J+5)=R8TT(J)                                        
             ENDDO                                                        
             RODATA(IGP+6)=RNTAPE                                         
             WRITE (9) (RODATA(J),J=1,IGP+6)                              
           ELSE                                                           
             WRITE(9)RKOUNT,RNTAPE,DAY,DOY,RIH,ADDZ,DDZ,RNTAPE            
             RIH=RIH+1                                                    
             WRITE(9)RKOUNT,RNTAPE,DAY,DOY,RIH,TT,RNTAPE                  
           ENDIF                                                          
         ENDIF                                                            
C      NTEMP=NAVRD                                                        
C      NAVRD=NAVWT                                                        
C      NAVWT=NTEMP                                                        
 300     CONTINUE                                                         
C                                                                         
C Apply dissipation and optional linear restoration.                      
C                                                                         
         CALL DIFUSE                                                      
C                                                                         
C Update spectral fields in the diabatic timestep.                        
C                                                                         
         CALL DSTEP                                                       
C                                                                         
C End of timestep                                                         
         DO J=1,IGB                                                       
           TAVE(j*2-1)=TAVE(J*2-1)+REAL(T(J))                             
           TAVE(j*2)=TAVE(J*2)+AIMAG(T(J))                                
         ENDDO                                                            
C                                                                         
         GO TO 1 ! time loop                                              
      ENDIF ! end of IF(KOUNT.LT.KTOTAL)                                  
C                                                                         
C Write the final restart record                                          
C                                                                         
CC      call end_netcdf                                                   
      RKOUNT=KOUNT                                                        
      WRITE(12) RKOUNT,RNTAPE,DAY,DOY,Z,D,T,TRA,SP,RNTAPE                 
     +     ,ZMI,DMI,TMI,TRAMI,SPMI,HTNET,RNTAPE                           
      IF(LMASCOR)THEN                                                     
         RREC=2.0                                                         
         WRITE(12) RKOUNT,RNTAPE,DAY,RREC,GMSP0,GMSPMI,RNTAPE             
      ENDIF                                                               
      WRITE(17) -999.999                                                  
      WRITE(17) RKOUNT,RNTAPE,DAY,DOY,TSTAR,TDEEP,SMSTAR,QSTAR,           
     $        HSNOW,SQSTAR,SALB,SBAL,TSTARO,TDEEPO,SNET,RNTAPE            
      WRITE(2,2000)12,RKOUNT,RNTAPE,DAY,DOY                               
      IF (LRESTIJ) THEN                                                   
         WRITE(13)RKOUNT,RNTAPE,DAY,DOY,TTRES,RNTAPE                      
         WRITE(2,2020)13,RKOUNT,RNTAPE,DAY,DOY                            
      ENDIF                                                               
CC      call end_graphics                                                 
C                          

!      IF(LPLOTMAP) THEN
C Close the graphics device.
c      CALL PGASK
C      CALL PGCLOS
!      ENDIF
                                               
      STOP                                                                
      END SUBROUTINE
