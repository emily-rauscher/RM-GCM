C**********************************************************               
C             SUBROUTINE RADIATION                                             
C**********************************************************               
      SUBROUTINE RADIATION(TROPHT,IH)                                             
C                                                                         
C     RADIATION SCHEME DERIVED FROM PREVIOUS CMORC.F AND THE
C     TOON CODES (TOON ET AL 1989). THE SCHEME IS CURRENTLY DOUBLE GRAY
C     AND APPLIES THE TWO-STREAM APPROXIMATION WITH QUADRATURE IN THE
C     VISIBLE AND HEMISPHERIC MEAN IN THE INFRARED. 
                          
C     It passes the pressure of the full sigma levels and the surface                                   
C     to the Radiation scheme temperatures from TG and TSTAR                                  
C                                                                         
C     Determines model resolution                                         
      use omp_lib                                                              
      include 'params.i'
      
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
C     Array ordering in GRIDP must correspond to that in SPECTR.          
C     Real arrays: multi-level arrays are 2-dimensional.                  
C     the variables have been renamed to coincide with                    
C     variable names in bgcm5 DGRMLT                                      
C                                                                         
C swapped round VLNG and UNLG to check                                    
      COMMON/GRIDP/ CHIG(IGC,NL),SFG(IGC,NL),UG(IGC,NL),VG(IGC,NL)        
     :              ,TTVD(IGC,NL),QTVD(IGC,NL),TG(IGC,NL)                 
     :              ,TRAG(IGC,NL,NTRAC)                                   
     :              ,PLG(IGC),TYBL(IGC),TXBL(IGC)                         
     :              ,SPG(IGC),VPG(IGC),TTRD(IGC,NL)                       
     :              ,TNLG(IGC,NL),TRANLG(IGC,NL,NTRAC),UNLG(IGC,NL)       
     :              ,VNLG(IGC,NL),TTLR(IGC,NL),UTRAG(IGC,NL,NTRAC)        
     :              ,TTCR(IGC,NL),VTRAG(IGC,NL,NTRAC)                     
     :              ,UTVD(IGC,NL),VTVD(IGC,NL)                            
     :         ,ASSBL(IGC),ASHBL(IGC),ASLBL(IGC),ARRCR(IGC),ARRLR(IGC)    
     :         ,arflux(igc,6),asfld(igc,6),acld(igc,4)                    
     :         ,SSBL(IGC),SHBL(IGC),SLBL(IGC),RRCR(IGC),RRLR(IGC)         
     :         ,rflux(igc,6),sfld(igc,6),cld(igc,4)                       
C                                                                         
C                                                                         
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 
      
       LOGICAL LPLOTMAP

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

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM, AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS

       CHARACTER(30) :: AEROSOLMODEL
       REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1)
       LOGICAL DELTASCALE
       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF

!       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
!     &   CLOUDTOP,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
!     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF
!       CHARACTER(20) :: AEROSOLMODEL
!       REAL TAUAEROSOL(nl+1,mg,2,jg)
!       REAL AERORPROF
!       LOGICAL DELTASCALE

!       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
!     &   CLOUDTOP,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
!     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,AERO4LAT,AEROPROF,
!     &   TAUAEROSOL
!       CHARACTER(20) :: AEROSOLMODEL
!       REAL AERO4LAT(NL+1,MG,2),AEROPROF(NL+1)
!       REAL TAUAEROSOL(nl+1,mg,2,jg)
!       LOGICAL DELTASCALE


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
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
      COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     :     ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     :     ,TTSW(IGC,NL),TTLW(IGC,NL)                           
      COMMON/GSG/GSG(IGC,JG)                                             
      REAL htnet                                                         
      COMMON /RADHT/ HTNET(NHEM,JG,MG,NL)                                
      REAL TAVE(IGP) 
                                                                          
       REAL PR(NL+1),T(NL+1),PRFLUX(nl+1),htlw(nl+1),htsw(nl+1)
       real PRB2T(NL+1),adum                                            


      integer ifirst                ! If =1, first time reading o3        
                                    ! and h2o (2 months' worth).          
                                                                          
      real amfrac                   ! fraction through month              
      integer ichange               ! =1 when in process of month change  
*---------------------------                                              
      integer ifirstcol             ! =1 first time through column        
                                    ! calculation (open new file).        
      real p0                                                             
      real ps                       ! sfc pressure (used in               
                                    ! interpolation from climatology      
                                    ! to model).                          
      integer im                    ! Pointer for array plg (for          
                                    ! getting sfc pressure).              
                                                                          
C     Array to hold fluxes at top and bottom of atmosphere                
C     1st index - flux 1=SW, 2=LW                                         
C     2nd index - Direction 1=DN, 2=UP                                    
C     3rd index - Where 1=TOP, 2=SURFACE                                  
      real fluxes(2,2,2) 

c     The following for parallel testing --MTR
      integer TID, NTHREADS             
      double precision test_wctime                              
                                                                        
      save                          ! Want to keep things like dcompl.    
                                                                          
                                                                          
                                                                          
      DATA IFIRST/1/                                                      
      data ifirstcol/1/                                                   
                                                                          
      RHSCL=288.0*GASCON/GA                                               
      CHRF=86400.*WW*CT   ! factor to non-dimensionalise heating rates     
c  Skipping part set up here.                                             
c  Radiation scheme only called every nskip longitudes                    
c  nskip must divide exactly into mg for longitude.                       
c  ntstep is the number of timesteps to skip.                             
                                                                          
C      ntstep=max(2,itspd/100000)       ! N.B. No point any more frequent as      
      ntstep=NTSTEP_IN
                         ! Morcrette doing diurnal averages anyway.       

!      write(*,*) 'ntstep',ntstep
      nskip=NSKIP_IN                                                       
C this sets the non-dimensional gridpoint temperature tendency (TTRD)     
C to get to TTRD from K/day divide the K/day heating rate by              
C     (86400*WW*CT)                                                       
C     TTRD(non-dim)=HTRT(k/day)/(86400*WW*CT)                             
c
   
                                           
 
c
c --------------------------------- Now start the radiation bit.          
c                                                                         
C loop over hemispheres                                                   
      IOFM=0                                                              
C
                                                                         
      DO 800 ihem=1,nhem                                                  
                                                                          
c calculates heating rates every ntstep time steps                        
         IF (mod(kount,ntstep).eq.1) then                                   

!          write(*,*) 'Starting Radiation scheme'
C                                                                         
C  Does do Radn scheme                                                    
C                                                                         
C loop over longitudes for radn calculation
!!!  NOTE, THIS PARALLELISM IS ONLY DOING THE LONGITUDES IN PARALLEL,
!    NOT THE LAITITUDES. THE |LATITUDE| LOOP STARTS IN DGRMLT, AND
!    HEMISPHERE LOOP (TO DO BOTH + AND - LAT VALUE) BEGINS ABOVE.  IF WE BEGIN
!    THE PARALLELISM OUTSIDE THE LATITUDE LOOP, THEN WE WOULD COMPUTE EACH
!    |LATITUDE| IN PARALLEL, BUT THE LONGITUDES AND HEMISPHERES WOULD BE DONE SERIALLY
!    WITHIN EACH THREAD, UNLESS NESTING IS POSSIBLE. 
!    THE NUMBER OF LONGITUDES EXCEEDS THE NUMBER OF
!    LATITUDES (MORE SO WHEN CONSIDERING HEMISPHERE SEPERATELY) SO THE
!    BENEFIT IS POTENTIALLY GREATER IF THE NUMBER OF CORES EXCEEDS THE
!    NUMBER OF LATITUDES. NESTING SHOULD SPEED THIS UP FURTHER. 
!      THE ARCHITECTURE IS AS FOLLOWS
!        DO_LATLOOP (IN DGRMLT)
!            DO_HEMLOOP
!              {START PARALLEL}
!                 DO_LONLOOP
!                      DO_VERTLOOP
!                          RADIATIVE TRANSFER
!                      ENDDO_VERTLOOP
!                 ENDDO_LONLOOP
!               {END PARALLEL}
!            ENDDO_HEMLOOP
!        ENDDO_LATLOOP
!
!            write(*,*) MAXVAL(AERO4LAT)
            ilast=0                                                             
c             write(*,*) 'line431'
c            TID = OMP_GET_MAX_THREADS()
c            write(*,*) 'MAXTHREADS:',TID, ' MG=',mg
c IMPLEMENTING PARALLEL PROCESSING! ENDS AT LN 642 --MTR

!  NOTE THAT THIS NEEDS TO BE UPDATED FOR THE NEW RADIATIVE TRANSFER
!  SCHEME!!!!!!!WARNING!
!$OMP PARALLEL DO schedule(guided), default(none), private(test_wctime),
!$OMP& private(im,idocalc,imp,PR,T,imm,alat1,cf,ic,SWALB,alon,htlw,
!$OMP& htsw,HTNETO,a,b),
!$OMP& shared(iofm,nskip,AMFRAC,h2omod2,ihem,h2omod1,o3mod2,o3mod1,TROPHT,
!$OMP& QG,SSLAT,SSLON,fluxes,CHRF),
!$OMP& firstprivate(ilast),
!$OMP& lastprivate(ilast),
!$OMP& shared(ABSLW, ABSSW, ACLD, AIOCT, AK,
!$OMP& AKAP, AKQC, AKQV, AKTC, AKTV, AKVV, ALAT, ALBSW1, ALP, ALPHA,
!$OMP& ALPJ, AQ, ARFLUX, ARRCR, ARRLR, ASFLD, ASHBL, ASLBL, ASSBL,
!$OMP& ASYMSW2, AW, BEGDAY, BEGDOY, BLA, BLCD, BLRH, BLVAD, BLVB, BM1,
!$OMP& BOTRELAXTIME, C, CBADJP, CBADJT, CCC, CCR, CD, CFRAC, CG, CHIG,
!$OMP& CLATNT, CLD, CLR, CPD, CQ, CS, CSSQ, CT, CTCR, CTLR, CTQ, CTQI,
!$OMP& CTRA, CUBMT, CURHM, CUT1, CUT2, CV, DALP, DALPJ, DAY, DELT,
!$OMP& DELT2, DELT2C, DOY, DRAG, DSIGMA, DTBUOY, EAM1, EAM2, ECCEN,
!$OMP& EPSIQ, ESCONA, ESCONB, EZ, FB, FBASEFLUX, FORCE1DDAYS, FRAD, FWS,
!$OMP& G, GA, GASCON, GSG, GWT, HSNOW, HTNET, ICFLAG, INLAT, INSPC,
!$OMP& ITSLL, ITSLO, ITSPD, JH, JINC, JL, JSKIPLAT, JSKIPLON, JZF, KITS,
!$OMP& KOLOUR, KOUNT, KOUNTE, KOUNTH, KOUNTP, KOUNTR, KOUTE, KOUTH,
!$OMP& KOUTP, KOUTR, KRUN, KSTART, KTOTAL, L1DZENITH, L22L, LBALAN, LBL,
!$OMP& LCBADJ, LCLIM, LCOND, LCR, LCSFCT, LCUBM, LDIUR, LFLUX,
!$OMP& LFLUXDIAG, LGPO, LLOGPLEV, LLR, LMINIH, LNNSK, LNOICE, LNOISE,
!$OMP& LOC, LOGICALLBL, LOGICALLLOGPLEV, LOLDBL, LOROG, LPERPET,
!$OMP& LPLOTMAP, LRD, LRESTIJ, LRSTRT, LSHIST, LSHORT, LSL, LSPO,
!$OMP& LSTRETCH, LTVEC, LVD, MF, MFP, NAVRD, NAVWT, NCOEFF, NCUTOP,
!$OMP& NEWTB, NEWTE, NF, NFP, NLAT, NLCR, NLPLOTMAP_IN, NLWMODEL,
!$OMP& NSKIP_IN, NSWMODEL, NTRACO, NTSTEP_IN, OBLIQ, OOM_IN,
!$OMP& OPACIR_POWERLAW, OPACIR_REFPRES, P0, PFAC, PLG, PNET, PNU, PNU2,
!$OMP& PNU21, PORB, PRFLUX, PRMIN, QC, QSTAR, QTDC, QTMC, QTVD, RADEA, RCON, RD,
!$OMP& RDLP, RDSIG, RFCOEFF_IN, RFLUX, RGG, RLP, RMG, RNTAPE, RNTAPO,
!$OMP& RRCR, RRFLUX, RRLR, RSQ, RSQR2, RV, SAICE, SALB, SASNOW, SBAL,
!$OMP& SCATSW2, SD1, SD2, SDSN, SDSND, SDW, SECSQ, SFG, SFLD, SHBL,
!$OMP& SHCI, SHCO, SHCS, SHCSN, SHCSP, SHSMAX, SHSSTAR, SI, SIGMA,
!$OMP& SIGMAH, SISQ, SK, SKAP, SKSE, SKSN, SLBL, SLHF, SMSTAR, SNET,
!$OMP& SOLC_IN, SPG, SQ, SQH, SQR2, SQSTAR, SSBL, SSMC, SVEGE, T0,
!$OMP& T01S2, TAU, TC, TDEEP, TDEEPO, TG, TKP, TNLG, TOAALB, TOUT1,
!$OMP& TOUT2, TRAG, TRANLG, TSLA, TSLB, TSLC, TSLD, TSTAR, TSTARO, TTCR,
!$OMP& TTDC, TTLR, TTLW, TTMC, TTRD, TTSW, TTVD, TXBL, TYBL, UG, UNLG,
!$OMP& UTRAG, UTVD, VG, VNLG, VPG, VTRAG, VTVD, WW)
            DO i=1,mg                                                         
         !      write(*,*) 'i (cmorc ln 433)',i  !MTR MODIF
c         write(*,*),mg
cM         write(*,*), 'Thread id=', omp_get_thread_num(),i
cm               test_wctime=omp_get_wtime()
               im=i+iofm   
               idocalc=0                                                        
               IF ((i.eq.1).or.(i-ilast.ge.nskip)) then                         
                  idocalc=1                                                       
               ELSE                                                             
                  IF (LNNSK) THEN
c                     write(*,*) 'line444 stop'
Cm                     stop                                                 
                     imp=im+1                                                       
                     IF (imp.gt.(mg+iofm)) imp=1+iofm                               
                     imm=im-1                                                       
                    IF (((gsg(im,jh).gt.0.).and.(gsg(imp,jh).eq.0.)).or.           
     $                  ((gsg(im,jh).eq.0.).and.(gsg(imp,jh).gt.0.)).or.              
     $                  ((gsg(im,jh).gt.0.).and.(gsg(imm,jh).eq.0.)).or.              
     $                ((gsg(im,jh).eq.0.).and.(gsg(imm,jh).gt.0.))) THEN            
                        idocalc=1                                                     
                     ENDIF                                                          
                  ENDIF                                                           
               ENDIF                                                            
               IF (idocalc.eq.1) then                                           
                                                                          
c ------------------------------ First set SURFACE VALUES                 
c Pressure, units of Pa 

!MTR                  write(*,*) 'PLG in rradiation', PLG
!MTR                  write(*,*) 'P0', P0
                                                  
!MTR                  PR(1)=PLG(im)*P0   

!MTR                  write(*,*) 'PR(1)=',PR(1) 
!MTR                  write(*,*) 'SIGMA',SIGMA
!MTR                  write(*,*) 'CT', CT                                            
c T, units of K                                                           
!          T(1)=TSTAR(IM,JH)*CT                                            
c          T(1)=288.00 !+((rrflux(IM,JH,1)+rrflux(IM,JH,3)-rrflux(IM,JH,2)-rrflux(IM,JH,4))/5.6704e-8)**0.25

C          IF(SNET(IM,JH).GE.0) THEN
C             T(1)=0.0 +(SNET(IM,JH)/5.6704e-8)**0.25

C      T(1)=((FBASEFLUX+(rrflux(IM,JH,1)+rrflux(IM,JH,3))*3.1416)
C     &         /5.6704e-8)**0.25    !PI adjustment

C      T(1)=((FBASEFLUX+rrflux(IM,JH,1)+rrflux(IM,JH,3))/5.6704e-8)**0.25
C     ER Modif, to be consistent with cnikos.f (no downward flux if f-diff)

!MTR                  T(1)=((FBASEFLUX+rrflux(IM,JH,1))/5.6704e-8)**0.25

c --------------------------------------- Now set rest of column.         
                  DO LD=1,NL    ! Start of loop over column.     
                     L=NL-LD+2  ! Reverse index (Morc goes bottom up).                     
                     PR(LD)=SIGMA(LD)*PLG(im)*P0 ! Pressure 
                     PRB2T(L)=PR(LD)                      
                     T(LD)=TG(im,ld)*CT ! Temperature 
                     AEROPROF(LD)=0.0 
                  ENDDO
                     AEROPROF(NL+1)=0.0
                   PRB2T(1)=PLG(im)*P0
                   PR(NL+1)=PLG(im)*P0
                   T(NL+1)=((FBASEFLUX+rrflux(IM,JH,1))/5.6704e-8)**0.25
! This could also just be the ground temperature... decision to be made
c ----------------------------------------------------- And alat1         
                                                                          
                  alat1=alat(JH)*REAL(-(ihem*2.)+3)
                  IF (KOUTP.EQ.KOUNTP-1) THEN
                     IF(JH.EQ.1.AND.IHEM.EQ.1.AND.I.EQ.1) THEN
                        REWIND(63) !! Rewind file for fluxes in nikosrad
                        REWIND(62) ! rwnd file for ancillary RT results 
                        IF (PORB.NE.0) THEN 
                           SSLON=(1./PORB-1.)*KOUNT*360./ITSPD
                           SSLON=MOD(SSLON,360.)
                        ELSE
                           SSLON=0.
                        ENDIF
                        SSLAT=ASIN(SIN(OBLIQ*PI/180.)
     +                       *SIN(PI2*KOUNT/ITSPD/PORB))*180./PI
                        WRITE(63,2021) DAY,SSLON,SSLAT
                        WRITE(62,2021) DAY,SSLON,SSLAT
 2021                  FORMAT('DAY:',F7.2,', SUBSTELLAR LON,LAT:',2F7.2)
                        WRITE(63,*)
                        WRITE(62,*)''
!       WRITE(62,*)'Pressure (bars) of SW clear gas tau = 2/3 @ mu0=1',
!     &               (2./3.)/(ABSSW*1e6/GA/100.)
!
!       WRITE(62,*)'Pressure (bars) of LW clear gas tau = 2/3 @ mu0=1',
!     &               (2./3.)/(ABSLW*1e6/GA/100.)
!       write(62,*)''
                     ENDIF
                  ENDIF
                                                                          
!! BOTTOM SHORT WAVE ALBEDO SET IN FORT.7 INISIMPRAD,
! Note to future modelers-this is not location or wavelenght dependent
             SWALB=ALBSW   
                  
! PR AND T ARE THE TEMPERATURE AND PRESSURE AT THE SIGMA LEVELS
! AND BOTTOM BOUNDARY, AS USED BY THE DYNAMICAL CODE. 
! TO COMPUTE HEATING RATES AT THESE CENTERS, WE NEED TO DEFINE
! LAYER EDGES AT WHICH FLUXES ARE COMPUTED, PRFLUX.
             DO LD    = 1,NL-1      
             PRFLUX(LD+1)=(pr(LD)+pr(LD+1))/2. 
             ENDDO
             PRFLUX(NL+1)=PR(NL+1)
!             write(*,*),'PR',PR
             PRFLUX(1)=pr(1)*0.5
!             write(*,*)'fluxes',fluxes
C cloud cf and ic passed. fluxes returned.                                
C which is net flux at TOA in profile                                     
C Call radiation scheme                                                   
                  alon=REAL(i-1)/REAL(mg)*360.0  

!     PR in pascals for layer boundaries (NL+1), T in Kelvin for layer
!     centers + one layer for the bottom boundary. The top is n=1, the
!     bottom is n=NL+1

            IF(AEROSOLS) THEN
!           AEROSOLS TIME!
!           Extract a single column from the array AER4LAT(NLEV,LON,HEM)
!            AEROPROF=AERO4LAT(:,mg,ihem)
              DO  LD=1,NL +1
              AEROPROF(LD)=TAUAEROSOL(LD,i,ihem,ih)
             ENDDO
!            write(*,*) alat1,alon
!            write(*,*) aeroprof
!            write(*,*) ''
!            WRITE(*,*) ''
!            write(*,*) AERO4LAT(:,10,1)
!            WRITE(*,*) ''
!            write(*,*) AERO4LAT(:,10,2)
            
!            write(*,*) 'TAUAEROSOL IN RRAD',TAUAEROSOL 
            ENDIF
            call calc_radheat(pr,t,prflux,alat1,alon,htlw,htsw,DOY,cf,           
     $                 ic,fluxes,swalb,kount,itspd)
                                 
c                  call nikosrad(pr,t,h2o,o3,alat1,htlw,htsw,DOY,cf,ic,            
c     $                 fluxes,swalb,alon,kount,itspd)
          
c          write(*,*) "Just called Nikosrad"
!           write(*,*) 'in radiation...'
!           write(*,*) 'htlw', htlw
!           write(*,*) 'htsw', htsw
c       flip indexing of pr to be consistent with bottom to top convention
           pr=prb2t   
!           write(*,*) 'pr',pr
!           write(*,*) 'fluxes',fluxes          
!           write(*,*) 'stop following calc_radheat'
          
                                                                
c store net flux in PNET                                                  
                  PNET(IM,JH)=fluxes(1,1,1)-fluxes(1,2,1)+fluxes(2,1,1)-          
     $                 fluxes(2,2,1)                                                  
                  SNET(IM,JH)=fluxes(1,1,2)-fluxes(1,2,2)+fluxes(2,1,2)-          
     $                 fluxes(2,2,2)                                                  
                  rrflux(im,jh,1)=fluxes(1,1,2)                                   
                  rrflux(im,jh,2)=fluxes(1,2,2)                                   
                  rrflux(im,jh,3)=fluxes(2,1,2)                                   
                  rrflux(im,jh,4)=fluxes(2,2,2)                                   
                  rrflux(im,jh,5)=fluxes(1,1,1)-fluxes(1,2,1)                     
                  rrflux(im,jh,6)=fluxes(2,2,1)                                   
                                                                          
                  DO l=nl,1,-1                                                    
c  bottom heating rate is zero in morecret                                
                     LD=NL+1-L                                                      
                     IM=I+IOFM                                                      
                     HTNETO=HTNET(IHem,JH,I,LD)                                     
C           htnet(ihem,jh,i,ld)=htlw(l+1)+htsw(l+1)                        
C ER Modif to turn off radiative heating for first 1 days (if not 1D and not PORB gt 0)
                     IF (((KOUNT/ITSPD).LE.1).AND.(.NOT.L1DZENITH)) THEN
                        IF (PORB.EQ.0) THEN
                           htnet(ihem,jh,i,ld)=0.
                        ENDIF
                     ELSE
                        htnet(ihem,jh,i,ld)=htlw(l+1)+htsw(l+1)                        
                     ENDIF

!KM To force net heating of bottom atmosphere layer to relax to flux temperature
C           if(ld.eq.nl) htnet(ihem,jh,i,ld)=(T(1)-T(2))/BOTRELAXTIME

!KM To force net heating of bottom atmosphere layer to be same as just above
!  if  BOTRELAXTIME < 0
C           if(ld.eq.nl.AND.BOTRELAXTIME.LT.0.0) 
C     &          htnet(ihem,jh,i,ld)=htnet(ihem,jh,i,ld-1)
                                                                          
c sets this heating rate                                                  
                     TTRD(IM,LD)=(HTNETO                                            
     $                    +HTNET(IHEM,JH,I,LD))/(CHRF*2.)                               
                                                                          
c  put in linear interpolation of heating rates between this              
c  longitude and last one calculated (i-nskip)                            
                     IF ((i-ilast.gt.1).and.(nskip.gt.0)) then                                         
c                        write(*,*),i,last
                        DO j=ilast+1,i-1                                              
                           a=REAL(j-ilast)/REAL(i-ilast)                               
                           b=1.-a
c                           write(*,*),IHEM,JH,ilast,ld
c                           write(*,*),HTNET(IHEM,JH,J,LD) 
                 write(*,*),'CANNOT SKIP LONGITUDES IN PARALLEL!! ABORT'
                 write(*,*),'Please set nskip=0 in fort.7'
                 STOP
                           HTNETO=HTNET(IHEM,JH,J,LD)                                  
                           htnet(ihem,jh,j,ld)=a*htnet(ihem,jh,i,ld)+                  
     $                          b*htnet(ihem,jh,ilast,ld)                  
                           im=j+iofm                                                   
                           TTRD(IM,LD)=(HTNETO                                         
     $                          +HTNET(IHEM,JH,J,LD))/(CHRF*2.)                  
                           IF (l.eq.nl) then                                           
                              pnet(im,jh)=a*pnet(i+iofm,jh)+                            
     $                             b*pnet(ilast+iofm,jh)                          
                              snet(im,jh)=a*snet(i+iofm,jh)+                            
     $                             b*snet(ilast+iofm,jh)                          
                              DO k=1,6                                                  
                                 rrflux(im,jh,k)=a*rrflux(i+iofm,jh,k)                   
     $                                +b*rrflux(ilast+iofm,jh,k)               
                              ENDDO                                                     
                           ENDIF                                                       
                        ENDDO                                                         
                     ENDIF                                                          
                  ENDDO                                                           

                  ilast=i                                                         
C end of conditional execution of morcrette code                          
               ENDIF 
C          PRINT *, 'threads number =',TID
C end of loop over longitudes                                            
c               test_wctime=omp_get_wtime()-test_wctime
c               write (*,*) 'DEBUG: thread ',omp_get_thread_num(),
c     +                     ' iteration i=',i,
c     +                     ' time=',test_wctime
            ENDDO
!$OMP END PARALLEL DO
C END OF PARALLEL PROCESSING REGION
c            write (*,*) 'DEBUG: parallel loop is over'
!           write(*,*) 'line 654 of cmorc'
c NOTE TO MIKE ROMAN... what is going on below? and why is ilast
c sometimes different?
            IF (nskip.ne.0) then
               write(*,*),'CANNOT SKIP LONGITUDES IN PARALLEL!! ABORT'
               write(*,*),'Please set nskip=0 in fort.7'
               STOP 
            ELSE
                ilast=mg
            ENDIF

            IF (ilast.ne.mg) then                                             
c               write(*,*) 'ilast',ilast
Cm                stop
                DO j=ilast+1,mg                                                  
                  a=REAL(j-ilast)/REAL(mg+1-ilast)                               
                  b=1.-a                                                         
                  im=j+iofm                                                      
                  DO l=nl,1,-1                                                   
                     ld=nl+1-l                                                    
                     HTNETO=HTNET(IHEM,JH,J,LD)                                   
c                     write(*,*),'l,dl,nl,ihem,jh,j',l,dl,nl,ihmem,jh,j
c                     write(*,*),'htneto,chrf,a,b',htneto,chrf,a,b
                     htnet(ihem,jh,j,ld)=a*htnet(ihem,jh,1,ld)+                   
     $                    b*htnet(ihem,jh,ilast,ld)                  
                     TTRD(IM,LD)=(HTNET(IHEM,JH,J,LD)                             
     $                    +HTNETO)/(CHRF*2.)                              
                     IF (l.eq.nl) then                                            
                        pnet(im,jh)=a*pnet(1+iofm,jh)+                             
     $                       b*pnet(ilast+iofm,jh)                          
                        snet(im,jh)=a*snet(1+iofm,jh)+                             
     $                       b*snet(ilast+iofm,jh)                          
                        DO k=1,6                                                   
                           rrflux(im,jh,k)=a*rrflux(1+iofm,jh,k)                     
     $                          +b*rrflux(ilast+iofm,jh,k)                
                        ENDDO                                                      
                     ENDIF                                                        
                  ENDDO                                                          
               ENDDO                                                            
            ENDIF                                                             
c       if (ihem.eq.2) print *, 'rad ',jh,(pnet(im1,jh),im1=1,IGC)        

         ELSE                   ! ntstep requirement
C                                                                         
c  Doesn't do rad scheme (simply uses old heating rates)                  
Cm           write(*,*) 'Ilast,mg',ilast,mg                                                              
            DO i=1,mg
Cm               write(*,*) 'line 691'                                                        
               DO LD=1,NL                                                     
                  im=i+IOFM                                                    
                  TTRD(im,LD)=(htnet(ihem,jh,i,ld))/CHRF                        
               ENDDO                                                          
            ENDDO                                                            
         ENDIF                                                              
         IOFM=MGPP                                                          
 800  CONTINUE                  ! end of loop over hemispheres                             
      IF (LSHORT.AND.(KOUNT.eq.1)) then                                   
         DO l=1,nl                                                         
            DO i=1,igc                                                      
               ttrd(i,l)=ttrd(i,l)*2.                                        
            ENDDO                                                           
         ENDDO                                                             
      ENDIF                                                               
                                                                          
      RETURN                                                              
      END                                                                 
