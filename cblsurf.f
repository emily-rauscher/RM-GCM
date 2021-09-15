C**********************************************************               
C             SUBROUTINE BLSURF                                           
C**********************************************************               
      subroutine BLSURF                                                   
C-----------------------------------------------------------------------  
C     Perform time-splitting integration of the boundary layer and        
C     Surface scheme together for a whole timestep                        
C     In the case of LSL or LOC being false, should forget the            
C     timesplitting and just do the BL stuff (leave the surface alone)    
C                                                                         
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
       COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     :               ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     :               ,TTSW(IGC,NL),TTLW(IGC,NL)                           
       COMMON/GSG/GSG(IGC,JG)                                             
C                                                                         
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
       REAL OCFLUX(130,16,12)  ! OCEAN HEAT FLUX MONTHLY-AVE              
C-----------------------------------------------------------------------  
C     set up values for reference pressure and various temperatures       
C     in non-dimensionalised form                                         
C-----------------------------------------------------------------------  
      parameter (tzc=273.15)                                              
      logical LSLH,LOCH                                                   
      DATA IFIRST/1/                                                      
C                                                                         
C Lookup table for vegetation to Z0: from Unified Model Doc. no. 70       
C Snow free lookup table                            MMJ April 2000        
C                                                                         
       DIMENSION SZLOOK(24)                                               
      save apfrac,mpth1,mpth2                                             
       DATA SZLOOK /0.001,1.0E-4,3.0E-4,1.0,1.2,1.0,                      
     $ 1.0,1.2,1.0,0.4,0.4,0.4,                                           
     $ 0.12,0.12,0.12,0.12,0.12,0.12,                                     
     $ 0.12,0.12,1.5,0.12,0.12,3.0E-3/                                    
C                                                                         
C   THIS BIT READS IN RUN 1.0 (10 year ocean flux)                        
C 30-3-98  PMF                                                            
      IF (IFIRST.EQ.1.AND.LOC) THEN                                       
        IFIRST=0                                                          
C OCEAN FLUX is Wm-2 +ve Downwards                                        
C FILE PRODUCED BY AVERAGING 10 YEAR UTF SURFACE BALANCE TERM             
C OF RUN 1.0                                                              
C                                                                         
        READ(41,*)                                                        
        READ(41,*) OCFLUX                                                 
        CLOSE(41)                                                         
        PRINT *,' READ IN OCEAN FLUXES'                                   
      ENDIF                                                               
c                                                                         
C  SORT OUT TIMINGS                                                       
c                                                                         
      IF (JH.EQ.1.AND.LOC) THEN                                           
         CALL CALNDR(DOY,Mpth1,ApFRAC)                                    
         MPth2=mpth1+1                                                    
         if (mpth2.eq.13) mpth2=1                                         
      ENDIF                                                               
      tp1=(tzc+1.)/ct                                                     
      tz=tzc/ct                                                           
      tm2=(tzc-2.)/ct                                                     
      tm3=(tzc-3.)/ct                                                     
                                                                          
C     Set up the number of timesplitting levels for land and ocean,       
C     also the timesteps.                                                 
                                                                          
      DELTL=DELT2C/ITSLL                                                  
      DELTO=DELT2C/ITSLO                                                  
                                                                          
                                                                          
C     Set up some more constants for the soil scheme                      
                                                                          
      ssmc23=ssmc*(2./3.)                                                 
      ssmc13=ssmc*(1./3.)                                                 
      srsum=sd1/skse+sd2/skse                                             
                                                                          
C     And some for the BL scheme                                          
                                                                          
      RCSJ=SECSQ(JH)                                                      
      SQRC=SQRT(RCSJ)                                                     
      RSIGF=2.*(1-SIGMA(NL))/(1+SIGMA(NL))                                
C     5cm roughness length over land - note includes factor of g          
      RGZ0L=2.2726E-6                                                     
C     1mm over ocean                                                      
      RGZ0O=4.5453E-8                                                     
                                                                          
                                                                          
C     Now begin looping over gridpoints                                   
                                                                          
      do ihem=1,nhem                                                      
         do i=1,mg                                                        
            j=(ihem-1)*mgpp+i                                             
                                                                          
C                                                                         
C Get land roughness lengths from vegetation index                        
C                                                                         
            RGZ0L=SZLOOK(NINT(SVEGE(J,JH)))/2.2E4                         
C                                                                         
C     Now do a single gridpoint                                           
                                                                          
            LSLH=(LSL.and.(gsg(j,jh).gt.0.))                              
            LOCH=(LOC.and.(gsg(j,jh).le.0.))                              
                                                                          
            if (gsg(j,jh).gt.0.) then                                     
               MAXITER=ITSLL                                              
               DELTU=DELTL                                                
               RGZ0=RGZ0L                                                 
               if (LSL) then                                              
                  TSCUR=TSTARO(J,JH)                                      
                  TDCUR=TDEEPO(J,JH)                                      
               else                                                       
                  TSCUR=TSTAR(J,JH)                                       
                  QCUR=QSTAR(J,JH)                                        
                  SQCUR=SQSTAR(J,JH)                                      
               endif                                                      
            else                                                          
               MAXITER=ITSLO                                              
               DELTU=DELTO                                                
               RGZ0=RGZ0O                                                 
               if (LOC) then                                              
                  TSCUR=TSTARO(J,JH)                                      
               else                                                       
                  TSCUR=TSTAR(J,JH)                                       
                  QCUR=QSTAR(J,JH)                                        
                  SQCUR=SQSTAR(J,JH)                                      
               endif                                                      
            endif                                                         
                                                                          
            QNLGS=0.                                                      
            TNLGS=0.                                                      
            UNLGS=0.                                                      
            VNLGS=0.                                                      
            SHBLS=0.                                                      
            SLBLS=0.                                                      
            BLCDS=0.                                                      
            BLVBS=0.                                                      
            TXBLS=0.                                                      
            TYBLS=0.                                                      
                                                                          
            do iter=1,MAXITER                                             
                                                                          
C     if we have interactive soil, calculate surface humidities           
                                                                          
               if (LSLH) then                                             
                  if (smstar(j,jh).gt.ssmc) then                          
                     srh=1.                                               
                     smstar(j,jh)=ssmc                                    
                  else if ((smstar(j,jh).gt.ssmc23)                       
     $                    .or.(hsnow(j,jh).gt.0.)) then                   
                     srh=1.                                               
                  else if (smstar(j,jh).gt.ssmc13) then                   
                     srh=(smstar(j,jh)/ssmc13)-1.                         
                  else                                                    
                     srh=0.                                               
                  endif                                                   
                                                                          
                  sqcur=pqsat(tscur)/plg(j)                               
                                                                          
                  qcur=srh*sqcur                                          
                                                                          
               endif                                                      
                                                                          
C     Ditto for interactive ocean                                         
                                                                          
               if (LOCH) then                                             
                  sqcur=pqsat(tscur)/plg(j)                               
                  qcur=sqcur                                              
               endif                                                      
                                                                          
               THNL=TG(J,NL)/SKAP(NL)                                     
               THBAR=0.5*(THNL+TSCUR)                                     
               RZZ0=RSIGF*THBAR/RGZ0                                      
               RLZZ02=(LOG(RZZ0))**2                                      
C     CD=(K/LN(Z/Z0))**2, K=0.41                                          
               CD=0.1681/RLZZ02                                           
               blcds=blcds+cd/maxiter                                     
               CTAU=CD/TSCUR                                              
               CSH=CTAU/AKAP                                              
               DRAGJ=CD*DRAG/TSCUR                                        
C     CALCULATE WIND SPEED, INCLUDING MINIMUM                             
               VM=SQRT(RCSJ*(UG(J,NL)*UG(J,NL)+VG(J,NL)*VG(J,NL)))        
               VMP=VM+BLVAD                                               
C     CALCULATE SURFACE STRESS AND LOWEST LEVEL FRICTION                  
               CVM=DRAGJ*VMP                                              
               UNLG(J,NL)=-CVM*UG(J,NL)                                   
               UNLGS=UNLGS+UNLG(J,NL)/MAXITER                             
               VNLG(J,NL)=-CVM*VG(J,NL)                                   
               VNLGS=VNLGS+VNLG(J,NL)/MAXITER                             
               CVM=CTAU*PLG(J)*VMP*SQRC                                   
               TXBL(J)=CVM*UG(J,NL)                                       
               TXBLS=TXBLS+TXBL(J)/MAXITER                                
               TYBL(J)=CVM*VG(J,NL)                                       
               TYBLS=TYBLS+TYBL(J)/MAXITER                                
C     CALCULATE SENSIBLE AND LATENT HEAT FLUXES                           
               QCUR=MIN(SQCUR,MAX(QCUR,QG(J,NL)))                         
               DTH=TSCUR-THNL                                             
               DQ=BLRH*(QCUR-QG(J,NL))                                    
               BLVB(J)=0.                                                 
               IF (DTH.GT.0.) THEN                                        
                  BLVB(J)=5.95*RLZZ02*SQRT(RGZ0*DTH/THBAR)                
C                  BLVB(J)=5.95*RLZZ02*SQRT(RGZ0L*DTH/THBAR)              
C                  if (gsg(j,jh).le.0.) then                              
C                     BLVB(J)=500./464.*SQRT(DTH/THBAR)                   
C                  endif                                                  
C     RESTRICT ENHANCEMENT TO A FACTOR OF 4 OVER NEUTRAL VALUES           
C                  IF (BLVB(J)/VM.GT.4.) BLVB(J)=VM*4.                    
               ENDIF                                                      
C                                                                         
C Cd_h is multiplied by 0.2: I think this underestimates                  
C Cd_h in unstable conditions SO                                          
C MAKE Cd_h = Cd, not 0.2*Cd in unstable conditions                       
C                                                                         
C               VMP=VM+BLVB(J)                                            
               VMP=VM+5.0*BLVB(J)                                         
               CVM=0.2*DRAGJ*VMP                                          
               TNLG(J,NL)=CVM*DTH                                         
               TNLGS=TNLGS+TNLG(J,NL)/MAXITER                             
               QNLG(J,NL)=CVM*DQ                                          
               QNLGS=QNLGS+QNLG(J,NL)/MAXITER                             
               CVM=0.2*CSH*PLG(J)*VMP                                     
               SHBL(J)=CVM*DTH                                            
               SHBLS=SHBLS+SHBL(J)/MAXITER                                
C     CHOOSE APPROPRIATE LATENT HEAT                                      
               IF (TSCUR.GT.TZ) THEN                                      
                  SLBL(J)=CVM*CTQ*DQ                                      
               ELSE                                                       
                  SLBL(J)=CVM*CTQI*DQ                                     
               ENDIF                                                      
               SLBLS=SLBLS+SLBL(J)/MAXITER                                
               UG(J,NL)=UG(J,NL)+DELTU*UNLG(J,NL)                         
               VG(J,NL)=VG(J,NL)+DELTU*VNLG(J,NL)                         
               TG(J,NL)=TG(J,NL)+DELTU*TNLG(J,NL)                         
               QG(J,NL)=QG(J,NL)+DELTU*QNLG(J,NL)                         
                                                                          
               if (LSLH) then                                             
                                                                          
C                  if ((j.eq.18).and.(jh.eq.10))then                      
C                     write (*,*) 'Initial TS=',tscur*                    
C     $                    752.-273.15                                    
C                  endif                                                  
                                                                          
C     CHECK SIGNS of shbl,slbl,snet                                       
                                                                          
                  SFC=-shbl(j)-slbl(j)+snet(j,jh)/(CV*P0)                 
                                                                          
                  rkappa=2./(srsum+hsnow(j,jh)/sksn)                      
                  pc1=shcs                                                
                  pc2=shcs                                                
                  if ((tscur.gt.tm3).and.(tscur.lt.tp1))                  
     $                 pc1=pc1+shcsp                                      
                  if ((tdcur.gt.tm3).and.(tdcur.lt.tp1))                  
     $                 pc1=pc1+shcsp                                      
                  hc1=pc1*sd1+min(hsnow(j,jh),shsmax)*shcsn               
                  hc2=pc2*sd2                                             
                                                                          
                  trans=rkappa*(tscur-tdcur)                              
                  rct1=(sfc-trans)/hc1                                    
                  rct2=trans/hc2                                          
                                                                          
                  tscur=tscur+rct1*deltu                                  
                  tdcur=tdcur+rct2*deltu                                  
                                                                          
                  if (tscur.gt.tz) then                                   
                     if (hsnow(j,jh).gt.0.) then                          
C     The factor of 2 is *REALLY* meant to be here, since we              
C     are dealing with a forward time single level variable (hsnow)       
C     within an (effectively) centre time timestep - thus each            
C     segment of time is covered twice in this routine...                 
                        hsnow(j,jh)=hsnow(j,jh)                           
     $                       -(tscur-tz)*hc1/(slhf*sdsn)/2.               
                        tscur=tz                                          
                     endif                                                
                  endif                                                   
                                                                          
               endif                                                      
                                                                          
               if (LOCH) then                                             
c                                                                         
C INTEROP between months                                                  
          OCF=(1.0-ApFRAC)*OCFLUX(J,JH,MpTH1)+                            
     $        ApFRAC*OCFLUX(J,JH,MpTH2)                                   
c                                                                         
                  if (tscur.lt.tm2) then !ice                             
              if (ocf.gt.0) Then  ! flux corr out of ice                  
                    ocf=0   !ice can't cool to ocean                      
              endif                                                       
                     hc=shci                                              
                  else                                                    
                     hc=shco                                              
                  endif                                                   
c                                                                         
C add OCF term to ocean heat flux                                         
                SFC=-shbl(j)-slbl(j)+(snet(j,jh)-OCF)/(CV*P0)             
c               SFC=-shbl(j)-slbl(j)+(snet(j,jh))/(CV*P0)                 
c                                                                         
                  dts=sfc/hc*DELTU                                        
c                                                                         
                  if ((tscur.lt.tm2).and.(tscur+dts                       
     $                 .ge.tm2)) then                                     
                     tscur=tm2+(tscur+dts-tm2)*shci/shco                  
                  else if ((tscur.ge.tm2).and.(                           
     $                    tscur+dts.lt.tm2)) then                         
                     tscur=tm2+(tscur+dts-tm2)*shco/shci                  
                  else                                                    
                     tscur=tscur+dts                                      
                  endif                                                   
C                                                                         
                  if (tscur.ge.tm2) then                                  
                     salb(j,jh)=sbal(j,jh)                                
                  else                                                    
                     salb(j,jh)=saice                                     
                  endif                                                   
                                                                          
               endif                                                      
            enddo                                                         
                                                                          
C     End of iteration process, done full double timestep.                
                                                                          
            if (LSLH) then                                                
               tdeepo(j,jh)=tdeep(j,jh)*(1.-2.*pnu)                       
     $              +pnu*(tdeepo(j,jh)+tdcur)                             
               tdeep(j,jh)=tdcur                                          
               tstaro(j,jh)=tstar(j,jh)*(1.-2.*pnu)                       
     $              +pnu*(tstaro(j,jh)+tscur)                             
               tstar(j,jh)=tscur                                          
               sqstar(j,jh)=sqcur                                         
               qstar(j,jh)=qcur                                           
            endif                                                         
                                                                          
            if (LOCH) then                                                
               tstaro(j,jh)=tstar(j,jh)*(1.-2.*pnu)                       
     $              +pnu*(tstaro(j,jh)+tscur)                             
               tstar(j,jh)=tscur                                          
               sqstar(j,jh)=sqcur                                         
               qstar(j,jh)=qcur                                           
            endif                                                         
                                                                          
            QNLG(J,NL)=QNLGS                                              
            TNLG(J,NL)=TNLGS                                              
            UNLG(J,NL)=UNLGS                                              
            VNLG(J,NL)=VNLGS                                              
            SHBL(J)=SHBLS                                                 
            SLBL(J)=SLBLS                                                 
            BLCD(J)=BLCDS                                                 
            BLVB(J)=BLVBS                                                 
            TXBL(J)=TXBLS                                                 
            TYBL(J)=TYBLS                                                 
                                                                          
         enddo                                                            
      enddo                                                               
                                                                          
      end                                                                 
