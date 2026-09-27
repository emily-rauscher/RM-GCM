C**********************************************************               
C             SUBROUTINE VDIFF                                            
C**********************************************************               
      SUBROUTINE VDIFF                                                    
C                                                                         
C     VERTICAL DIFFUSION OF MOMENTUM, HEAT AND MOISTURE.                  
C     CONSTANT DIFFUSION COEFFICIENTS SET IN INITAL.                      
C     ASSUMES ZERO FLUX AT UPPER AND LOWER BOUNDARIES.                    
C     SURFACE FLUX DEALT WITH SEPARATELY IN BLAYER.                       
C                                                                         
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'rcommons.h'
C      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
C     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
C                                                                         
C                                                                         
C     Sets basic constants, especially those needed for array dimensions  
C                                                                         
!       PARAMETER(MH=2,PI=3.14159265359,PI2=2.0*PI                          
!      +,NNP=NN+1,MGPP=MG+2,JGP=JG+1,JGG=JG*NHEM,JGGP=JGG+1,MJP=NWJ2+NWJ2   
!      +,NLM=NL-1,NLP=NL+1,NLPP=NL+2,NLA=NL+3,NLB=NL+4,NL2=NL*NL            
!      +,IDA=(MG+MG+MG)/2+1,IDB=NWJ2*NL,IDC=IDB+IDB,IDD=MGPP*NL             
!      +,IDE=NL2*NN,IDF=NCRAY*(MG+1),IDG=JG*NL,IDH=JG*MG                    
!      +,IDI=NNP/2,IDJ=IDI*IDI,IDK=NL*IDI,IDL=MGPP/2,IDM=NNP/2,IDN=IDM*NL   
!      +,NWW=1+(MM-1)/MOCT)                                                 
!       PARAMETER(IGA=NWJ2*NHEM,IGB=IDB*NHEM,IGC=MGPP*NHEM,IGD=IDD*NHEM     
!      +,IGG=IDG*NHEM,IGL=IDL*NHEM,IGM=IDM*NHEM,IGN=IDN*NHEM                
!      +,IGO=IGA+IGA,IGP=IGB+IGB,NFTWG=(5+NTRAC)*NL+3                       
!      +,NFTGW=(6+3*NTRAC)*NL+2,NFTGD=(3+NTRAC)*NL,NLTR=NL*NTRAC)           
C     Number of 2D (surface) output fields. This value is                 
C     Doubled due to averaged and instantaneous fields.                   
      ! PARAMETER (N2DFLD=21,NGRPAD=N2DFLD*2*IGC)                           
C                                                                         
C                                                                         
C     Basic planetary parameters for run plus information about           
C     vertical grid structure                                             
C                                                                         
C     Note that RD and GASCON are identical and CPD is set from RD,AKAP.  
!       COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)                
!      +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)     
!      +              ,TKP(NL),C(NL2),SQH(NNP)                              
!      +              ,MF,MFP,JZF,NF                                    
!      +              ,AKAP,GA,GASCON,RADEA,WW,PFAC,EZ,AIOCT             
!      +              ,RD,RV,CPD,CLATNT                                     
!      +              ,P0,LRSTRT,LSHORT,LTVEC,LSTRETCH                         
!      +              ,LFLUX                                                
!      +              ,LBALAN,LRESTIJ                                       
!      +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
!      +              ,LNOISE,NFP                                               
!       COMPLEX EZ,AIOCT                                                    
!       LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ                 
!      +       ,LFLUX,LNOISE                                                
!      +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                           
C                                                                         
C           
C                                                                         
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),F(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
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

!        COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
!      & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
!      & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
!      & PORB, OBLIQ, ECCEN 
      
!        LOGICAL LPLOTMAP

!        COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
!      & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
!      & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM(3), AEROSOLS,ABSSW, ABSLW,
!      & ALBSW, NEWTB, NEWTE,RAYPERBARCONS(3),with_TiO_and_VO, opacity_method

!        LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
!      + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS
!        CHARACTER(len=6) :: opacity_method

       COMMON/MAG/ BFIELD,TDRAG_MIN,RAMPUP,LBDRAG
       LOGICAL LBDRAG
                                                   
C                                                                         
      REAL BVP(MG),BQP(MG),BTP(MG),BVM(MG),BQM(MG),BTM(MG)                
      REAL FA(NL)                                                         
C                                                                         
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
!KM Modif: Add new arrays to add new tendencies
      REAL*8 TTVD_RD(IGC,NL), BTP_RD(MG), BTM_RD(MG) ! Radiative diffusion
C      REAL*8 TTVD_NR(IGC,NL) ! Newtonian Relaxation
C      REAL*8 EKTH, EKV, TINT, TIRR, RAMPUP
! Magnetic drag (U, V, and T tendencies)
      REAL*8 UTVD_BDRAG(IGC,NL), VTVD_BDRAG(IGC,NL), TTVD_BDRAG(IGC,NL)
      REAL*8 RHO_CGS, ANG_FACTOR ! in Gauss

C                                                                         
      DO 10 L=1,NL                                                        
   10 FA(L)=1.0/DSIGMA(L)                                                 
C                                                                         
      DO 800 IHEM=1,NHEM                                                  
        IOFM=(IHEM-1)*MGPP                                                
        DO 20 I=1,MG                                                      
          J=I+IOFM                                                        
          TLPH=TG(J,2)+TG(J,1)                                            
          FTSQ=FB(1)/(TLPH*TLPH)                                          
          BVP(I)=AKVV*FTSQ                                                
          BQP(I)=AKQV*FTSQ                                                
          BTP(I)=AKTV*FTSQ   
          BTP_RD(I)=0.!*5.e2/(RADEA*RADEA*WW)*FTSQ !KM


C Scheme for radiation diffusion turned off by zeroing BTP_RD and FBOT
          FBOT= 0.!*1e-6 !KM

          UTVD(J,1)=FA(1)*BVP(I)*(UG(J,2)-UG(J,1))                        
          VTVD(J,1)=FA(1)*BVP(I)*(VG(J,2)-VG(J,1))                        
          QTVD(J,1)=FA(1)*BQP(I)*(QG(J,2)-QG(J,1))                        
          TTVD(J,1)=FA(1)*BTP(I)*(SK(1)*TG(J,2)-TG(J,1))                 
          TTVD_RD(J,1)=FA(1)*BTP_RD(I)*(TG(J,2)-TG(J,1)) - 1e-2*FBOT !KM 
C         TTVD_NR(J,1)=
C     &(TGUI(P0*SIGMA(1),TINT,GA*1.e2,TIRR,EKTH,EKV,1.0,0.5)/CT
C     & -TG(J,1))/
C     & TAURAD(P0*SIGMA(1),CT*TG(J,1),GA*1.e2,EKTH,GASCON/AKAP)/WW !KM
! NB: P0 is used above, rather than actual surface pressure.Check SPG values

! Force T-P profile for 1 day, on relaxation time of 0.1 day, then stop
C           IF ((KOUNT/ITSPD).LE.1) THEN
C         TTVD_NR(J,1)=
C     &(TGUI(P0*SIGMA(1),TINT,GA*1.e2,TIRR,EKTH,EKV,1.0,0.5)/CT
C     & -TG(J,1))/0.1
C          ELSE
C         TTVD_NR(J,1)=0.0     
C           ENDIF

! Magnetic Drag Part
           IF((LBDRAG).AND.((KOUNT/ITSPD).GT.RAMPUP)) THEN
              ALAT1=ALAT(JH)*REAL(-(ihem*2.)+3)
              ANG_FACTOR=abs(cos(3.141592/2.0 *(1.0- ALAT1/90.0)))
              RFAC=MIN(1., (KOUNT/ITSPD)/RAMPUP-1.)

!      write(*,*) JH, ALAT1, cos(3.141592/2.0/90.0*ALAT1), ANG_FACTOR
!      pause
               
              RHO_CGS=P0*(SPG(J)+1.)*SIGMA(1)/(GASCON*CT*TG(J,1))*1.0e-3
              UTVD_BDRAG(J,1)=-UG(J,1)/MAX(TMAG(RHO_CGS,CT*TG(J,1),
     & BFIELD,GASCON,METALLICITY)/ANG_FACTOR*WW,TDRAG_MIN)*RFAC
              VTVD_BDRAG(J,1)=-VG(J,1)/MAX(TMAG(RHO_CGS,CT*TG(J,1),
     & BFIELD,GASCON,METALLICITY)/ANG_FACTOR*WW,TDRAG_MIN)*RFAC*0.0 !No Vdrag
              TTVD_BDRAG(J,1)=(-UTVD_BDRAG(J,1)*UG(J,1)
     & -VTVD_BDRAG(J,1)*VG(J,1))*(CG/(GASCON/AKAP)/CT)  !RFAC already here
           ELSE
              UTVD_BDRAG(J,1)=0.0
              VTVD_BDRAG(J,1)=0.0
              TTVD_BDRAG(J,1)=0.0
           ENDIF

          UTVD(J,1)= UTVD(J,1) + UTVD_BDRAG(J,1) !KM
          VTVD(J,1)= VTVD(J,1) + VTVD_BDRAG(J,1) !KM
          TTVD(J,1)=TTVD(J,1)+TTVD_RD(J,1)+TTVD_BDRAG(J,1)!+TTVD_NR(J,1)!KM

   20   CONTINUE                                                          
        DO 30 L=2,NLM                                                     
          LP=L+1                                                          
          LM=L-1                                                          
          DO 30 I=1,MG                                                    
            J=I+IOFM                                                      
            BVM(I)=BVP(I)                                                 
            BQM(I)=BQP(I)                                                 
            BTM(I)=BTP(I) 
            BTM_RD(I)=BTP_RD(I) !KM 
                                                
            TLPH=TG(J,LP)+TG(J,L)                                         
            FTSQ=FB(L)/(TLPH*TLPH)                                        
            BVP(I)=AKVV*FTSQ                                              
            BQP(I)=AKQV*FTSQ                                              
            BTP(I)=AKTV*FTSQ 

C Use simple T^5/P^2 scaling for radiative diffusion coefficient
C            BTP_RD(I)=(5.e1*((TG(J,LP)+TG(J,L))/TG(J,2)+TG(J,1))**5.0 /
C     &   (SIGMA(L)/SIGMA(1))**2.0)/(RADEA*RADEA*WW)*FTSQ !KM

           BTP_RD(I)=(0.*5.e2*((TG(J,LP)+TG(J,L))/TG(J,2)+TG(J,1))**0.0/
     &   (SIGMA(L)/SIGMA(1))**0.0)/(RADEA*RADEA*WW)*FTSQ !KM

C            write(*,*) L, I, BTP_RD(I)

C            IF (L.EQ.NLM) THEN                                            
C              if (LBL.and.(.not.LOLDBL)) THEN                             
!KM Modif: remove extra momentum diffusion (x5 bottom layer)
C                BVP(I)=BVP(I)*5.                                          
CC                BQP(I)=BQP(I)*5.                                         
CC                BTP(I)=BTP(I)*5.                                         
C              endif                                                       
C            ENDIF                                                         

C            IF (L.LT.(NLM/2)) THEN                                            
C!KM Modif: remove diffusion in top half of domain
C                BVP(I)=BVP(I)*0.                                          
C                BQP(I)=BQP(I)*0.                                         
C                BTP(I)=BTP(I)*0.                                         
C            ENDIF                                                         

            UTVD(J,L)=FA(L)*( BVP(I)*(UG(J,LP)-UG(J,L))                   
     :                 -BVM(I)*(UG(J,L)-UG(J,LM)))                        
            VTVD(J,L)=FA(L)*( BVP(I)*(VG(J,LP)-VG(J,L))                   
     :                 -BVM(I)*(VG(J,L)-VG(J,LM)))                        
            QTVD(J,L)=FA(L)*( BQP(I)*(QG(J,LP)-QG(J,L))                   
     :                 -BQM(I)*(QG(J,L)-QG(J,LM)))                        
            TTVD(J,L)=FA(L)*( BTP(I)*(SK(L)*TG(J,LP)-TG(J,L))             
     :                 -BTM(I)*(TG(J,L)-TG(J,LM)/SK(LM)))                 

            TTVD_RD(J,L)=FA(L)*( BTP_RD(I)*(TG(J,LP)-TG(J,L)) !KM            
     :                 -BTM_RD(I)*(TG(J,L)-TG(J,LM)))  

C         TTVD_NR(J,L)=
C     & (TGUI(P0*SIGMA(L),TINT,GA*1.e2,TIRR,EKTH,EKV,1.0,0.5)/CT
C     &  -TG(J,L))/
C     & TAURAD(P0*SIGMA(L),CT*TG(J,L),GA*1.e2,EKTH,GASCON/AKAP)/WW !KM
! NB: P0 is used above, rather than actual surface pressure.Check SPG values

! Force T-P profile for 1 day, on relaxation time of 0.1 day, then stop
C           IF ((KOUNT/ITSPD).LE.1) THEN
C         TTVD_NR(J,L)=
C     &(TGUI(P0*SIGMA(L),TINT,GA*1.e2,TIRR,EKTH,EKV,1.0,0.5)/CT
C     & -TG(J,L))/0.1
C          ELSE
C         TTVD_NR(J,L)=0.0     
C           ENDIF

! Magnetic Drag Part
            IF((LBDRAG).AND.((KOUNT/ITSPD).GT.RAMPUP)) THEN
               ALAT1=ALAT(JH)*REAL(-(ihem*2.)+3)
               ANG_FACTOR=abs(cos(3.141592/2.0 *(1.0- ALAT1/90.0)))
               RFAC=MIN(1., (KOUNT/ITSPD)/RAMPUP-1.)
              
             RHO_CGS=P0*(SPG(J)+1.)*SIGMA(L)/ (GASCON*CT*TG(J,L))*1.0e-3
               UTVD_BDRAG(J,L)=-UG(J,L)/MAX(TMAG(RHO_CGS,CT*TG(J,L),
     & BFIELD,GASCON,METALLICITY)/ANG_FACTOR*WW,TDRAG_MIN)*RFAC
              VTVD_BDRAG(J,L)=-VG(J,L)/MAX(TMAG(RHO_CGS,CT*TG(J,L),
     & BFIELD,GASCON,METALLICITY)/ANG_FACTOR*WW,TDRAG_MIN)*RFAC*0.0 !No Vdrag
              TTVD_BDRAG(J,L)=(-UTVD_BDRAG(J,L)*UG(J,L)
     & -VTVD_BDRAG(J,L)*VG(J,L))*(RADEA*WW*RADEA*WW/(GASCON/AKAP)/CT)
           ELSE
              UTVD_BDRAG(J,L)=0.0
              VTVD_BDRAG(J,L)=0.0
              TTVD_BDRAG(J,L)=0.0
           ENDIF

c           if ((j.eq.120).and.(l.eq.10).and.((kount/itspd).gt.rampup))
c     &          then
c              write(*,*) j,ug(j,l),vg(j,l),tg(j,l)
c              write(*,*) utvd(j,l),vtvd(j,l),ttvd(j,l)
c              write(*,*) utvd_bdrag(j,l),vtvd_bdrag(j,l),ttvd_bdrag(j,l)
c           endif

          UTVD(J,L)= UTVD(J,L) + UTVD_BDRAG(J,L) !KM
          VTVD(J,L)= VTVD(J,L) + VTVD_BDRAG(J,L) !KM
          TTVD(J,L)=TTVD(J,L)+TTVD_RD(J,L)+TTVD_BDRAG(J,L)!+TTVD_NR(J,L)!KM


               
C            TTVD(J,L)=TTVD(J,L)+TTVD_RD(J,L) + TTVD_NR(J,L)!KM




   30   CONTINUE                                                          
        DO 40 I=1,MG                                                      
          J=I+IOFM                                                        
          UTVD(J,NL)=-FA(NL)*BVP(I)*(UG(J,NL)-UG(J,NLM))                  
          VTVD(J,NL)=-FA(NL)*BVP(I)*(VG(J,NL)-VG(J,NLM))                  
          QTVD(J,NL)=-FA(NL)*BQP(I)*(QG(J,NL)-QG(J,NLM))                  
          TTVD(J,NL)=-FA(NL)*BTP(I)*(TG(J,NL)-TG(J,NLM)/SK(NLM)) 

         
          TTVD_RD(J,NL)=-FA(NL)*BTP_RD(I)*(TG(J,NL)-TG(J,NLM)) +FBOT !KM
C       TTVD_NR(J,NL)=
C     & (TGUI(P0*SIGMA(NL),TINT,GA*1.e2,TIRR,EKTH,EKV,1.0,0.5)/CT
C     & -TG(J,NL))
C     & /TAURAD(P0*SIGMA(NL),CT*TG(J,NL),GA*1.e2,EKTH,GASCON/AKAP)/WW !KM
! NB: P0 is used above, rather than actual surface pressure.Check SPG values
 
! Force T-P profile for 1 day, on relaxation time of 0.1 day, then stop
C           IF ((KOUNT/ITSPD).LE.1) THEN
C         TTVD_NR(J,NL)=
C     &(TGUI(P0*SIGMA(NL),TINT,GA*1.e2,TIRR,EKTH,EKV,1.0,0.5)/CT
C     & -TG(J,NL))/0.1
C          ELSE
C         TTVD_NR(J,NL)=0.0     
C           ENDIF


! Magnetic Drag Part
            IF((LBDRAG).AND.((KOUNT/ITSPD).GT.RAMPUP)) THEN
               ALAT1=ALAT(JH)*REAL(-(ihem*2.)+3)
               ANG_FACTOR=abs(cos(3.141592/2.0 *(1.0- ALAT1/90.0)))
               RFAC=MIN(1., (KOUNT/ITSPD)/RAMPUP-1.)

              RHO_CGS=P0*(SPG(J)+1.)*SIGMA(NL)/ 
     &              (GASCON*CT*TG(J,NL))*1.0e-3
              UTVD_BDRAG(J,NL)=-UG(J,NL)/MAX(TMAG(RHO_CGS,CT*TG(J,NL),
     & BFIELD,GASCON,METALLICITY)/ANG_FACTOR*WW,TDRAG_MIN)*RFAC
              VTVD_BDRAG(J,NL)=-VG(J,NL)/MAX(TMAG(RHO_CGS,CT*TG(J,NL),
     & BFIELD,GASCON,METALLICITY)/ANG_FACTOR*WW,TDRAG_MIN)*RFAC*0.0 !No Vdrag
              TTVD_BDRAG(J,NL)=(-UTVD_BDRAG(J,NL)*UG(J,NL)
     & -VTVD_BDRAG(J,NL)*VG(J,NL))*(RADEA*WW*RADEA*WW/(GASCON/AKAP)/CT)
           ELSE
              UTVD_BDRAG(J,NL)=0.0
              VTVD_BDRAG(J,NL)=0.0
              TTVD_BDRAG(J,NL)=0.0
           ENDIF


          UTVD(J,NL)= UTVD(J,NL) + UTVD_BDRAG(J,NL) !KM
          VTVD(J,NL)= VTVD(J,NL) + VTVD_BDRAG(J,NL) !KM
          TTVD(J,NL)=TTVD(J,NL)+TTVD_RD(J,NL)+TTVD_BDRAG(J,NL) !KM
C     &         +TTVD_NR(J,NL)


C          TTVD(J,NL)=TTVD(J,NL)+TTVD_RD(J,NL)  +TTVD_NR(J,NL)!KM
         
   40   CONTINUE                                                          
        DO 50 L=1,NL                                                      
          DO 50 I=1,MG                                                    
            J=I+IOFM                                                      
            UG(J,L)=UG(J,L)+DELT2C*UTVD(J,L)                              
            VG(J,L)=VG(J,L)+DELT2C*VTVD(J,L)                              
            QG(J,L)=QG(J,L)+DELT2C*QTVD(J,L)                              
            TG(J,L)=TG(J,L)+DELT2C*TTVD(J,L)    

                  
   50   CONTINUE                                                          
  800 CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
      subroutine abort                                                    
      stop 1                                                              
      end                                                                 


C------ Supplementary functions

CC Function returning Guillot T-P profiles
C      REAL*8 FUNCTION TGUI(press, tint, g, tirr, kth, kv, cosa, f)
C      
C      real*8 press, tint, g, tirr, kth, kv, cosa, f
C
CC      write(*,'(6e13.5)') press, tint, g, tirr, kth, kv, cosa, f
C
C
C
C      gamma=kv/kth
C      press=press/1e5 !Switch from Pascals to bars
C
C      tgui=( 0.75*Tint**4*(2./3. + (kth/g)*press*1.e6) 
C     &   + 0.75*Tirr**4*f*(2./3. + 1./(gamma*sqrt(3.)) + 
C     &      (gamma/sqrt(3.) - 1./(gamma*sqrt(3.))) 
C     &      *exp(-gamma*kth/g*press*1.e6*sqrt(3.))))**0.25
C
CC      write(*,'(6e13.5)') tgui
CC      pause
C      RETURN
C      END

C Function returning radiative relaxation time (diffusion approx)
C      REAL*8 FUNCTION TAURAD(press, temp, gra, kth, Cp)
      
C      real*8 press, temp, gra, kth, Cp


C      Cp=Cp*1e4 !Scale from MKS to CGS
C      press=press/1e5 !Scale from Pascals to bars

C       taurad= 3.*((press*1.e6)**2.0)*kth*Cp/
C     & (16*5.67e-5*(temp**3.0)*(gra**2.0))

C       taurad=max(taurad,1.e-1)
CC       taurad=1e0
CC       taurad=3e-5

CC      write(*,'(6e13.5)')  press, temp, gra, kth, Cp, taurad
      
C      RETURN
C      END

C Function returning magnetic drag time
      REAL*8 FUNCTION TMAG(dens, temp, bfield,gascon,METALLICITY)
 
      real*8 dens,temp,bfield ! in CGS
      real*8 gascon           ! in SI
      real*8 eta
      real*8 METALLICITY
      
C      eta=find_etasimp(dens,temp)
      eta=find_eta(dens,temp,gascon,METALLICITY)
      tmag=12.566371*dens*eta/bfield/bfield

C      write(*,*)'TMAG(s)=',tmag

      RETURN
      END



C Function returning magnetic resistivity
      REAL*8 FUNCTION FIND_ETAsimp(dens, temp)
 
      real*8 dens,temp ! in CGS
      real*8 xe

! For Potassium only
      xe=6.47e-13 * sqrt(1e-7/1e-7) *(temp/1e3)**0.75 * 
     & sqrt(2.4e15/dens*1.6726e-24) * exp(-25188./temp)/1.15e-11

      find_etasimp=230*sqrt(temp)/xe


      RETURN
      END


C Function returning magnetic resistivity
      REAL*8 FUNCTION FIND_ETA(Rho, Tref,gascon,METALLICITY)
 
      real*8 Rho,Tref ! in CGS
      real*8 xe
      real*8 mu, k, h3, me, mp, e, Cl, Pi, ndens, twopime3_2
      real*8 f(28), c(28), of(28)  ! f:abundances H-Ni ; c:first ioniz. potent in eV
      real*8 ni(28), Ki(28), xi(28)
      real*8 METALLICITY

      data f/2.66e10, 1.8e9, 60, 1.2, 9, 1.11e7, 2.31e6, 1.84e7, 
     & 780, 2.6e6, 6.0e4, 1.06e6, 8.5e4, 1.0e6, 6500, 5.0e5, 4740, 
     & 1.06e5, 3500, 6.25e4, 31, 2400, 254, 1.27e4, 9300, 9.0e5, 
     & 2200, 4.78e4/  ! solar abundances

      data c/13.6, 24.6, 5.4, 9.32, 8.30, 11.26, 14.54, 13.61, 
     & 17.42, 21.56, 5.14, 7.64, 5.98, 8.15, 10.55, 10.36, 13.01, 
     & 15.76, 4.34, 6.11, 6.56, 6.83, 6.74, 6.76, 7.43, 7.90, 
     & 7.86, 7.63/ ! first ionization potentials (eV)

      ! If non-solar metallicity, adjust the abundances accordingly
      IF (f(3) .eq. 60) THEN
            IF (METALLICITY .NE. 0.0) THEN
            write(*,*) 'using non-solar abundances for magnetic drag:'
            DO I=28, 3, -1
                  f(i) = f(i) * 10**METALLICITY
                  f(i) = f(i) / f(1)
            ENDDO
            f(2) = f(2)/f(1)
            f(1) = 1.0
            WRITE(*,*) f
            ELSE
            write(*,*) 'using solar abundances for magnetic drag:'
            DO I=28, 1, -1
                  f(i) = f(i) / f(1)
            ENDDO
            write(*,*) f
            ENDIF
            DO I=28, 1, -1
            c(i) = c(i) * 11605 ! converts from eV to deg K, but only once
            ENDDO
      ENDIF

C      data f/2.66e10, 1.8e9, 180., 3.6, 27., 3.33e7, 6.93e6, 5.52e7, 
C     & 2340., 7.8e6, 1.8e5, 3.18e6, 2.55e5, 3.0e6, 1.95e4, 1.5e6, 1.422e4, 
C     & 3.18e5, 1.05e4, 1.875e5, 93., 7200., 762., 3.81e4, 2.79e4, 2.7e6, 
C     & 6600., 1.434e5/  ! 3 x solar abundances

      mu = 8.3145/gascon*1.e3/6.022e23  ! cgs mass of particle
      k = 1.38e-16  ! erg/K
      h3 = 2.9011753e-79
      me = 9.11e-28
      twopime3_2 = 4.330593e-40
      mp = 1.67e-24  ! not used
      e = 4.8e-10
      Cl = 3e10 ! speed o light 
      Pi = 3.141592

      ndens=Rho/mu

      xe=0.0
      DO i=28, 1, -1
C         c(i)=c(i)*11605 ! ev -> deg K
         ni(i)=ndens*f(i)
         Ki(i)=(1.0/(ni(i)))*(twopime3_2)*
     & ((k*Tref)**(3.0/2.0)) *(exp(-c(i)/Tref)) / h3
         xi(i)=(Ki(i)/(1.+Ki(i)))**0.5 
C         write(*,*) i,f(i),xi(i),230*sqrt(Tref)/xi(i)
         xe=xe+f(i)*xi(i)
C         write(*,'(i4,6e13.5)') i,c(i),f(i),ni(i),Ki(i),xi(i),xe
       ENDDO

C         pause

      find_eta=230*sqrt(Tref)/xe

C      write(*,'(6e13.5)') Rho, Tref, 230*sqrt(Tref)/xe

C      write(*,*) Rho, Tref, xe, 6.47e-13 * sqrt(1e-7/1e-7) 
C     & *(tref/1e3)**0.75 * 
C     & sqrt(2.4e15/rho*1.6726e-24) * exp(-25188/tref)/1.15e-11
C      pause

C      test_xe=6.47e-13 * sqrt(1e-7/1e-7) 
C     & *(tref/1e3)**0.75 * 
C     & sqrt(2.4e15/rho*1.6726e-24) * exp(-25188/tref)/1.15e-11
C      IF(xe.gt.10.0*test_xe.OR.xe.lt.0.1*test_xe) THEN
C        write(*,'(6e13.5)') Rho, Tref, xe, test_xe
C         pause
C      ENDIF

      RETURN
      END


