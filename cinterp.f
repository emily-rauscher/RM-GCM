C**********************************************************               
C             SUBROUTINE INTERPF                                          
C**********************************************************               
      subroutine interpf(fclim,fmod,ps)                                   
c                                                                         
c This subroutine takes in a profile from climatology (15 levels)         
c and interpolates to model vertical grid (nl.gt.15), passing back        
c profile fmod.                                                           
c N.B. The fclim profile must correspond to the climatology pressure      
c levels given in this routine (ocliml(15)).                              
c  SMR 03-06-97                                                           
c                                                                         
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
      COMMON/BATS/  BM1(IDE),AK(NNP),AQ(NL2),G(NL2),TAU(NL2)              
     +              ,KOUNT,KITS,KSTART,KTOTAL,KRUN,BEGDAY,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,CTRA(NTRAC),KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C                                                                         
      COMMON/PHYS/  LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ                      
     :              ,TSTAR(IGC,JG),QSTAR(IGC,JG),FRAD(JG,NHEM)            
     :              ,TSTARO(IGC,JG),TDEEPO(IGC,JG),smstar(igc,jg)         
     :              ,tdeep(igc,jg),hsnow(igc,jg),sqstar(igc,jg)           
     :              ,SALB(IGC,JG),SBAL(IGC,JG),BLCD(IGC)                  
     :              ,SVEGE(IGC,JG),CD,DRAG,BLVAD,BLA,BLRH,BLVB(IGC)       
     :              ,AKVV,AKTV,AKQV,ESCONA,ESCONB,EPSIQ,CTQ,CCC           
     : ,ctqi,sdsn,shcs,shcsp,shcsn,skse,sksn,slhf,sd1,sd2,sdw             
     :        ,ssmc,sdsnd,LSL,sasnow,saice,shsstar,shsmax                 
     :     ,LOC,SHCO,SHCI,LNOICE,LOLDBL,LCOND,LNNSK,ITSLL,ITSLO           
     :              ,CCR,RCON,DTBUOY,TSLA,TSLB,TSLC,TSLD,CUT1,CUT2        
     :              ,NLCR,NCUTOP,CURHM,AKTC,AKQC,CUBMT,CBADJT,CBADJP      
     :              ,SKAP(NL),SK(NLM),FWS(NL),CLR(NL),FB(NLM)             
     :              ,TTDC(NL),QTDC(NL),TTMC(NL),QTMC(NL),TC(NL),QC(NL)    
     :              ,CTCR(NL,NHEM),CTLR(NL,NHEM),NAVRD,NAVWT,DELT2C       
      LOGICAL LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ,LSL,LOC                    
     :       ,LNOICE,LOLDBL,LCOND,LNNSK                                   
C                                                                         
C                                                                         
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
                                                                          
      real fclim(15)                ! Climatology profile. Must           
                                    ! correspond to                       
                                    ! pressure levels given below!        
      real fmod(nl)                 ! Model (interpolated) profile.       
      integer l                     ! model level counter.                
      integer lo                    ! climatology level counter.          
      real pressu                   ! pressure in mb of current           
                                    ! level being interpolated.           
      real ocliml(15)               ! ozone climatology levels in mb.     
      real ps                       ! Surface pressure.                   
                                                                          
                                                                          
      save                                                                
      data ocliml/1.0,3.0,10.0,30.0,50.0,70.0,100.0,150.0,200.0,          
     &              250.0,300.0,400.0,500.0,700.0,850.0/                  
                                                                          
                                                                          
                                                                          
                                                                          
c                                                                         
c -------------------------- Start of interpolation                       
                                                                          
      do 224 l=1,nl                  ! Loop over model levels.            
         pressu=sigma(l)*ps/100.0    ! Pressure in mb of model level      
                                     ! where ozone required.              
         do 226 lo=1,14              ! Loop over climatology levels.      
                                                                          
            if (pressu.lt.1.0) then                 ! Top                 
              fmod(l)=fclim(1)                                            
            elseif (pressu.eq.ocliml(lo)) then   ! if on a                
              fmod(l)=fclim(lo)                  ! climatology level.     
            elseif (pressu.gt.ocliml(lo).and.pressu.lt.ocliml(lo+1))      
     &        then                               ! need to interpolate    
              fmod(l)=exp(                                                
     &             log(fclim(lo)) +                                       
     &             (log(fclim(lo+1)) -                                    
     &              log(fclim(lo)))*                                      
     &             (log(pressu) - log(ocliml(lo))) /                      
     &             (log(ocliml(lo+1)) - log(ocliml(lo))))                 
            elseif (pressu.ge.850.0) then         ! Bottom                
              fmod(l)=fclim(15)                                           
            endif                                                         
                                                                          
 226  continue                       ! End of loop over climatology.      
 224  continue                       ! End of loop over model levels.     
c                                                                         
c --------------------------- End of interpolation                        
c                                                                         
      return                                                              
      end                                                                 
