C**********************************************************               
C             SUBROUTINE INISURF                                          
C**********************************************************               
      subroutine INISURF                                                  
C-----------------------------------------------------------------------  
C     Subroutine to initialise the surface model from history on channel  
C     18.                                                                 
C-----------------------------------------------------------------------  
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
       COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     :               ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     :               ,TTSW(IGC,NL),TTLW(IGC,NL)                           
C                                                                         
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
       COMMON/GSG/GSG(IGC,JG)                                             
C                                                                         
C Lookup table for vegetation to albedo: from Unified Model Doc. no. 70   
C MMJ April 2000                                                          
C Snow free lookup table                                                  
       DIMENSION SALOOK(24)                                               
       DATA SALOOK /0.1, 0.75, 0.06, 0.14, 0.12, 0.13,                    
     $ 0.13, 0.13, 0.13, 0.17, 0.16, 0.16,                                
     $ 0.19, 0.2, 0.2, 0.12, 0.17, 0.19,                                  
     $ 0.19, 0.25, 0.18, 0.15, 0.12, 0.35/                                


                                                                          
      DAYNEAR=0.                                                          
 101  read(18,end=200) RCHECK                                             
      if (abs(RCHECK+999.999).gt.0.0001) then                             
          write (2,*) 'Wrong type of surface restart ',                   
     :                   rcheck,abs(rcheck+999.999)                       
          call abort                                                      
      endif                                                               
      read(18) rrkount,rm1tape,rday,CDOY,TSTAR,TDEEP,SMSTAR,QSTAR,        
     $        HSNOW,SQSTAR,SALB,SBAL,TSTARO,TDEEPO,SNET,RM2TAPE           
      if (RM1TAPE.ne.RM2TAPE) then                                        
        write (2,*) 'Surface (SOIL) restart record wrong.'                
        write (2,*) 'RM1TAPE != RM2TAPE',RM1TAPE,RM2TAPE                  
        call abort                                                        
      endif                                                               
      IF (ABS(RDAY-BEGDAY).LE.1.0E-2)THEN                                 
         GOTO 1000                                                        
      ELSE                                                                
         IF (ABS(RDAY-BEGDAY).LT.(ABS(DAYNEAR-BEGDAY)))THEN               
            DAYNEAR=RDAY                                                  
         ENDIF                                                            
         GOTO 101                                                         
      ENDIF                                                               
 200  write (2,*) 'EOF on surface restart file'                           
      write (2,*) 'Looking for day ',begday,' but closest was ',          
     $     daynear                                                        
      stop                                                                
 1000 continue                                                            
C                                                                         
C Read in vegetation data                                                 
C                                                                         
      READ(31)SVEGE                                                       
C                                                                         
      if (CDOY.NE.DOY) then                                               
        write (2,*) 'Surface DOY not equal to dynamical DOY.'             
        write (2,*) cdoy,doy                                              
        write (2,*) 'Value from dynamics used.'                           
      endif                                                               
      PRINT *, 'resetting albedo to SBAL'                                 
      IF(lshort)then                                                      
        DO I=1,IGC                                                        
          DO IH=1,JG                                                      
C Albedo field initialised from vegetation index                          
C The .AND. is because the vegetation index does not contain              
C Antarctica for some unknown reason.                                     
C Find latitude of Antarctica south of 60S                                
                                                                          
           if(i .lt. igc/2)svlat=SI(IH)                                   
           if(i .gt. igc/2)svlat=SI(JGG-IH+1)                             
                                                                          
           if(nint((gsg(i,ih)/(gsg(i,ih)+0.1))) .eq. 1 .and.              
     :        svlat .gt. -0.866)                                          
     :        SBAL(I,IH)=SALOOK(NINT(SVEGE(I,IH)))                        
C Change albedo to cool surface- this is because Tg is                    
C too high (probably because of our clouds) and mucking                   
C around with albedo is better than changing the solar constant           
                                                                          
           SBAL(I,IH) = SBAL(I,IH)+0.05                                   
           SALB(I,IH) = SBAL(I,IH)                                        
        end do                                                            
      end do                                                              
                                                                          
      print*,'Albedos set from Vegetation index'                          
      print*,'Desert value set to',SALOOK(24)                             
        DO j=1,jg                                                         
          DO i=1,igc                                                      
            TSTARO(i,j)=tstar(i,j)                                        
            tdeepo(i,j)=tdeep(i,j)                                        
          ENDDO                                                           
        ENDDO                                                             
      ELSE                                                                
        print*,'Albedos read from RESTART record'                         
      ENDIF                                                               
                                                                          
                                                                          
                                                                          
                                                                          
                                                                          
      END                                                                 
