C**********************************************************               
C             SUBROUTINE SFCT                                             
C**********************************************************               
C   Subroutine to find new T and Q based on DOY and climatology           
      SUBROUTINE SFCT(PLG,JH,IFIRST,TROPHT)                               
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
       COMMON/GSG/GSG(IGC,JG)                                             
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
      INTEGER JH,IFIRST                                                   
      REAL TSTAR1(IGC,JG),TSTAR2(IGC,JG),PLG(IGC)                         
      integer cmth                                                        
      character*3 mn(13)                                                  
      character*31 sstfile                                                
      integer mnlen(13)                                                   
      data mn/'dec','jan','feb','mar','apr','may','jun','jul','aug'       
     &        ,'sep','oct','nov','dec'/                                   
      data mth1/0/                                                        
      real shdum                                                          
      real tropht(mg,nhem,jg)                                             
      save mth1                                                           
c                                                                         
C                                                                         
      RHSCL=288.0*GASCON/GA                                               
      CALL CALNDR(DOY,MTH2,AMFRAC)                                        
      IF (IFIRST.EQ.1) THEN                                               
C                                                                         
       IF (JH.EQ.1) THEN                                                  
        cmth=mth2+1                                                       
        IF (cmth.eq.13) cmth=1                                            
        write(sstfile,'(a,i3.3,a,i3.3,a)')                                
     &      'tqpap_'//mn(cmth)//'_jg',jg,'_mg',mg,'_nhem2.dat'            
        open(unit=45,file='climdata/'//sstfile,status='old')              
c                                                                         
       ENDIF                                                              
C read into tstar2 because changes below to tstar1                        
       DO i2=1,mg                        ! NH                             
         read(45,*)sst,ssq,ssp,shdum,tropht(i2,1,jh)                      
         tstar2(i2,jh)=sst/ct                                             
       ENDDO                                                              
c                                                                         
       DO i2=1,mg                        ! SH                             
          IF (nhem.eq.2) then                                             
            read(45,*)sst,ssq,ssp,shdum,tropht(i2,2,jh)                   
            tstar2(i2+mgpp,jh)=sst/ct                                     
          ELSE                                                            
            read(45,*)         ! because climatology is global            
          ENDIF                                                           
       ENDDO                                                              
c                                                                         
       IF (JH.EQ.JG) THEN                                                 
         CLOSE(45)                                                        
         write(2,*)'CLIMATOLOGICAL TROPOPAUSE HEIGHT'                     
         do k=1,nhem                                                      
           write(2,*)(tropht(1,k,j),j=1,jg)                               
         enddo                                                            
       ENDIF                                                              
                                                                          
       MTH1=CMTH                                                          
      ENDIF  ! IFIRST                                                     
c                                                                         
      IF (MTH2.NE.MTH1) THEN   !time to change month                      
                                                                          
       IF (JH.EQ.1) THEN                                                  
        cmth=mth2+1                                                       
        if (cmth.eq.13) cmth=1                                            
        write(sstfile,'(a,i3.3,a,i3.3,a)')                                
     &      'tqpap_'//mn(cmth+1)//'_jg',jg,'_mg',mg,'_nhem2.dat'          
        open(unit=46,file='climdata/'//sstfile,status='old')              
       ENDIF                                                              
                                                                          
c transfer tstar2 to tstar1 for latitude                                  
       DO i2=1,igc                                                        
         tstar1(i2,jh)=tstar2(i2,jh)                                      
       ENDDO                                                              
c read in new tstar2                                                      
       DO i2=1,mg                        ! NH                             
         read(46,*)sst,ssq,ssp,shdum,tropht(i2,1,jh)                      
         tstar2(i2,jh)=sst/ct                                             
       ENDDO                                                              
c                                                                         
       DO i2=1,mg                        ! SH                             
         IF (nhem.eq.2) then                                              
           read(46,*)sst,ssq,ssp,shdum,tropht(i2,2,jh)                    
           tstar2(i2+mgpp,jh)=sst/ct                                      
         ELSE                                                             
           read(46,*)  ! because climatology is global                    
         ENDIF                                                            
       ENDDO                                                              
       IF (JH.EQ.JG) THEN                                                 
         CLOSE(46)                                                        
         MTH1=MTH2  ! stop doing loop until next month                    
       ENDIF                                                              
      ENDIF   ! MTH2.NE.MTH1                                              
c                                                                         
c                                                                         
c  set TSTAR and QSTAR                                                    
C                                                                         
C                                                                         
      IF (.NOT.LSL) then                                                  
        DO ihem=1,nhem                                                    
          iofm=(ihem-1)*mgpp                                              
          DO i2=1,mg                                                      
            j=i2+iofm                                                     
            IF (gsg(j,jh).gt.0.) then                                     
              tstar(j,jh)=tstar2(j,jh)*amfrac+                            
     $                 tstar1(j,jh)*(1.0-amfrac)                          
              IF (LPERPET) THEN                                           
                escon=1./exp(-gsg(j,jh)/RHSCL)                            
              ELSE                                                        
                escon=1./plg(j)                                           
              ENDIF                                                       
              sqstar(j,jh)=escon*pqsat(tstar(j,jh))                       
              qstar(j,jh)=sqstar(j,jh)*0.75                               
            endif                                                         
          enddo                                                           
        enddo                                                             
      endif                                                               
      IF (.NOT.LOC) THEN                                                  
        DO ihem=1,nhem                                                    
          iofm=(ihem-1)*mgpp                                              
          DO i2=1,mg                                                      
            j=i2+iofm                                                     
            IF (gsg(j,jh).le.0.) then                                     
              tstar(j,jh)=tstar2(j,jh)*amfrac+                            
     $              tstar1(j,jh)*(1.0-amfrac)                             
              IF (LPERPET) THEN                                           
                escon=1.                                                  
              ELSE                                                        
                escon=1./plg(j)                                           
              ENDIF                                                       
              sqstar(j,jh)=escon*pqsat(tstar(j,jh))                       
              qstar(j,jh)=sqstar(j,jh)                                    
            endif                                                         
          enddo                                                           
        enddo                                                             
      ENDIF                                                               
      RETURN                                                              
      END                                                                 
