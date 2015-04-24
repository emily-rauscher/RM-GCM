C**********************************************************               
C             SUBROUTINE MASCOR                                           
C**********************************************************               
      SUBROUTINE MASCOR                                                   
C                                                                         
C     MASCOR - Correction of global average surface pressure.             
C     No longer requires KOUNTE=1 to work correctly.                      
C                                                                         
C     Mike Blackburn  15/05/2000.  Original version for IGCM2.            
C                                                                         
C     Purpose.                                                            
C     --------                                                            
C                                                                         
C        This subroutine modifies the m=0,n=0 coefficient of ln(ps)       
c     to correct the mean surface pressure (total mass) at that at the    
c     start of the run, preventing drift in extended integrations.        
C                                                                         
C     Interface.                                                          
C     ----------                                                          
C                                                                         
C     MASCOR is called from the main program MLTRI.                       
C     Timestep counters are stored in common BATS.                        
C     Mean surface pressures and switches are stored in common STATS.     
C     Spectral arrays SP,SPA are stored in common SPECTR.                 
C                                                                         
C     Method.                                                             
C     -------                                                             
C                                                                         
C        Rescale the global mean surface pressure to its initial value    
c     by adding the required constant to the m=0,n=0 ln(ps) coefficient.  
C     Use the time-lagged mean pressure to avoid numerical instability.   
C                                                                         
C        It is assumed there are no sources or sinks in the mass          
C     continuity equation (i.e. will not work if SPA changes in the diabat
C     part of the time step).                                             
C     The correction should be applied to the                             
C     mass of dry air if precip/evap mass variations were included.       
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
      COMMON/STATS/GMSP0,GMSPMI,LMASCOR,LMASOLD,LMASPRT                   
      LOGICAL LMASCOR,LMASOLD,LMASPRT                                     
C                                                                         
C-----------------------------------------------------------------------  
C                                                                         
C     ------------------------------------------------------------------  
C                                                                         
C     Extract the global mass.                                            
C                                                                         
      GMSP=1.0+REAL(SPA(1))/SQRT(2.0)                                     
C                                                                         
C     ------------------------------------------------------------------  
C                                                                         
C     Preserve the mass at (re)start as required.                         
C     Note that the timestep counter has already been incremented.        
C Pressure at start of initial run. Otherwise read in from restart file.  
C                                                                         
      IF ((KOUNT.EQ.KSTART+1).AND.GMSP0.EQ.0.0) THEN                      
         GMSP0=GMSP                                                       
         GMSPMI=GMSP                                                      
         WRITE(2,6900) GMSP                                               
      ENDIF                                                               
C                                                                         
C     ------------------------------------------------------------------  
C                                                                         
C     Compute correction.                                                 
C                                                                         
      ZSPRAT=GMSPMI/GMSP0                                                 
      ZLNRAT=LOG(ZSPRAT)                                                  
      SP(1)=SP(1)-SQRT(2.)*CMPLX(ZLNRAT,0.)                               
C                                                                         
C     ------------------------------------------------------------------  
C                                                                         
C     Print diagnostics, correcting for (KITS-1) offset.                  
C                                                                         
      IF (KITS.GT.0) THEN                                                 
         KTEMP=KOUNT+1-KITS                                               
      ELSE                                                                
         KTEMP=KOUNT                                                      
      ENDIF                                                               
      IF (LMASPRT.AND.                                                    
     :    (KTEMP.LE.(KSTART+ITSPD).OR.MOD(KTEMP,ITSPD).EQ.0)) THEN        
         ZSPDIF=GMSPMI-GMSP0                                              
         WRITE(2,6901) GMSPMI,GMSP0,ZSPDIF,ZLNRAT,REAL(SP(1))             
      ENDIF                                                               
C                                                                         
C     ------------------------------------------------------------------  
C                                                                         
C     Update the time-lagged mass.                                        
C                                                                         
      GMSPMI=GMSP                                                         
C                                                                         
C     ------------------------------------------------------------------  
C                                                                         
C     Formats.                                                            
C                                                                         
 6900 FORMAT(' GLOBAL MASS RESET TO CURRENT RESTART VALUE =',1PE15.8)     
 6901 FORMAT('  PMEAN=',1PE15.8,' PM0=',1PE15.8                           
     :,' DIFF=',1PE15.8,' LN(PS) CORR=',1PE15.8,' LN(PS)=',1PE15.8)       
C                                                                         
C     ------------------------------------------------------------------  
C                                                                         
      RETURN                                                              
      END                                                                 
