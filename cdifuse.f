C**********************************************************               
C             SUBROUTINE DIFUSE                                           
C**********************************************************               
      SUBROUTINE DIFUSE                                                   
C                                                                         
C     Calculates spectral tendencies from restoration (if included)       
C     and biharmonic diffusion.                                           
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
C
       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM(3), AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS(3), with_TiO_and_VO, opacity_method

       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS
      CHARACTER(len = 6) opacity_method
C
C
C ER Modif
      REAL TDAMP
      real rfcoeff(nlevrf)                                                
      logical ifirstrf                                                    
      save ifirstrf, rfcoeff                                              
      data ifirstrf /.true./                                              
C                                                                         
C Add newtonian cooling and surface Rayleigh friction.                    
C

                                                                         
      IF (LRESTIJ) THEN                                                   
C ER Modif, use Newtonian cooling with trad=0.1 days to set initial profiles,
C     unless PORB EQ 0!
         IF (((KOUNT/ITSPD).LE.1).AND.(PORB.EQ.0)) THEN
C           (only do this if 3D run)
            IF (.NOT.L1DZENITH) THEN
C              ddamp(l) is set as 1./(2.*pi*restim)
               TDAMP=1./(PI2*0.1)
               DO 41 IHEM=1,NHEM                                                   
                  DO 42 L=1,NL                                                     
                     I=(IHEM-1)*NWJ2+(L-1)*IGA                                     
                     DO 43 J=1,NWJ2                                                
                        TT(I+J)=TT(I+J)-TDAMP*(T(I+J)-TTRES(I+J))               
                        ZT(I+J)=ZT(I+J)-TFRC(L)*Z(I+J)                             
                        DT(I+J)=DT(I+J)-TFRC(L)*D(I+J)                             
 43                  CONTINUE                                                       
 42               CONTINUE                                                         
 41            CONTINUE                               
            ENDIF
         ELSE
            DO 21 IHEM=1,NHEM                                                   
               DO 22 L=1,NL                                                     
                  I=(IHEM-1)*NWJ2+(L-1)*IGA                                     
                  DO 23 J=1,NWJ2                                                
                     TT(I+J)=TT(I+J)-DDAMP(L)*(T(I+J)-TTRES(I+J))               
                     ZT(I+J)=ZT(I+J)-TFRC(L)*Z(I+J)                             
                     DT(I+J)=DT(I+J)-TFRC(L)*D(I+J)                             
 23               CONTINUE                                                       
 22            CONTINUE                                                         
 21         CONTINUE                                         
         ENDIF
C                                                                         
C       No friction on planetary vorticity                                
C                                                                         
      DO 24 L=1,NL                                                        
         I=(L-1)*IGA+1                                                    
         ZT(I)=ZT(I)+TFRC(L)*EZ                                           
 24   CONTINUE                                                            
      ELSE                                                                
      IF(DAMP.GT.0.0) THEN                                                
         I=0                                                              
         IR=0                                                             
         DO 800 IHEM=1,NHEM                                               
            DO 20 L=1,NL                                                  
               DO 10 J=1,IDM                                              
                  I=I+1                                                   
                  IR=IR+1                                                 
                  ZT(I)=ZT(I)-DAMP*(Z(I)-ZRES(IR))                        
                  DT(I)=DT(I)-DAMP*(D(I)-DRES(IR))                        
                  TT(I)=TT(I)-DAMP*(T(I)-TRES(IR))                        
   10          CONTINUE                                                   
               I=I+IGA-IDM                                                
               IR=IR+IGM-IDM                                              
   20       CONTINUE                                                      
            I=NWJ2                                                        
            IR=IDM                                                        
  800    CONTINUE                                                         
      ENDIF                                                               
      ENDIF                                                               
C                                                                         
C     Add in biharmonic diffusion if required                             
C                                                                         
      IF (AK(2).GT.0.0) THEN
         DO 820 IHEM=1,NHEM
            DO 821 L=1,NL                                                 
               JOFF=(L-1)*IGA+(IHEM-1)*NWJ2                               
               J=JOFF                                                     
               DO 822 MP=1,MFP,MOCT                                       
                  DO 823 IN=MP,NFP,MH                                     
                     J=J+1                                                
                     AKZ =AK(IN+2-IHEM)                                   
                     AKDT=AK(IN-1+IHEM)                                   
                     ZT(J)=ZT(J)-AKZ *Z(J)                                
                     DT(J)=DT(J)-AKDT*D(J)                                
                     TT(J)=TT(J)-AKDT*T(J)                                
823               CONTINUE                                                
822            CONTINUE                                                   
               DO 824 KK=1,NTRAC                                          
                  J=JOFF                                                  
                  DO 825 MP=1,MFP,MOCT                                    
                     DO 826 IN=MP,NFP,MH                                  
                        J=J+1                                             
                        AKDT=AK(IN-1+IHEM)                                
                        TRAT(J,KK)=TRAT(J,KK)-AKDT*TRA(J,KK)              
826                  CONTINUE                                             
825               CONTINUE                                                
824            CONTINUE                                                   
821         CONTINUE                                                      
820      CONTINUE                                                         
C                                                                         
C        No diffusion on EZ (planetary vorticity)                         
C                                                                         
         I=1                                                              
         DO 30 L=1,NL                                                     
            ZT(I)=ZT(I)+AK(2)*EZ                                          
            I=I+IGA                                                       
30       CONTINUE                                                         
      ENDIF                                                               
C                                                                         
c                                                                         
c 11-08-97 Very crude Rayleigh Friction in the top NLEVRF levels.         
c          NLEVRF is set in the parameter line.                           
c 29-06-98 Changed to use RFCOEFF instead of TFRC.  SMR                   
C 11-05-00 NB RFCOEFF AND TFRC ARE BOTH RAYLEIGH FRICTION                 
C          COEFFICIENT ARRAYS.                      SH                    
c                                                                         
      if (ifirstrf) then                                                  
        ifirstrf=.false.                                                  
        do l=1,nlevrf                                                     
! The 0.0625 is a tunable parameter. 0.1875 recommended for L26.          
!          rfcoeff(l)=((1.0-((l-1)*(1.0/float(nlevrf))))/pi2)*0.0625 
          RFCOEFF(l)=((1.0-((l-1)*(1.0/float(nlevrf))))/pi2)*RFCOEFF_IN 
! KM Modif to adjust top layer Rayleigh friction from input
        enddo                                                             
      endif                                                               
      do ihem=1,nhem                                                      
        do l=1,nlevrf                                                     
          i=(ihem-1)*nwj2+(l-1)*iga                                       
          do j=1,nwj2                                                     
            zt(i+j)=zt(i+j)-rfcoeff(l)*z(i+j)                             
            dt(i+j)=dt(i+j)-rfcoeff(l)*d(i+j)                             
          enddo                                                           
        enddo                                                             
      enddo                                                               
      do l=1,nlevrf                                                       
        i=(l-1)*iga+1                                                     
        zt(i)=zt(i)+rfcoeff(l)*ez                                         
      enddo                                                               
      RETURN                                                              
      END                                                                 
