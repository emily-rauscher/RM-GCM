C**********************************************************               
C             SUBROUTINE RADSW                                            
C**********************************************************               
      SUBROUTINE RADSW (ZTAVE,ZWV,PQOF,PCLFR,PQLWP,AMU0,                  
     :                  ZCARDI,ZSCT,PDP,PR,SWALB,FLS,GA,GASCON)           
C                                                                         
C inputs                                                                  
C ZTAVE - misd slab temps(K) all vertical -level 1 ground                 
C ZWV - water MMR                                                         
C PQOF -ozone mmr                                                         
C PCLFR - cloud fractions                                                 
C PQLWP - cloud Liquid water path  g/m^2                                  
C PR - flux pressures (mb)                                                
C -sWALB  surface albedo                                                  
C -AMU0    cosine                                                         
C
C ZSCT - solar constant * day length * earth-sun factor(Wm2)              
C                                                                         
C Outputs                                                                 
C FLS - NET SW fluxes level 1 ground                                      
C PDP - shortwave DPs                                                     
C                                                                         
C        ORIGINAL : 23-2-98                                               
C        Calls SW scheme and returns fluxes                               
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,ABSMMRLW,
     & JSKIPLON,JSKIPLAT,NSWMODEL,NLWMODEL,ABSSW1,ABSSTRAT,PRMIN,ALBSW1,
     & ABSSW2,SCATSW2,ASYMSW2,ABSLW1,NEWTB,NEWTE, with_TiO_and_VO, opacity_method
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR
       CHARACTER(len=6) :: opacity_method

      INTEGER NLON,NLEV,NLM,KDLON,JLON,JK,IMP,IABS,IAER,KFLEV             
     : ,JAE,NLP                                                           
      PARAMETER(NLON=1,NLEV=NL+1)                                         
      PARAMETER(NLM=NLEV-1,KDLON=NLON                                     
     :          ,KFLEV=NLEV,NLP=KFLEV+1)                                  
      REAL PRMU0(NLON),SWALB,FLS(NLP)                                     
      REAL PQLWP(NLON,NLEV),AMU0                                          
      REAL ZEPSCO,ZEPSC,ZEPSCQ,ZEPSCT,ZEPSCW,ZEELOG                       
      REAL CDAY,CCLWMR,DAYL                                               
      REAL ZALBSU(NLON,2),ZFSUP(NLON,NLP),ZFSDWN(NLON,NLP)                
     :,     ZPMB(NLON,NLP)                                                
      REAL ZTAVE(NLON,NLEV),ZOZ(NLON,NLEV),ZCLDSW(NLON,NLEV)              
     :,    ZTAU(NLON,2,NLEV),ZCG(NLON,2,NLEV),ZHEATR(NLON,NLEV)           
     :,    ZOMEGA(NLON,2,NLEV),PAER(NLON,NLEV,5)                          
c      REAL ZLWGKG,ZFCC                                                   
      REAL ZFLWP,PIERSLWP,ZRADEF                                          
     :,    ZTAUEQ                                                         
C      REAL ZFCCA,ZFCCB                                                   
C      REAL PPSOL(NLON),PR(NLON,NLEV),                                     
      REAL PR(NLON,NLEV),                                     
     :      PQOF(NLON,NLEV),PCLFR(NLON,NLEV),ZWV(NLON,NLEV)               
      REAL PDP(NLON,NLEV)                                                 
      REAL ZZRMUZ                                                         
      REAL ZSCT,ZCARDI                                                    

C---------- KM Modified Version. One needs PDPs for heating rate calculation

      write(*,*) 'In cradsw, shouldnt be here I think'
      write(*,*) 'I wrote this flag because I (Malsky) redefined some variables here'
      stop

      DO 113 JLON = 1 , KDLON                                             
         ZALBSU(JLON,1)=SWALB                                             
         ZALBSU(JLON,2)=SWALB                                             
         ZFSUP(JLON,KFLEV+1) = 0.                                         
         ZFSDWN(JLON,KFLEV+1) = 0.                                    
C         PPSOL(JLON)=100.*PR(JLON,1)   ! surface Pa                       
C         ZPMB(JLON,1)=PPSOL(JLON)/100.0  ! mb surface p                   
         ZPMB(JLON,1)=PR(JLON,1)
 113  CONTINUE                                                            
C        ZPMB(1,NLP)=0.0                                                   
C        ZPMB(1,2)=ZPMB(1,1)  ! bottom level has no thickness              
        DO JK=KFLEV,1,-1                                                  
         ZPMB(1,JK+1) = PR(1,JK)                                          
        ENDDO                                                             
        DO JK=2,KFLEV                                                     
         PDP(1,JK)=100.0*(ZPMB(1,JK)-ZPMB(1,JK+1))                        
        ENDDO                                                             
         PDP(1,1)=0.0                                                     

C      DO 503 JK = NLEV+1,1,-1                                             
CC            FLS(JK) = ZSCT* AMU0 * EXP(- 0.3 * PR(1,JK)/PR(1,1))       
C         FLS(JK) = ZSCT* AMU0 * EXP(- ABSSW1 * PR(1,JK)/PR(1,1))      
C 503  CONTINUE   

C NSWMODEL=1 is SW clear sky with surface reflection specified by albedo
C See p. 843 of Stephens (1984) RT review
         IF(NSWMODEL.EQ.1) THEN

C            PRMIN=0.15          ! Pressure level above which extra absorption applies
C            ABSSTRAT=0.2  ! Extra stratospheric absorption, in cm^2/g
C            ABSSTRAT=P0/GA*ABSSTRAT/10.0 ! Non-dimensionalize to match code (see inisimprad.f) 

            DO 602 JK = NLEV+1,1,-1                                             
               DO 601 JLON = 1 , KDLON                                          
                  IF (AMU0.GT.0) THEN
C                    Gradual absorption going down
                     ZFSDWN(JLON,JK)= ZSCT* AMU0 * EXP(- ABSSW1/AMU0 
     &                    * PR(1,JK)/PR(1,1))
                     IF (ABSSTRAT.NE.0) THEN ! Extra stratospheric absorption
                        IF ((PR(1,JK)/PR(1,1)).LT.PRMIN) THEN
                           ZFSDWN(JLON,JK)=ZFSDWN(JLON,JK)
     &                          *EXP(-ABSSTRAT/AMU0*PR(1,JK)/PR(1,1))
                        ELSE
                           ZFSDWN(JLON,JK)=ZFSDWN(JLON,JK)
     &                          *EXP(-ABSSTRAT/AMU0*PRMIN)
                        ENDIF
                     ENDIF
                     IF (ALBSW1.GT.0) THEN
                        WRITE(*,*) 'ERROR: non-zero surface albedo'
                        STOP
C     Non-absorbed flux at the bottom, times surface albedo, goes up and is
C     further absorbed, with an extra  1.66 boosting factor for diffusive regime 
                        ZFSUP(JLON,JK)= ALBSW1* ZSCT* AMU0 * 
     &                       EXP(-ABSSW1/AMU0-ABSSTRAT/AMU0)* 
     &                  EXP(- ABSSW1* 1.66 * (PR(1,1)-PR(1,JK))/PR(1,1))  
C     Extra stratospheric absorption above some pressure level
                        IF((PR(1,JK)/PR(1,1)).LT.PRMIN) THEN
                           ZFSUP(JLON,JK)=ZFSUP(JLON,JK) 
     &                  *EXP(-ABSSTRAT/PRMIN*(PRMIN-(PR(1,JK)/PR(1,1))))
                        ENDIF
                     ELSE
                        ZFSUP(JLON,JK)=0.0
                     ENDIF
                  ELSE
                     ZFSDWN(JLON,JK)=0.0
                     ZFSUP(JLON,JK)=0.0
                  ENDIF

                  FLS(JK) = ZFSDWN(JLON,JK) - ZFSUP(JLON,JK)                    
               
 601           CONTINUE                                                         
 602        CONTINUE     

         ENDIF

C NSWMODEL=2 is SW semi-infinite scattering atmosphere (e.g Schneider & Liu 08)
      IF(NSWMODEL.EQ.2) THEN

! Could be optimized by calculating fixed parameters once and storing results
         PAR1=SQRT(1.0-SCATSW2)
         PAR2=SQRT(1.0-SCATSW2*ASYMSW2)
         PARGAMMA=2.0*PAR1*PAR2
         PARRINF=(PAR2-PAR1)/(PAR2+PAR1)

      DO 611 JK = NLEV+1,1,-1                                             
         DO 610 JLON = 1 , KDLON                                          
            FLS(JK) =  ZSCT* AMU0/ 3.141592* (1.0-PARRINF) *
     & EXP(-PARGAMMA* ABSSW2 * PR(1,JK)/PR(1,1))                 
 610     CONTINUE                                                         
 611  CONTINUE                                                            


      ENDIF

C NSWMODEL=3 should become a modified version of Morcrette / Fouquart-Bonnel  
      IF(NSWMODEL.EQ.3) THEN

         WRITE(2,*) ' NEED IMPLEMENTATION OF MODIFIED MORCRETTE SW SCHEME'
         CALL ABORT
      ENDIF



      RETURN                                                              
      END                                                                 


