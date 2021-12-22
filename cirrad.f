C***********************************************************************  
C*                         SUBROUTINE IRRAD                            *  
C***********************************************************************  
      SUBROUTINE IRRAD(GAS,CKD,NLEV,PRES,PFLUX,TEMP,MMR,CL,ACL,LWCCL,     
     $                BASCL,TOPCL,FUP,FDWN,FNET,DP,LLBLM,GR,TLEV)              
C 3-1-98                                                                  
C Piers commented out DP calc - it is now passed                          
C 18-1-00 SH                                                              
C Fuller gas tables and choice between LBLM or NBM gas tables.            
C                                                                         
C Subroutine IRRAD calculates the the upward, downward and net            
C irradiances at each level of a given atmospheric profile, for a         
C given absorber                                                          
C                                                                         
C INPUT:  Selection of the absorbers (GAS)                                
C         Water vapour continuum index (CKD)                              
C         Number of levels in the atmosphere without the surface (NLEV)   
C         Pressure at each level (PRES)                                   
C         Pressure at flux levels (PFLUX)                                 
C         Temperature at each level (TEMP)                                
C         Absorber mass mixing ratio at each level (MMR)                  
C         Selection of cloud types (CL)                                   
C         Cloud amount as a fraction (ACL)                                
C         Liquid water content (LWCCL)                                    
C         Cloud base and top (BASCL,TOPCL)                                
C         Slab pressure DP                                                
C         Gas table switch LLBLM                                          
C OUTPUT: Upward, Downward and Net irradiances at each mid-level          
C         (FUP,FDWN,FNET)                                                 
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                         P A R R A Y                                 *  
C*                                                                     *  
C***********************************************************************  
      include 'params.i'


       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 
      
       LOGICAL LPLOTMAP

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,ABSMMRLW,
     & JSKIPLON,JSKIPLAT,NSWMODEL,NLWMODEL,ABSSW1,ABSSTRAT,PRMIN,ALBSW1,
     & ABSSW2,SCATSW2,ASYMSW2,ABSLW1,NEWTB,NEWTE, with_TiO_and_VO
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR


C-----------------------------------------------------------------------  
C In this part of the code the user can set the dimensions of arrays      
C used for the calculations, by setting the values in the paramter        
C statements                                                              
C-----------------------------------------------------------------------  
                                                                          
      INTEGER MXGAS,MXLEV,MXBAND,MXCL                                     
                                                                          
      PARAMETER(MXGAS=8)     ! Maximum number of gases                    
                                                                          
      PARAMETER(MXLEV=NL)     ! Maximum number of levels in the            
                             ! atmosphere (not including the surface)     
                                                                          
      PARAMETER(MXBAND=9)    ! Maximum number of spectral bands (not      
                             ! including the whole spectrum, 0-3000cm-1)  
                                                                          
      PARAMETER(MXCL=3)      ! Maximum number of cloud types              
                                                                          
C-----------------------------------------------------------------------  
C-----------------                                                        
C Input Variables                                                         
C-----------------                                                        
      INTEGER GAS(MXGAS),CKD                                              
      INTEGER NLEV                                                        
      REAL PRES(0:MXLEV),PFLUX(0:MXLEV)                                   
      REAL TEMP(0:MXLEV),MMR(MXGAS,0:MXLEV)                               
      INTEGER CL(0:MXCL)                                                  
      REAL ACL(0:MXCL),LWCCL(MXCL),BASCL(0:MXCL),TOPCL(0:MXCL)            
      REAL GR      ! Constants  (grav. acceleration)                      
C------------------                                                       
C Output Variables                                                        
C------------------                                                       
      REAL FUP(0:MXLEV),FDWN(0:MXLEV),FNET(0:MXLEV)                       
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER IGAS,IBND,IHAL,ILAY,ICOUNT,IUP,IDWN,ICL                     
      INTEGER CLOUD   ! Clear/Cloudy sky index                            
                                                                          
      REAL ULAY(MXGAS,MXLEV) ! Abs. amount between 2 successive levels    
      REAL H2O(MXLEV)     ! Water vapour mid-point mixing ratio           
      REAL DP(MXLEV)      ! Thickness of the atmospheric layers           
      REAL PMID(MXLEV)    ! Mid-point pressures                           
      REAL TMID(MXLEV)    ! Mid-point temperatures                        
      REAL QULAY(MXLEV)   ! Abs. amount in sub-layers                     
      REAL QH2O(MXLEV)    ! Water vapour mixing ratio in sub-layers       
      REAL QDP(MXLEV)     ! Thickness of the sub-layers                   
      REAL QP(MXLEV)      ! Pressure in sub-layers                        
      REAL QT(MXLEV)      ! Temperature in sub-layers                     
      REAL TLEV(0:MXLEV)  ! Temperature at mid-levels                     
      REAL TSURF          ! Surface temperature                           
                                                                          
      REAL UPATH       ! Absorber amount in an atmospheric path           
      REAL WVEFF       ! Effective water vapour abs. amount               
      REAL PEFF,TEFF   ! Effective pressure & effective temprature        
                                                                          
      REAL WVCNT(0:MXLEV-1,MXLEV) ! Effective H2O amount for              
                                    ! continuum absorption                
      REAL UH2O(0:MXLEV-1,MXLEV)  ! H2O amount                            
      REAL PH2O(0:MXLEV-1,MXLEV)  ! H2O effective pressure                
      REAL TH2O(0:MXLEV-1,MXLEV)  ! H2O effective pressure                
                                                                          
      REAL TRANS    ! Transmittance                                       
      REAL TRCNT    ! Continuum transmittance                             
      REAL HALTRAN  ! Halocarbon transmittance                            
                                                                          
      REAL SLP,INTCPT  ! Variables used for interpolation                 
                                                                          
      REAL TRIN(0:MXLEV-1,MXLEV)  ! Input/Output trans. used for          
      REAL TROUT(0:MXLEV-1,MXLEV) ! the treatment of spectral overlap     
                                                                          
      REAL PLANCK ! Planck function                                       
      REAL PF(MXLEV),DPF(0:MXLEV) ! Matrices used for the calculation     
                                    ! of the irradiances                  
      REAL TRCL(0:MXLEV-1,MXLEV)  ! Cloud transmittance                   
      REAL MTRCL(0:MXLEV-1,MXLEV) ! Modified tran. case                   
                                                                          
C Pre-computed tables for gas absorption                                  
                                                                          
      REAL AH2O(MXBAND,48,28),AGAS(MXBAND-1,48,28)                        
      REAL BH2O(MXBAND,48,28),BGAS(MXBAND-1,48,28)                        
      REAL CH2O(MXBAND,48,28),CGAS(MXBAND-1,48,28)                        
      REAL DH2O(MXBAND,48,28),DGAS(MXBAND-1,48,28)                        
      REAL AN03(48,28),BN03(48,28),CN03(48,28),DN03(48,28)                
      REAL AC02(48,28),BC02(48,28),CC02(48,28),DC02(48,28)                
      REAL AC06(48,28),BC06(48,28),CC06(48,28),DC06(48,28)                
      REAL AC07(48,28),BC07(48,28),CC07(48,28),DC07(48,28)                
                                                                          
      REAL FGUP(0:MXLEV)     ! Irradiances in bands of absorbers          
      REAL FGDWN(0:MXLEV)    ! (other than the water vapour               
                               !  0-3000 cm-1 band)                       
      REAL FUPCL(0:MXCL,MXLEV)  ! Irradiances at cloud levels in          
      REAL FDWNCL(0:MXCL,MXLEV) ! the 0-3000 cm-1 band                    
                                                                          
      REAL FBACK  ! Background upward irradiance                          
      LOGICAL LLBLM                                                       
      INTEGER IFIRST                                                      
      DATA IFIRST/0/  !! KM Modif  
C MTR COMMENT: By defining IFIRST to be zero in the above line, the 
C following IF statement conditional is skipped and no LLBLM gas files
C are read, I think                                                  
                                                                        
C     READ in tables                                                      
      IF (IFIRST.EQ.1) THEN                                               
        IFIRST=0                                                          
        IF(LLBLM)THEN                                                     
         READ(35) ah2o,bh2o,ch2o,dh2o,agas,bgas,cgas,dgas,ac02,bc02,      
     :        cc02,dc02,an03,bn03,cn03,dn03,ac06,bc06,cc06,dc06,          
     :        ac07,bc07,cc07,dc07                                         
         WRITE (2,*)'read in LBLM gas table'                              
        ELSE                                                              
          READ(36) ah2o,bh2o,ch2o,dh2o,agas,bgas,cgas,dgas,ac02,bc02,     
     :        cc02,dc02,an03,bn03,cn03,dn03,ac06,bc06,cc06,dc06,          
     :        ac07,bc07,cc07,dc07                                         
          WRITE (2,*)'read in NBM gas table'                              
        ENDIF                                                             
        WRITE(2,*)'gases on:',gas(1),gas(2),gas(3),gas(4),gas(5),         
     :                     gas(6),gas(7),gas(8)                           
       ENDIF                                                              
C Calculate pressure, temperature absorber amount and water vapour        
C mixing ratio between the first two mid-levels (PMID,TMID,ULAY,H2O),     
C split the layer into 2 and calculate the same for the upper sub-layer   
C (QLAY,QH2O,QP,QT). Sub-layers are required for calculations involving   
C modified transmittance. Also calculate the  temperature at              
C  mid-levels (TLEV)                                                      

       IF(LLOGPLEV) THEN
!        WRITE(2,*) 'NEED IMPLEMENTATION OF LOG-LEVELS SUB-LAYERS'
!        CALL ABORT
       ENDIF
      
!      write(*,*) 'line 181 cirrad.f'
!      write(*,*) pres                                                                          
      PMID(1)=0.5*(0.5*(PRES(1)+PRES(2))+PRES(0))                         
      QP(1)=0.5*(PMID(1)+0.5*(PRES(1)+PRES(2)))                           
 
!     ER Modif to prevent calculation from Tsurf ~500K
C      SLP=(TEMP(0)-TEMP(2))/(PRES(0)-PRES(2))                             
C      INTCPT=PRES(0)*TEMP(2)-PRES(2)*TEMP(0)                              
C      INTCPT=INTCPT/(PRES(0)-PRES(2))                                     
      SLP=(TEMP(1)-TEMP(2))/(PRES(1)-PRES(2))                             
      INTCPT=PRES(1)*TEMP(2)-PRES(2)*TEMP(1)                              
      INTCPT=INTCPT/(PRES(1)-PRES(2))                                     
      TMID(1)=SLP*PMID(1)+INTCPT                                          
      QT(1)=SLP*QP(1)+INTCPT                                              
                                                                          
C   Temperature at the surface and the first mid-level                    
      TLEV(0)=TEMP(0)                                                     
      TLEV(1)=0.5*SLP*(PRES(1)+PRES(2))+INTCPT                            
                                                                          
      DO IGAS=1,MXGAS                                                     
                                                                          
         IF (GAS(IGAS).EQ.1) THEN                                         
                                                                          
            IF ((IGAS.EQ.6).AND.(GAS(2).EQ.1)) THEN                       
               ULAY(IGAS,1)=ULAY(2,1)                                     
               GOTO 1972                                                  
            END IF                                                        
                                                                          
            IF ((IGAS.EQ.7).AND.(GAS(3).EQ.1)) THEN                       
               ULAY(IGAS,1)=ULAY(3,1)                                     
               GOTO 1972                                                  
            END IF                                                        
                                      
C     ER Modif, for slp and intcpt use levels 1 and 2 (instead of 0,2)
            SLP=(MMR(IGAS,1)-MMR(IGAS,2))/(PRES(1)-PRES(2))               
            INTCPT=PRES(1)*MMR(IGAS,2)-PRES(2)*MMR(IGAS,1)                
            INTCPT=INTCPT/(PRES(1)-PRES(2))                               
            ULAY(IGAS,1)=SLP*PMID(1)+INTCPT                               
            IF (IGAS.EQ.1) THEN                                           
               H2O(1)=ULAY(IGAS,1)                                        
            END IF                                                        
            ULAY(IGAS,1)=ULAY(IGAS,1)*DP(1)/GR                            
                                                                          
            IF (IGAS.EQ.1) THEN                                           
                                                                          
               QULAY(1)=SLP*QP(1)+INTCPT                                  
               QH2O(1)=QULAY(1)                                           
C     ER Modif: use PMID(1) instead of PRES(0), because DPF(1) is now that
               QDP(1)=0.25*(2*PMID(1)-PRES(1)-PRES(2))                    
               QULAY(1)=QULAY(1)*QDP(1)/GR                                
                                                                          
            END IF                                                        
                                                                          
 1972       CONTINUE                                                      
                                                                          
         END IF                                                           
      END DO                                                              
                                                                          
C Calculate pressures and temperatures and absorber amounts               
C between two successive mid-levels and in the sub-layers                 
                                                                          
      DO ILAY=2,NLEV-1                                                    
                                                                          
C      Mid-point & sub-layer pressures                                    
                                                                          
         PMID(ILAY)=PRES(ILAY)                                            
C     ER Modif/hack to weigh at deeper pressures:
C         QP(ILAY)=0.5*(PMID(ILAY)+0.5*(PRES(ILAY)+PRES(ILAY-1)))          
C     original:
         QP(ILAY)=0.5*(PMID(ILAY)+0.5*(PRES(ILAY)+PRES(ILAY+1)))          

                                                                          
C      Mid-point & sub-layer temperatures                                 
                                                                          
         TMID(ILAY)=TEMP(ILAY)                                            
                                                                          
         SLP=(TEMP(ILAY)-TEMP(ILAY+1))/(PRES(ILAY)-PRES(ILAY+1))          
         INTCPT=PRES(ILAY)*TEMP(ILAY+1)-PRES(ILAY+1)*TEMP(ILAY)           
         INTCPT=INTCPT/(PRES(ILAY)-PRES(ILAY+1))                          
         QT(ILAY)=SLP*QP(ILAY)+INTCPT                                     
                                                                          
C      Temperature at mid-levels                                          
         TLEV(ILAY)=0.5*SLP*(PRES(ILAY)+PRES(ILAY+1))+INTCPT              

!! KM Try different temperature weigth
!!         TLEV(ILAY)=0.5*(TEMP(ILAY)+TEMP(ILAY+1))
                                                                          
C      Absorber amounts                                                   
                                                                          
C         DP(ILAY)=0.5*(PRES(ILAY-1)-PRES(ILAY+1))                        
C MTRTEST  
                                                                                 
         DO IGAS=1,MXGAS                                                  
                                                                          
            IF (GAS(IGAS).EQ.1) THEN                                      
                                                                          
               IF ((IGAS.EQ.6).AND.(GAS(2).EQ.1)) THEN                    
                  ULAY(IGAS,ILAY)=ULAY(2,ILAY)                            
                  GOTO 1973                                               
               END IF                                                     
                                                                          
               IF ((IGAS.EQ.7).AND.(GAS(3).EQ.1)) THEN                    
                  ULAY(IGAS,ILAY)=ULAY(3,ILAY)                            
                  GOTO 1973                                               
               END IF                                                     
                                                                          
               ULAY(IGAS,ILAY)=MMR(IGAS,ILAY)*DP(ILAY)/GR                 
               IF (IGAS.EQ.1) THEN                                        
                  H2O(ILAY)=MMR(IGAS,ILAY)                                
               END IF                                                     
                                                                          
               IF (IGAS.EQ.1) THEN                                        
                                                                          
                  SLP=(MMR(IGAS,ILAY)-MMR(IGAS,ILAY+1))/                  
     $                (PRES(ILAY)-PRES(ILAY+1))                           
                  INTCPT=PRES(ILAY)*MMR(IGAS,ILAY+1)-                     
     $                   PRES(ILAY+1)*MMR(IGAS,ILAY)                      
                  INTCPT=INTCPT/(PRES(ILAY)-PRES(ILAY+1))                 
                  QULAY(ILAY)=SLP*QP(ILAY)+INTCPT                         
                  QH2O(ILAY)=QULAY(ILAY)                                  
                  QDP(ILAY)=0.25*(PRES(ILAY-1)-PRES(ILAY+1))              
                  QULAY(ILAY)=QULAY(ILAY)*QDP(ILAY)/GR                    
               END IF                                                     
                                                                          
 1973          CONTINUE                                                   
                                                                           
            END IF                                                        
                                                                          
         END DO                                                           
                                                                          
      END DO                                                              
                                                                          
C Calculate pressure, temperature and absorber amount for                 
C the upper atmospheric layer and its sub-layers                          
                                                                          
      PMID(NLEV)=0.5*(PRES(NLEV)+0.5*PRES(NLEV-1))                        
      QP(NLEV)=0.5*PMID(NLEV)                                             
                                                                          
      TMID(NLEV)=TEMP(NLEV)                                               
      SLP=(TEMP(NLEV-1)-TEMP(NLEV))/(PRES(NLEV-1)-PRES(NLEV))             
      INTCPT=PRES(NLEV-1)*TEMP(NLEV)-PRES(NLEV)*TEMP(NLEV-1)              
      INTCPT=INTCPT/(PRES(NLEV-1)-PRES(NLEV))                             
      QT(NLEV)=SLP*QP(NLEV)+INTCPT                                        
                                                                          
C   Temperature at the upper mid-level                                    
      TLEV(NLEV)=0.5*SLP*PRES(NLEV)+INTCPT                                
                                                                          
      DO IGAS=1,MXGAS                                                     
                                                                          
         IF (GAS(IGAS).EQ.1) THEN                                         
                                                                          
            IF ((IGAS.EQ.6).AND.(GAS(2).EQ.1)) THEN                       
               ULAY(IGAS,NLEV)=ULAY(2,NLEV)                               
               GOTO 1974                                                  
            END IF                                                        
                                                                          
            IF ((IGAS.EQ.7).AND.(GAS(3).EQ.1)) THEN                       
               ULAY(IGAS,NLEV)=ULAY(3,NLEV)                               
               GOTO 1974                                                  
            END IF                                                        
                                                                          
            ULAY(IGAS,NLEV)=MMR(IGAS,NLEV)*DP(NLEV)/GR                    
            IF (IGAS.EQ.1) THEN                                           
               H2O(NLEV)=MMR(IGAS,NLEV)                                   
            END IF                                                        
                                                                          
            IF (IGAS.EQ.1) THEN                                           
                                                                          
               SLP=(MMR(IGAS,NLEV-1)-MMR(IGAS,NLEV))/                     
     $             (PRES(NLEV-1)-PRES(NLEV))                              
               INTCPT=PRES(NLEV-1)*MMR(IGAS,NLEV)-                        
     $                PRES(NLEV)*MMR(IGAS,NLEV-1)                         
               INTCPT=INTCPT/(PRES(NLEV-1)-PRES(NLEV))                    
                                                                          
               QULAY(NLEV)=SLP*QP(NLEV)+INTCPT                            
               QH2O(NLEV)=QULAY(NLEV)                                     
               QDP(NLEV)=0.25*PRES(NLEV-1)                                
               QULAY(NLEV)=QULAY(NLEV)*QDP(NLEV)/GR                       
                                                                          
            END IF                                                        
 1974       CONTINUE                                                      
                                                                          
         END IF                                                           
                                                                          
      END DO                                                              
                                                                          
C Set cloud index                                                         
                                                                          
      CLOUD=0                                                             
      DO ICL=0,MXCL                                                       
         IF (CL(ICL).EQ.1) THEN                                           
!!            CLOUD=1                                                       
         END IF                                                           
      END DO                                                              
                                                                          
C If cloudy skies then call the cloud scheme                              
                                                                          
      IF (CLOUD.EQ.1) THEN                                                
         CALL CLDTRAN(NLEV,PFLUX,DP,CL,ACL,LWCCL,BASCL,TOPCL,             
     $                TRCL,MTRCL)                                         
      END IF                                                              
                                                                          
C---------------------------------                                        
C Calculation of  the irradiances                                         
C---------------------------------                                        
                                                                          
C--- Background upward irradiance                                         
                                                                          
      TSURF=TLEV(0)                                                       
      FBACK=PLANCK(0,TSURF)                                               
                                                                          
C--- Initialise irradiances                                               
                                                                          
      DO ILAY=0,NLEV                                                      
         FUP(ILAY)=FBACK                                                  
         FDWN(ILAY)=0.0                                                   
         FNET(ILAY)=FBACK                                                 
      END DO                                                              

                                                                          
C--- Calculate the irradiances for water vapour absorption.               
C    The wideband covers the whole thermal IR spectrum (0-3000cm-1)       
                                                                          
      IF (GAS(1).EQ.1) THEN                                               
                                                                          
         CALL H2OFLUX(CKD,NLEV,ULAY,H2O,DP,PMID,TMID,QULAY,QH2O,QDP,      
     $                QP,QT,TLEV,CLOUD,TRCL,MTRCL,WVCNT,UH2O,PH2O,        
     $                TH2O,FUP,FDWN)                                      
                                                                          
C--- Avoid smoothing at cloud levels                                      
                                                                          
         IF (CLOUD.EQ.1) THEN                                             
            DO ICL=0,MXCL                                                 
               IF (CL(ICL).EQ.1) THEN                                     
                  DO ILAY=1,NLEV                                          
                     IF ((100.*PFLUX(ILAY).LE.(BASCL(ICL)+1.D4)).AND.     
     $                   (100.*PFLUX(ILAY).GE.(TOPCL(ICL)-1.D4))) THEN    
                        FUPCL(ICL,ILAY)=FUP(ILAY)                         
                        FDWNCL(ICL,ILAY)=FDWN(ILAY)                       
                     ELSE                                                 
                        FUPCL(ICL,ILAY)=0.0                               
                        FDWNCL(ICL,ILAY)=0.0                              
                     END IF                                               
                  END DO                                                  
               END IF                                                     
            END DO                                                        
         END IF                                                           
                                                                          
C--- Smoothing of the irradiance profiles                                 
                                                                          
!         CALL SMFLUX(NLEV,PFLUX,FUP)                                      
!         CALL SMFLUX(NLEV,PFLUX,FDWN)                                     
                                                                          
         IF (CLOUD.EQ.1) THEN                                             
            DO ICL=0,MXCL                                                 
               IF (CL(ICL).EQ.1) THEN                                     
                  DO ILAY=1,NLEV                                          
                     IF (FUPCL(ICL,ILAY).NE.0.0) THEN                     
                        FUP(ILAY)=FUPCL(ICL,ILAY)                         
                        FDWN(ILAY)=FDWNCL(ICL,ILAY)                       
                     END IF                                               
                  END DO                                                  
               END IF                                                     
            END DO                                                        
         END IF                                                           
                                                                          
         FUP(0)=PLANCK(0,TLEV(0))                                         
         FDWN(NLEV)=0.0                                                   

         DO ILAY=0,NLEV                                                   
            FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                               
         END DO                                                           
                                                                          
      END IF                                                              

      DO IBND=1,MXBAND                                                    

      GOTO 1975   !! KM Modif to avoid all bands
                                                                          
C Skip the loop, if there is no gas that absorbs in the band              
                                                                          
         IF ((IBND.EQ.1).AND.(GAS(2).EQ.0)) THEN                          
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.2).AND.(GAS(3).EQ.0).AND.(GAS(6).EQ.0)) THEN        
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.3).AND.(GAS(4).EQ.0).AND.(GAS(5).EQ.0)) THEN        
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.4).AND.(GAS(5).EQ.0)) THEN                          
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.5).AND.(GAS(6).EQ.0)) THEN                          
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.6).AND.(GAS(7).EQ.0)) THEN                          
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.7).AND.(GAS(5).EQ.0)) THEN                          
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.8).AND.(GAS(5).EQ.0)) THEN                          
            GOTO 1975                                                     
         END IF                                                           
         IF ((IBND.EQ.9).AND.(GAS(8).EQ.0)) THEN                          
            GOTO 1975                                                     
         END IF                                                           
                                                                          
C Initialise input transmittances (used for overlap treatment)            
                                                                          
         DO IUP=1,NLEV                                                    
            DO IDWN=0,NLEV-1                                              
               IF (CLOUD.EQ.1) THEN                                       
                  TRIN(IDWN,IUP)=TRCL(IDWN,IUP) ! Cloudy skies            
               ELSE                                                       
                  TRIN(IDWN,IUP)=1.0            ! Clear skies             
               END IF                                                     
            END DO                                                        
         END DO                                                           
                                                                          
C Calculate arrays PF and DPF used in the calculation of the fluxes       
C PF:  Integrated Planck function at mid-points                           
C DPF: Differences of the PFs                                             
                                                                          
         DO ILAY=1,NLEV                                                   
            PF(ILAY)=PLANCK(IBND,TMID(ILAY))                              
         END DO                                                           
                                                                          
         DPF(0)=PLANCK(IBND,TSURF)-PLANCK(IBND,TMID(1))                   
         DPF(NLEV)=PLANCK(IBND,TMID(NLEV))                                
                                                                          
         DO ILAY=1,NLEV-1                                                 
            DPF(ILAY)=PLANCK(IBND,TMID(ILAY))-PLANCK(IBND,TMID(ILAY+1))   
         END DO                                                           
                                                                          
                                                                          
C--- If overlap with H2O then calculate H2O irradiances in the            
C    absorber's band and subtract them from the total                     
                                                                          
         IF (GAS(1).EQ.1) THEN                                            
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                                                                          
                  UPATH=UH2O(IDWN,IUP)                                    
                  PEFF=PH2O(IDWN,IUP)                                     
                  TEFF=TH2O(IDWN,IUP)                                     
                                                                          
                  CALL GASSEARCH(1,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,        
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                  TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                    
                                                                          
C          Calculate self-broadened continuum transmittance               
                  IF (CKD.EQ.1) THEN                                      
                     WVEFF=WVCNT(IDWN,IUP)                                
                     CALL GASCNT(IBND,WVEFF,TEFF,TRCNT)                   
                     TROUT(IDWN,IUP)=TROUT(IDWN,IUP)*TRCNT                
                  END IF                                                  
                                                                          
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
                                                                          
               END DO                                                     
            END DO                                                        
                                                                          
            IF ((IBND.EQ.6).AND.(GAS(2).EQ.1)) THEN                       
               CONTINUE                                                   
            ELSE                                                          
               IF ((IBND.EQ.7).AND.(GAS(2).EQ.1)) THEN                    
                  CONTINUE                                                
               ELSE                                                       
                  DO ILAY=0,NLEV                                          
                     FUP(ILAY)=FUP(ILAY)-FGUP(ILAY)                       
                     FDWN(ILAY)=FDWN(ILAY)-FGDWN(ILAY)                    
                     FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                      
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,NLEV-1                                           
                  TRIN(IDWN,IUP)=TROUT(IDWN,IUP)                          
               END DO                                                     
            END DO                                                        
                                                                          
         END IF                                                           
                                                                          
C Set halocarbon band index (for different overlap cases)                 
                                                                          
         IF (IBND.EQ.9) THEN                                              
            IHAL=9                                                        
            IF (GAS(2).EQ.1) THEN                                         
               IHAL=10                                                    
            END IF                                                        
            IF (GAS(3).EQ.1) THEN                                         
               IHAL=11                                                    
            END IF                                                        
            IF ((GAS(2).EQ.1).AND.GAS(3).EQ.1) THEN                       
               IHAL=12                                                    
            END IF                                                        
         END IF                                                           
                                                                          
C------------------------                                                 
C  BAND 1: 540-820 cm-1                                                   
C------------------------                                                 
                                                                          
         IF (IBND.EQ.1) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        

                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                  UPATH=0.0                                               
                  DO ILAY=IDWN+1,IUP                                      
                     UPATH=UPATH+ULAY(2,ILAY)                             
                  END DO                                                  
                                                                          
C          Calculate the effective pressure of the path                   
                  PEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     PEFF=PEFF+PMID(ILAY)*ULAY(2,ILAY)                    
                  END DO                                                  
                  PEFF=PEFF/UPATH                                         
                                                                          
C          Calculate the effective temperature of the path                
                  TEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     TEFF=TEFF+TMID(ILAY)*ULAY(2,ILAY)                    
                  END DO                                                  
                  TEFF=TEFF/UPATH                                         
                                                                          
                  CALL GASSEARCH(2,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,        
     $                     CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,       
     $                     CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,       
     $                     CC06,DC06,AC07,BC07,CC07,DC07,TRANS)           
                  TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                    
                                                                          
               END DO                                                     
            END DO                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
               END DO                                                     
            END DO                                                        
                                                                          
            IF (GAS(1).EQ.1) THEN                                         
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            ELSE                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FBACK-PLANCK(IBND,TSURF)+FGUP(ILAY)           
                  FDWN(ILAY)=FGDWN(ILAY)                                  
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
                                                                          
C------------------------                                                 
C  BAND 2: 980-1100 cm-1                                                  
C------------------------                                                 
                                                                          
         IF (IBND.EQ.2) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
            IF (GAS(3).EQ.1) THEN   ! O3 absorption                       
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(3,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(3,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(3,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(3,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,NLEV-1                                        
                     TRIN(IDWN,IUP)=TROUT(IDWN,IUP)                       
                  END DO                                                  
               END DO                                                     
                                                                          
            END IF                                                        
                                                                          
            IF (GAS(6).EQ.1) THEN   ! CO2 absorption                      
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(6,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(6,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(6,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(6,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
            END IF                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
               END DO                                                     
            END DO                                                        
                                                                          
            IF (GAS(1).EQ.1) THEN                                         
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            ELSE                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)       
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
                                                                          
C-------------------------                                                
C  BAND 3: 1200-1400 cm-1                                                 
C-------------------------                                                
                                                                          
         IF (IBND.EQ.3) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
            IF (GAS(4).EQ.1) THEN   ! CH4 absorption                      
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(4,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(4,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(4,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(4,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,NLEV-1                                        
                     TRIN(IDWN,IUP)=TROUT(IDWN,IUP)                       
                  END DO                                                  
               END DO                                                     
                                                                          
            END IF                                                        
                                                                          
            IF (GAS(5).EQ.1) THEN   ! N2O absorption                      
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(5,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(5,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(5,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(5,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
            END IF                                                        
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
               END DO                                                     
            END DO                                                        
                                                                          
            IF (GAS(1).EQ.1) THEN                                         
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            ELSE                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)       
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
                                                                          
C-------------------------                                                
C  BAND 4: 2100-2300 cm-1                                                 
C-------------------------                                                
                                                                          
         IF (IBND.EQ.4) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
                                    ! N2O absorption                      
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                  UPATH=0.0                                               
                  DO ILAY=IDWN+1,IUP                                      
                     UPATH=UPATH+ULAY(5,ILAY)                             
                  END DO                                                  
                                                                          
C          Calculate the effective pressure of the path                   
                  PEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     PEFF=PEFF+PMID(ILAY)*ULAY(5,ILAY)                    
                  END DO                                                  
                  PEFF=PEFF/UPATH                                         
                                                                          
C          Calculate the effective temperature of the path                
                  TEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     TEFF=TEFF+TMID(ILAY)*ULAY(5,ILAY)                    
                  END DO                                                  
                  TEFF=TEFF/UPATH                                         
                                                                          
                  CALL GASSEARCH(5,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,        
     $                     CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,       
     $                     CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,       
     $                     CC06,DC06,AC07,BC07,CC07,DC07,TRANS)           
                  TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                    
                                                                          
               END DO                                                     
            END DO                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
               END DO                                                     
            END DO                                                        
                                                                          
            IF (GAS(1).EQ.1) THEN                                         
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            ELSE                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)       
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
                                                                          
C-----------------------                                                  
C  BAND 5: 900-980 cm-1                                                   
C-----------------------                                                  
                                                                          
         IF (IBND.EQ.5) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        

                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                  UPATH=0.0                                               
                  DO ILAY=IDWN+1,IUP                                      
                     UPATH=UPATH+ULAY(6,ILAY)                             
                  END DO                                                  
                                                                          
C          Calculate the effective pressure of the path                   
                  PEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     PEFF=PEFF+PMID(ILAY)*ULAY(6,ILAY)                    
                  END DO                                                  
                  PEFF=PEFF/UPATH                                         
                                                                          
C          Calculate the effective temperature of the path                
                  TEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     TEFF=TEFF+TMID(ILAY)*ULAY(6,ILAY)                    
                  END DO                                                  
                  TEFF=TEFF/UPATH                                         
                                                                          
                  CALL GASSEARCH(6,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,        
     $                     CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,       
     $                     CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,       
     $                     CC06,DC06,AC07,BC07,CC07,DC07,TRANS)           
                  TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                    
                                                                          
               END DO                                                     
            END DO                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
               END DO                                                     
            END DO                                                        
                                                                          
            IF (GAS(1).EQ.1) THEN                                         
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            ELSE                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)       
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
                                                                          
C-----------------------                                                  
C  BAND 6: 650-770 cm-1                                                   
C-----------------------                                                  
                                                                          
         IF (IBND.EQ.6) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
            IF (GAS(2).NE.1) THEN                                         
                                    ! O3 absorption only                  
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(7,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(7,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(7,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(7,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                     FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)        
                     FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)     
                  END DO                                                  
               END DO                                                     
                                                                          
               IF (GAS(1).EQ.1) THEN                                      
                  DO ILAY=0,NLEV                                          
                     FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                       
                     FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                    
                     FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                      
                  END DO                                                  
               ELSE                                                       
                  DO ILAY=0,NLEV                                          
                     FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)    
                     FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                    
                     FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                      
                  END DO                                                  
               END IF                                                     
                                                                          
            ELSE                                                          
                                    ! CO2 + O3 absorption                 
                                    ! (a) CO2                             
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(2,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(2,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(2,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(2,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                     FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)        
                     FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)     
                  END DO                                                  
               END DO                                                     
                                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)-FGDWN(ILAY)                       
               END DO                                                     
                                                                          
                                    ! (b) O3                              
                                                                          
C Boundary Conditions                                                     
               FGUP(0)=PLANCK(IBND,TSURF)                                 
               FGDWN(NLEV)=0.0                                            
                                                                          
C Initialise the irradiances                                              
               DO ILAY=1,NLEV                                             
                  FGUP(ILAY)=PF(ILAY)                                     
                  FGDWN(ILAY-1)=PF(ILAY)                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
                     TRIN(IDWN,IUP)=TROUT(IDWN,IUP)                       
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(7,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(7,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(7,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(7,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                     FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)        
                     FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)     
                  END DO                                                  
               END DO                                                     
                                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
                                                                          
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
                                                                          
C-----------------------                                                  
C  BAND 7: 550-630 cm-1                                                   
C-----------------------                                                  
                                                                          
         IF (IBND.EQ.7) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
            IF (GAS(2).NE.1) THEN                                         
                                    ! N2O absorption only                 
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(5,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(5,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(5,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(5,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                     FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)        
                     FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)     
                  END DO                                                  
               END DO                                                     
                                                                          
               IF (GAS(1).EQ.1) THEN                                      
                  DO ILAY=0,NLEV                                          
                     FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                       
                     FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                    
                     FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                      
                  END DO                                                  
               ELSE                                                       
                  DO ILAY=0,NLEV                                          
                     FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)    
                     FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                    
                     FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                      
                  END DO                                                  
               END IF                                                     
                                                                          
            ELSE                                                          
                                    ! CO2 + N2O absorption                
                                    ! (a) CO2                             
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(2,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(2,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(2,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(2,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,     
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                     FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)        
                     FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)     
                  END DO                                                  
               END DO                                                     
                                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)-FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
                                                                          
                                    ! (b) N2O                             
                                                                          
C Boundary Conditions                                                     
               FGUP(0)=PLANCK(IBND,TSURF)                                 
               FGDWN(NLEV)=0.0                                            
                                                                          
C Initialise the irradiances                                              
               DO ILAY=1,NLEV                                             
                  FGUP(ILAY)=PF(ILAY)                                     
                  FGDWN(ILAY-1)=PF(ILAY)                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                                                                          
                     TRIN(IDWN,IUP)=TROUT(IDWN,IUP)                       
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                     UPATH=0.0                                            
                     DO ILAY=IDWN+1,IUP                                   
                        UPATH=UPATH+ULAY(5,ILAY)                          
                     END DO                                               
                                                                          
C          Calculate the effective pressure of the path                   
                     PEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        PEFF=PEFF+PMID(ILAY)*ULAY(5,ILAY)                 
                     END DO                                               
                     PEFF=PEFF/UPATH                                      
                                                                          
C          Calculate the effective temperature of the path                
                     TEFF=0.0                                             
                     DO ILAY=IDWN+1,IUP                                   
                        TEFF=TEFF+TMID(ILAY)*ULAY(5,ILAY)                 
                     END DO                                               
                     TEFF=TEFF/UPATH                                      
                                                                          
                     CALL GASSEARCH(5,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,  
     $                        CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,    
     $                        CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,    
     $                        CC06,DC06,AC07,BC07,CC07,DC07,TRANS)        
                     TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                 
                                                                          
                  END DO                                                  
               END DO                                                     
                                                                          
               DO IUP=1,NLEV                                              
                  DO IDWN=0,IUP-1                                         
                     FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)        
                     FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)     
                  END DO                                                  
               END DO                                                     
                                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
                                                                          
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
C-------------------------                                                
C  BAND 8: 1100-1200 cm-1                                                 
C-------------------------                                                
                                                                          
         IF (IBND.EQ.8) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
                                    ! N2O absorption                      
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                  UPATH=0.0                                               
                  DO ILAY=IDWN+1,IUP                                      
                     UPATH=UPATH+ULAY(5,ILAY)                             
                  END DO                                                  
                                                                          
C          Calculate the effective pressure of the path                   
                  PEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     PEFF=PEFF+PMID(ILAY)*ULAY(5,ILAY)                    
                  END DO                                                  
                  PEFF=PEFF/UPATH                                         
                                                                          
C          Calculate the effective temperature of the path                
                  TEFF=0.0                                                
                  DO ILAY=IDWN+1,IUP                                      
                     TEFF=TEFF+TMID(ILAY)*ULAY(5,ILAY)                    
                  END DO                                                  
                  TEFF=TEFF/UPATH                                         
                                                                          
                  CALL GASSEARCH(5,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,        
     $                     CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,       
     $                     CN03,DN03,AC02,BC02,CC02,DC02,AC06,BC06,       
     $                     CC06,DC06,AC07,BC07,CC07,DC07,TRANS)           
                  TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                    
                                                                          
               END DO                                                     
            END DO                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
               END DO                                                     
            END DO                                                        
                                                                          
            IF (GAS(1).EQ.1) THEN                                         
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            ELSE                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)       
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            END IF                                                        
                                                                          
            GOTO 1975                                                     
                                                                          
         END IF                                                           
                                                                          
C-----------------------                                                  
C  BAND 9: 820-900 cm-1                                                   
C-----------------------                                                  
                                                                          
         IF (IBND.EQ.9) THEN                                              
                                                                          
C Boundary Conditions                                                     
            FGUP(0)=PLANCK(IBND,TSURF)                                    
            FGDWN(NLEV)=0.0                                               
                                                                          
C Initialise the irradiances                                              
            DO ILAY=1,NLEV                                                
               FGUP(ILAY)=PF(ILAY)                                        
               FGDWN(ILAY-1)=PF(ILAY)                                     
            END DO                                                        
                                                                          
                                    ! Halocarbon absorption               
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
                  UPATH=0.0                                               
                  DO ILAY=IDWN+1,IUP                                      
                     UPATH=UPATH+ULAY(8,ILAY)                             
                  END DO                                                  
                                                                          
                  TRANS=HALTRAN(IHAL,UPATH)                               
                  TROUT(IDWN,IUP)=TRANS*TRIN(IDWN,IUP)                    
                                                                          
               END DO                                                     
            END DO                                                        
                                                                          
            DO IUP=1,NLEV                                                 
               DO IDWN=0,IUP-1                                            
                  FGUP(IUP)=FGUP(IUP)+TROUT(IDWN,IUP)*DPF(IDWN)           
                  FGDWN(IDWN)=FGDWN(IDWN)-TROUT(IDWN,IUP)*DPF(IUP)        
               END DO                                                     
            END DO                                                        
                                                                          
            IF (GAS(1).EQ.1) THEN                                         
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)+FGUP(ILAY)                          
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            ELSE                                                          
               DO ILAY=0,NLEV                                             
                  FUP(ILAY)=FUP(ILAY)-PLANCK(IBND,TSURF)+FGUP(ILAY)       
                  FDWN(ILAY)=FDWN(ILAY)+FGDWN(ILAY)                       
                  FNET(ILAY)=FUP(ILAY)-FDWN(ILAY)                         
               END DO                                                     
            END IF                                                        
                                                                          
         END IF                                                           
                                                                          
 1975    CONTINUE                                                         
                                                                          
      END DO                                                              
                                                                          
                                                                          
      END                                                                 
C***********************************************************************  
C*                                                                     *  
C*                     FUNCTION PLANCK                                 *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      REAL FUNCTION PLANCK (IBND,TEMP)                                    
                                                                          
C This function returns the Planck function, integrated over a wide       
C band at specific temperature                                            
C                                                                         
C INPUT:  Band index (IBND)                                               
C         Temperature (TEMP)                                              
C OUTPUT: Planck Function, integrated over band IBND at                   
C         temperature TEMP (PLANCK)                                       
C The integrated Planck function is tabulated for each band in the        
C array B in plfunc, over a range of temperatures from 151 to 350 K,      
C with a step of 1 K. A simple linear interpolation in temperature is     
C employed.                                                               
C--------------------------------------------------                       
C WARNING: The function works only for temperatures                       
C          between 151 and 350 K                                          
C--------------------------------------------------                       
                                                                          
      IMPLICIT NONE                                                       
C***********************************************************************  
C*                                                                     *  
C*                         P A R R A Y                                 *  
C*                                                                     *  
C***********************************************************************  
C-----------------------------------------------------------------------  
C In this part of the code the user can set the dimensions of arrays      
C used for the calculations, by setting the values in the paramter        
C statements                                                              
C-----------------------------------------------------------------------  

      INTEGER NN,MM,NHEM,NL,MOCT,MG,JG,NWJ2,NCRAY,JGL,NTRAC,NLEVRF
      include 'params.i'
                           
      INTEGER MXGAS,MXLEV,MXBAND,MXCL                                     
                                                                          
      PARAMETER(MXGAS=8)     ! Maximum number of gases                    
                                                                          
      PARAMETER(MXLEV=NL)     ! Maximum number of levels in the            
                             ! atmosphere (not including the surface)     
                                                                          
      PARAMETER(MXBAND=9)    ! Maximum number of spectral bands (not      
                             ! including the whole spectrum, 0-3000cm-1)  
                                                                          
      PARAMETER(MXCL=3)      ! Maximum number of cloud types              
                                                                          
C-----------------------------------------------------------------------  
                                                                          
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER IBND,ICOUNT                                                 
      REAL TEMP,TEMP1,TEMP2                                               
                                                                          
                                                                          
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER IPF            ! Table element index                        
      REAL SLP,INTCPT      ! Variables used for interpolation             
      REAL B(0:MXBAND,200) ! Planck function pre-computed values          
                                                                          
      INTEGER IFIRST                                                      
      DATA IFIRST/0/                                                      
                                                                          
C     READ in planck function                                             
      IF (IFIRST.EQ.1) THEN                                               
         IFIRST=0                                                         
         read(37) B                                                       
         print *,'read in Planck function'                                
       ENDIF                                                              
                                                                          
C-------------------------                                                
C Band 0:    0-3000 cm-1                                                  
C Band 1:   540-820 cm-1                                                  
C Band 2:  980-1100 cm-1                                                  
C Band 3: 1200-1400 cm-1                                                  
C Band 4: 2100-2300 cm-1                                                  
C Band 5:   900-980 cm-1                                                  
C Band 6:   650-770 cm-1                                                  
C Band 7:   550-630 cm-1                                                  
C Band 8: 1100-1200 cm-1                                                  
C Band 9:   820-900 cm-1                                                  
C-------------------------                                                
                                                                          
C Catch temperatures that are out of bounds                               
!      IF(TEMP.LT.151.0.OR.TEMP.GT.350.0)THEN                              
!         print *,'WARNING: TEMP out of bounds for Planck function'        
!         print *,'model will crash. TEMP= ',TEMP                          
!      ENDIF                                                               
                                                                          
      IPF=INT(TEMP-150.)                                                  
      TEMP1=AINT(TEMP)                                                    
      TEMP2=TEMP+1.                                                       
                                                                          
      IF ((TEMP-TEMP1).LE.0.01) THEN  ! No need to interpolate            
         PLANCK=B(IBND,IPF)                                               
      ELSE                                 ! Linear interpolation         
         SLP=(B(IBND,IPF)-B(IBND,IPF+1))/(TEMP1-TEMP2)                    
         INTCPT=TEMP1*B(IBND,IPF+1)-TEMP2*B(IBND,IPF)                     
         INTCPT=INTCPT/(TEMP1-TEMP2)                                      
         PLANCK=SLP*TEMP+INTCPT                                           
      END IF                                                              
               

      PLANCK= 5.6704e-8*(TEMP**4) !No division by PI: integrated over angles 
C      PLANCK= 1.8050e-8*(TEMP**4) ! division by PI

C     ER Modif to correct for blackbody curve lost from wavelength region
c      IF (TEMP.GT.1000.) THEN
c         PLANCK=PLANCK*(1.28-TEMP*2.8E-4)
c      END IF
                                                           
      END                                                                 
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                         SUBROUTINE H2OFLUX                          *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      SUBROUTINE H2OFLUX(CKD,NLEV,ULAY,H2O,DP,PMID,TMID,QULAY,QH2O,QDP,   
     $                   QP,QT,TLEV,CLOUD,TRCL,MTRCL,WVCNT,UH2O,PH2O,     
     $                   TH2O,FUP,FDWN)                                   
                                                                          
C Subroutine H2OFLUX calculates the irradiances (upward, downward and     
C net) for water vapour line absorption. Wideband transmittance is used   
C for the upward, and modified wideband transmittance is used for the     
C downward irradiances, which cover the whole thermal IR                  
C spectrum (0-3000)                                                       
C                                                                         
C INPUT:  Water vapour continuum index (CKD)                              
C         Number of levels in the atmosphere without the surface (NLEV)   
C         Abs. amount between 2 successive levels (ULAY)                  
C         Water vapour mid-point mixing ratio (H2O)                       
C         Thickness of the atmospheric layers (DP)                        
C         Mid-point pressures (PMID)                                      
C         Mid-point temperatures(TMID)                                    
C         Abs. amount in sub-layers (QULAY)                               
C         Water vapour mixing ratio in sub-layers (QH2O)                  
C         Thickness of the sub-layers (QDP)                               
C         Pressure in sub-layers (QP)                                     
C         Temperature in sub-layers (QT)                                  
C         Temperature at mid-levels (TLEV)                                
C         Cloud index (CLOUD)                                             
C         Cloud transmittance (TRCL)                                      
C         Cloud transmittance for the modified transm. case (MTRCL)       
C                                                                         
C OUTPUT: Effective water vapour amount (WVCNT)                           
C         Water vapour - path parameters (UH2O,PH2O,TH2O)                 
C         Upward and Downward irradiances at each mid-level(FUP,FDWN)     
                                                                          
      IMPLICIT NONE                                                       
C***********************************************************************  
C*                                                                     *  
C*                         P A R R A Y                                 *  
C*                                                                     *  
C***********************************************************************  
C-----------------------------------------------------------------------  
C In this part of the code the user can set the dimensions of arrays      
C used for the calculations, by setting the values in the paramter        
C statements                                                              
C-----------------------------------------------------------------------  
                      
      INTEGER NN,MM,NHEM,NL,MOCT,MG,JG,NWJ2,NCRAY,JGL,NTRAC,NLEVRF
      include 'params.i'

      LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR
      REAL ABSMMRLW,ABSSW1,ALBSW1,ABSSTRAT,PRMIN,ABSSW2,SCATSW2,ASYMSW2,
     &     ABSLW1, with_TiO_and_VO
      INTEGER JSKIPLON,JSKIPLAT,NSWMODEL,NLWMODEL,NEWTB,NEWTE

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,ABSMMRLW,
     & JSKIPLON,JSKIPLAT,NSWMODEL,NLWMODEL,ABSSW1,ABSSTRAT,PRMIN,ALBSW1,
     & ABSSW2,SCATSW2,ASYMSW2,ABSLW1,NEWTB,NEWTE, with_TiO_and_VO

       LOGICAL LPLOTMAP
       REAL OOM_IN, RFCOEFF_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 
       INTEGER NLPLOTMAP_IN, NTSTEP_IN, NSKIP_IN

       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 
      

                                                    
      INTEGER MXGAS,MXLEV,MXBAND,MXCL                                     
                                                                          
      PARAMETER(MXGAS=8)     ! Maximum number of gases                    
                                                                          
      PARAMETER(MXLEV=NL)     ! Maximum number of levels in the            
                             ! atmosphere (not including the surface)     
                                                                          
      PARAMETER(MXBAND=9)    ! Maximum number of spectral bands (not      
                             ! including the whole spectrum, 0-3000cm-1)  
                                                                          
      PARAMETER(MXCL=3)      ! Maximum number of cloud types              
                                                                          
C-----------------------------------------------------------------------  
                                                                          
C-----------------                                                        
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER CKD,NLEV,CLOUD                                              
      REAL ULAY(MXGAS,MXLEV)                                              
      REAL H2O(MXLEV)                                                     
      REAL DP(MXLEV)                                                      
      REAL PMID(MXLEV)                                                    
      REAL TMID(MXLEV)                                                    
      REAL QULAY(MXLEV)                                                   
      REAL QH2O(MXLEV)                                                    
      REAL QDP(MXLEV)                                                     
      REAL QP(MXLEV)                                                      
      REAL QT(MXLEV)                                                      
      REAL TLEV(0:MXLEV)                                                  
      REAL TRCL(0:MXLEV-1,MXLEV)                                          
      REAL MTRCL(0:MXLEV-1,MXLEV)                                         
                                                                          
C------------------                                                       
C Output Variables                                                        
C------------------                                                       
                                                                          
      REAL WVCNT(0:MXLEV-1,MXLEV)                                         
      REAL UH2O(0:MXLEV-1,MXLEV)                                          
      REAL PH2O(0:MXLEV-1,MXLEV)                                          
      REAL TH2O(0:MXLEV-1,MXLEV)                                          
      REAL FUP(0:MXLEV),FDWN(0:MXLEV)                                     
                                                                          
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER ICOUNT,ILAY,IUP,IDWN                                        
      INTEGER MOD   ! Index for transmittance type                        
                    ! (1 : Modified, 2 : Standard)                        
                                                                          
      REAL PF(0:MXLEV),DPF(MXLEV) ! Matrices used for the calculation     
                                    ! of the irradiances                  
                                                                          
      REAL UPATH       ! Absorber amount in an atmospheric path           
      REAL WVEFF       ! Effective water vapour abs. amount               
      REAL PEFF,TEFF   ! Effective pressure & effective temprature        
      REAL CWVC        ! Parameter 1/(296K)                               
      PARAMETER(CWVC=3.3783783D-03)                                       
                                                                          
      REAL PLANCK      ! Planck function integrated over a band           
                                                                          
      REAL TRANS       ! Transmittance of a layer                         
                                                                          
C Pre-computed tables for H2O-line absorption                             
      REAL AWVL(2,34,28),BWVL(2,34,28)                                    
      REAL CWVL(2,34,28),DWVL(2,34,28)                                    
                                                                          
C Pre-computed tables for H2O-continuum absorption                        
      REAL AWV(2,41,29,22),BWV(2,41,29,22)                                
      REAL CWV(2,41,29,22),DWV(2,41,29,22)                                
                                                                          
      INTEGER IFIRST                                                      
      DATA IFIRST/0/  !!KM Modif                                                      
C KM Modif: absorption coeff integrated along UPATH
      REAL ABSCOEFF, CABSLW1
         
C Scale ABSLW1 from cm^2/g to the Pascal and MKS units used in cirrad.f
C ULAY in Pascal/GA (m /s^2) is 10 dyne/cm^2 divided per 100 cm/s^2 hence
C it is 0.1 ULAY in CGS (cm^2 /g). Thus, divide ABSLW1 / 10
       CABSLW1=ABSLW1/10.0
                                                      
C     READ in tables                                                      
      IF (IFIRST.EQ.1) THEN                                               
         IFIRST=0                                                         
         read(33) AWVL,BWVL,CWVL,DWVL,AWV,BWV,CWV,DWV                     
         print *,'read in water'                                          
       ENDIF                                                              
                               
                                           
C-----------------------------------------------------------------------  
C Calculate irradiances                                                   
C-----------------------------------------------------------------------  
                                                                          
C Boundary Conditions                                                     
                                                                          
      FUP(0)=PLANCK(0,TLEV(0))                                            
      FDWN(NLEV)=0.0                                                      
                                                                          
C-----------------------------------------                                
C Calculations for the upward irradiances                                 
C-----------------------------------------                                
                                                                          
      MOD=1                                                               
                                                                          
C Calculate arrays PF and DPF used in the calculation of the up. fluxes   
C PF:  Integrated Planck function at atmospheric levels                   
C DPF: Differences of the PFs                                             
                                                                          
      PF(0)=PLANCK(0,TLEV(0))                                             
                                                                          
      DO ILAY=1,NLEV                                                      
         PF(ILAY)=PLANCK(0,TLEV(ILAY))                                    
      END DO                                                              

C     ER Modif: use TMID(1) instead of TLEV(0) for DPF(1)
      DPF(1)=PLANCK(0,TLEV(1))-PLANCK(0,TMID(1))
C      DO ILAY=1,NLEV                                                      
      DO ILAY=2,NLEV                                                      
         DPF(ILAY)=PLANCK(0,TLEV(ILAY))-PLANCK(0,TLEV(ILAY-1))            
      END DO                                                              
                                                                        
C Initialise the irradiances                                              
                                                                          
      DO ILAY=1,NLEV                                                      
         FUP(ILAY)=PF(ILAY)                                               
      END DO                                                              
                                                                          
C For each layer, calculate the transmittance and update the irradiances  

      DO IUP=1,NLEV                                                       
                                                                          
         DO IDWN=0,IUP-1                                                  
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
            UPATH=QULAY(IDWN+1)  !*0.0                                         
            IF ((IUP-IDWN).GT.1) THEN                                     
               DO ILAY=IDWN+2,IUP                                         
                  UPATH=UPATH+ULAY(1,ILAY)                                
               END DO                                                     
            END IF 

C KM Modif: Integration Absorption along UPATH. Diffusive factor 1.66 applied
            IF(NLWMODEL.EQ.1) THEN
               
                IF (OPACIR_POWERLAW.eq.0) THEN
                   ABSCOEFF=CABSLW1*QULAY(IDWN+1)
                ELSE IF (OPACIR_POWERLAW.eq.1) THEN
                   ABSCOEFF=CABSLW1*QULAY(IDWN+1)
     &             *MAX(1e-6,(PMID(IDWN+1)/OPACIR_REFPRES))
                ELSE IF (OPACIR_POWERLAW.eq.2) THEN
                   ABSCOEFF=CABSLW1*QULAY(IDWN+1)     
     &           *MAX(1e-6,(PMID(IDWN+1)/OPACIR_REFPRES))
     &           *MAX(1e-6,(PMID(IDWN+1)/OPACIR_REFPRES))
                ELSE IF (OPACIR_POWERLAW.eq.3) THEN
                   ABSCOEFF=CABSLW1*QULAY(IDWN+1)     
     &           *MAX(1e-6,(PMID(IDWN+1)/OPACIR_REFPRES))
     &           *MAX(1e-6,(PMID(IDWN+1)/OPACIR_REFPRES))
     &           *MAX(1e-6,(PMID(IDWN+1)/OPACIR_REFPRES))
                ELSE  
                   ABSCOEFF=CABSLW1*QULAY(IDWN+1) 
     &           *MAX(1e-6,(PMID(IDWN+1)/OPACIR_REFPRES)**OPACIR_POWERLAW)

               ENDIF        
!KM Modif for Saturated limit (Hourdin): use total UPATH below, outside loop
               IF ((IUP-IDWN).GT.1) THEN                                     
                  DO ILAY=IDWN+2,IUP                                         
                     IF (OPACIR_POWERLAW.eq.0) THEN  !MTR Modif
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)
                     ELSE IF (OPACIR_POWERLAW.eq.1) THEN
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)  
     &           *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES)) 
                     ELSE IF (OPACIR_POWERLAW.eq.2) THEN
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)
     &           *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
     &           *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
                     ELSE IF (OPACIR_POWERLAW.eq.3) THEN
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)
     &           *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
     &           *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
     &           *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
                     ELSE 
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)
     &           *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES)**OPACIR_POWERLAW)
 
                     ENDIF

             
!KM Modif for Saturated limit (Hourdin): use total UPATH below, outside loop
                  END DO                                        
               END IF 
!KM Modif for Saturated limit (Hourdin): outside loop since no integral done
!            ABSCOEFF=CABSLW1*SQRT(UPATH*(PMID(IUP)+PMID(IDWN+1))/2e5) 
!!C           IF(IDWN.EQ.0) ABSCOEFF=CABSLW1*SQRT(UPATH*(PMID(IUP)+1e5)/2e5) 
            ENDIF

            IF(NLWMODEL.EQ.2) THEN
               WRITE(2,*) 'NEED IMPLEMENTATION OF PLANCK MEAN OPACITIES'
               CALL ABORT
               ABSCOEFF=CABSLW1*QULAY(IDWN+1)           
               IF ((IUP-IDWN).GT.1) THEN                                     
                  DO ILAY=IDWN+2,IUP                                         
                     ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)   
                  END DO                                                     
               END IF 
            ENDIF

C          Calculate the effective pressure of the path                   
C            PEFF=QP(IDWN+1)*QULAY(IDWN+1)                                 
C            IF ((IUP-IDWN).GT.1) THEN                                     
C               DO ILAY=IDWN+2,IUP                                         
C                  PEFF=PEFF+PMID(ILAY)*ULAY(1,ILAY)                       
C               END DO                                                     
C            END IF                                                        
C            PEFF=PEFF/UPATH                                               
                                                                          
C          Calculate the effective temperature of the path                
C            TEFF=QT(IDWN+1)*QULAY(IDWN+1)                                 
C            IF ((IUP-IDWN).GT.1) THEN                                     
C               DO ILAY=IDWN+2,IUP                                         
C                  TEFF=TEFF+TMID(ILAY)*ULAY(1,ILAY)                       
C               END DO                                                     
C            END IF                                                        
C            TEFF=TEFF/UPATH                                               
                                                                          
            IF (CKD.NE.1) THEN                                            
C--- Line absorption                                                      
!!               CALL LINSEARCH(MOD,UPATH,PEFF,TEFF,AWVL,BWVL,              
!!     $                        CWVL,DWVL,TRANS)                            
            ELSE                                                          
                                                                          
C--- Include continuum, if required                                       
                                                                          
C          Calculate effective water vapour amount                        
C               IF (QH2O(IDWN+1).GT.(1.0D-04)) THEN                        
C                  WVEFF=(QH2O(IDWN+1)**2)*PMID(IDWN+1)*QDP(IDWN+1)*       
C     $            EXP(1800.0*((1.0/TMID(IDWN+1))-CWVC))                   
C               ELSE                                                       
C                  WVEFF=0.0                                               
C               END IF                                                     
C               IF ((IUP-IDWN).GT.1) THEN                                  
C                  DO ILAY=IDWN+2,IUP                                      
C                     IF (H2O(ILAY).GT.(1.0D-04)) THEN                     
C                        WVEFF=WVEFF+(H2O(ILAY)**2)*PMID(ILAY)*DP(ILAY)*   
C     $                  EXP(1800.0*((1.0/TMID(ILAY))-CWVC))               
C                     END IF                                               
C                  END DO                                                  
C               END IF                                                     
C               WVEFF=(1.6178235D-06)*WVEFF                                
C               WVCNT(IDWN,IUP)=WVEFF                                      
                                                                          
               IF ((WVEFF.NE.0.0).AND.(LOG10(WVEFF).GT.(-4.5))) THEN      
!!                  CALL CONSEARCH(MOD,UPATH,PEFF,TEFF,WVEFF,AWV,BWV,       
!!     $                           CWV,DWV,TRANS)                           
               ELSE                                                       
!!                  CALL LINSEARCH(MOD,UPATH,PEFF,TEFF,AWVL,BWVL,           
!!     $                           CWVL,DWVL,TRANS)                         
               END IF                                                     
                                                                          
            END IF                                                        
 
            TRANS=EXP(-1.66*ABSCOEFF) !! KM  Modif

!            write(*,*) IUP, IDWN, PMID (IDWN+1), ABSCOEFF
!          write (*,*) IDWN, IUP, TRANS, UPATH 
!          write (*,*) PMID(IDWN),PMID(IUP)

C          Include cloud transmittance for cloudy skies                   
            IF (CLOUD.EQ.1) THEN                                          
               TRANS=TRANS*MTRCL(IDWN,IUP)                                
            END IF                                                        

!!            TRANS=EXP(-UPATH * 1.66*ABSLW1) !! KM  Modif
                                                                          
            FUP(IUP)=FUP(IUP)-TRANS*DPF(IDWN+1)                           
                                                                          
         END DO                                                           
      END DO                                                              
                                                                          
                                                                          
C--------------------------------------                                   
C Calculations for the downward fluxes                                    
C--------------------------------------                                   
                                                                          
      MOD=2                                                               
                                                                          
C Calculate arrays PF and DPF used in the calculation of the dn. fluxes   
C PF:  Integrated Planck function at mid-points                           
C DPF: Differences of the PFs                                             

      DO ILAY=1,NLEV                                                      
         PF(ILAY)=PLANCK(0,TMID(ILAY))                                    
      END DO   
                                                                          
      DPF(NLEV)=PLANCK(0,TMID(NLEV))                                      
                                                                          
      DO ILAY=1,NLEV-1                                                    
         DPF(ILAY)=PLANCK(0,TMID(ILAY))-PLANCK(0,TMID(ILAY+1))            
      END DO                                                              
                                                                          
C Initialise the irradiances                                              
                                                                          
      DO ILAY=1,NLEV                                                      
         FDWN(ILAY-1)=PF(ILAY)                                            
      END DO 

                                                                          
C For each layer, calculate the transmittance and update the irradiances  
                                                                          
      DO IUP=1,NLEV                                                       
                                                                          
         DO IDWN=0,IUP-1                                                  
                                                                          
C          Calculate absorber amount in the path between IDWN and IUP     
            UPATH=0.0                                                     
            DO ILAY=IDWN+1,IUP                                            
               UPATH=UPATH+ULAY(1,ILAY)                                   
            END DO                                                        
            UH2O(IDWN,IUP)=UPATH                                          
                                                                          
C KM Modif: Integrating Absorption along UPATH. Diffusive factor 1.66 applied
            IF(NLWMODEL.EQ.1) THEN   

            ABSCOEFF=0.0
               DO ILAY=IDWN+1,IUP                                         
                     IF (OPACIR_POWERLAW.eq.0) THEN !MTR Modif to avoid exponent
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY) 
                     ELSE IF (OPACIR_POWERLAW.eq.1) THEN
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)
     &                 *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES)) 
                     ELSE IF (OPACIR_POWERLAW.eq.2) THEN
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)
     &                 *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
     &                 *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES)) 
                     ELSE IF (OPACIR_POWERLAW.eq.3) THEN
                       ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY)
     &                 *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
     &                 *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
     &                 *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES))
                     ELSE 
                        ABSCOEFF=ABSCOEFF+CABSLW1*ULAY(1,ILAY) 
     &                *MAX(1e-6,(PMID(ILAY)/OPACIR_REFPRES)**OPACIR_POWERLAW) 

                     ENDIF
!KM Modif for Saturated limit (Hourdin): use total UPATH outside integral 
               END DO                                                     
!KM Modif for Saturated limit (Hourdin): outside loop since no integral done
!            ABSCOEFF=CABSLW1*SQRT(UPATH*(PMID(IUP)+PMID(IDWN+1))/2e5)     
!!C           IF(IDWN.EQ.0) ABSCOEFF=CABSLW1*SQRT(UPATH*(PMID(IUP)+1e5)/2e5) 

            ENDIF

            IF(NLWMODEL.EQ.2) THEN

               WRITE(2,*) 'NEED IMPLEMENTATION OF PLANCK MEAN OPACITIES'
               CALL ABORT

               ABSCOEFF=0.0
               DO ILAY=IDWN+1,IUP                                         
                  ABSCOEFF=ABSCOEFF+ CABSLW1*ULAY(1,ILAY)   
               END DO                                                     

            ENDIF

C          Calculate the effective pressure of the path                   
C            PEFF=0.0                                                      
C            DO ILAY=IDWN+1,IUP                                            
C               PEFF=PEFF+PMID(ILAY)*ULAY(1,ILAY)                          
C            END DO                                                        
C            PEFF=PEFF/UPATH                                               
C            PH2O(IDWN,IUP)=PEFF                                           
                                                                          
C          Calculate the effective temperature of the path                
C            TEFF=0.0                                                      
C            DO ILAY=IDWN+1,IUP                                            
C               TEFF=TEFF+TMID(ILAY)*ULAY(1,ILAY)                          
C            END DO                                                        
C            TEFF=TEFF/UPATH                                               
C            TH2O(IDWN,IUP)=TEFF                                           
                                                                          
            IF (CKD.NE.1) THEN                                            
                                                                          
!!               CALL LINSEARCH(MOD,UPATH,PEFF,TEFF,AWVL,BWVL,              
!!     $                        CWVL,DWVL,TRANS)                            
                                                                          
            ELSE                                                          
                                                                          
C--- Include continuum if required                                        
                                                                          
C          Calculate effective water vapour amount                        
C               WVEFF=0.0                                                  
C               DO ILAY=IDWN+1,IUP                                         
C                  WVEFF=WVEFF+(H2O(ILAY)**2)*PMID(ILAY)*DP(ILAY)*         
C     $                  EXP(1800.0*((1.0/TMID(ILAY))-CWVC))               
C               END DO                                                     
C               WVEFF=(1.6178235D-06)*WVEFF                                
                                                                          
               IF (WVEFF.NE.0.0) THEN                                     
                  IF ((LOG10(WVEFF).LT.(-4.5)).AND.                       
     $               (UPATH.LT.1.0D-02)) THEN                             
!!                     CALL LINSEARCH(MOD,UPATH,PEFF,TEFF,AWVL,BWVL,        
!!     $                              CWVL,DWVL,TRANS)                      
                  ELSE                                                    
!!                     CALL CONSEARCH(MOD,UPATH,PEFF,TEFF,WVEFF,AWV,BWV,    
!!     $                              CWV,DWV,TRANS)                        
                  END IF                                                  
               ELSE                                                       
!!                  CALL LINSEARCH(MOD,UPATH,PEFF,TEFF,AWVL,BWVL,           
!!     $                           CWVL,DWVL,TRANS)                         
               END IF                                                     
                                                                          
            END IF                                                        
                                                                          
            TRANS=EXP(-1.66*ABSCOEFF) !! KM  Modif

C          Include cloud transmittance for cloudy skies                   
            IF (CLOUD.EQ.1) THEN                                          
               TRANS=TRANS*TRCL(IDWN,IUP)                                 
            END IF                                                        

!!            TRANS=EXP(-UPATH * 1.66*ABSLW1) !!KM Modif
                                  
!!            if (IDWN.eq.0.and.IUP.eq.NLEV) 
!!          write (*,*) IDWN, IUP, TRANS, UPATH 
                                        
            FDWN(IDWN)=FDWN(IDWN)-TRANS*DPF(IUP)                          
                                                                          
         END DO                                                           
      END DO                                                              
                                                                          
      END                                                                 
                                                                          
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                         SUBROUTINE LINSEARCH                        *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      SUBROUTINE LINSEARCH(MOD,UPATH,PEFF,TEFF,AWVL,BWVL,                 
     $                     CWVL,DWVL,TRANS)                               
                                                                          
C Subroutine LINSEARCH calculates the water vapour line transmittane of   
C a homogeneous  atmospheric path, using the data of pre-computed         
C tables stored in 'wvlin.f'                                              
C                                                                         
C INPUT:  Index for transmittance type (MOD)                              
C         Absorber amount (UPATH)                                         
C         Pressure (PEFF)                                                 
C         Pre-computed tables (AWVL,BWVL,CWVL,DWVL)                       
C         Temperature (TEFF)                                              
C                                                                         
C OUTPUT: Transmittance (TRANS)                                           
                                                                          
      IMPLICIT NONE                                                       
                                                                          
C-----------------                                                        
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER MOD                                                         
      REAL UPATH,PEFF,TEFF                                                
      REAL AWVL(2,34,28),BWVL(2,34,28)                                    
      REAL CWVL(2,34,28),DWVL(2,34,28)                                    
                                                                          
C-----------------                                                        
C Output Variable                                                         
C-----------------                                                        
                                                                          
      REAL TRANS                                                          
                                                                          
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER IABS,IPRES,NABS,NPRES                                       
      PARAMETER(NABS=34,NPRES=28)                                         
      REAL USTR(NABS),PSTR(NPRES)  ! Stored values of absorber            
                                     ! amount and pressure                
                                                                          
      REAL TRU,DTU,TRP,DTP         ! Variables used for interpolation     
                                                                          
                                                                          
C Stored values of absorber amount (kg m-2)                               
                                                                          
      DATA(USTR(IABS),IABS=1,NABS)/  1.0D-09, 1.0D-08, 5.0D-08,           
     $                               1.0D-07, 5.0D-07, 1.0D-06,           
     $                               5.0D-06, 1.0D-05, 5.0D-05,           
     $                               1.0D-04, 5.0D-04, 1.0D-03,           
     $                               5.0D-03, 1.0D-02, 2.0D-02,           
     $                               4.0D-02, 6.0D-02, 8.0D-02,           
     $                               1.0D-01, 2.0D-01, 4.0D-01,           
     $                               6.0D-01, 8.0D-01, 1.0, 2.0,          
     $                               4.0, 6.0, 8.0, 10.0, 20.0,           
     $                               40.0, 60.0, 80.0, 100.0/             
                                                                          
                                                                          
C Stored values of pressure (Pa)                                          
                                                                          
      DATA(PSTR(IPRES),IPRES=1,NPRES)/ 1.0, 1.0D+01, 2.0D+01, 4.0D+01,    
     $                                 6.0D+01, 8.0D+01, 1.0D+02,         
     $                                 2.0D+02, 4.0D+02, 6.0D+02,         
     $                                 8.0D+02, 1.0D+03, 2.0D+03,         
     $                                 4.0D+03, 6.0D+03, 8.0D+03,         
     $                                 1.0D+04, 2.0D+04, 3.0D+04,         
     $                                 4.0D+04, 5.0D+04, 6.0D+04,         
     $                                 7.0D+04, 8.0D+04, 9.0D+04,         
     $                                 1.0D+05, 1.1D+05, 1.2D+05/         
                                                                          
                                                                          
C Find the table index for the absorber amount                            
                                                                          
      IABS=INT(NABS/2.)                                                   
                                                                          
      IF (UPATH.NE.USTR(IABS)) THEN                                       
                                                                          
         IF (UPATH.LT.USTR(IABS)) THEN                                    
                                                                          
            IABS=INT(IABS/2.)                                             
                                                                          
            IF (UPATH.LT.USTR(IABS)) THEN                                 
               IABS=1                                                     
               IF (UPATH.GT.USTR(1)) THEN                                 
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (UPATH.NE.USTR(IABS)) THEN                              
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         ELSE                                                             
                                                                          
            IABS=IABS+INT(NABS/4.)                                        
                                                                          
            IF (UPATH.GT.USTR(IABS)) THEN                                 
               IF (UPATH.GE.USTR(NABS)) THEN                              
                  IABS=NABS                                               
               ELSE                                                       
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (UPATH.NE.USTR(IABS)) THEN                              
                  IABS=INT(NABS/2.)                                       
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
C Find the table index for the pressure                                   
                                                                          
      IPRES=INT(NPRES/2.)                                                 
                                                                          
      IF (PEFF.NE.PSTR(IPRES)) THEN                                       
                                                                          
         IF (PEFF.LT.PSTR(IPRES)) THEN                                    
                                                                          
            IPRES=INT(IPRES/2.)                                           
                                                                          
            IF (PEFF.LE.PSTR(IPRES)) THEN                                 
               IPRES=2                                                    
               IF (PEFF.GT.PSTR(1)) THEN                                  
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (PEFF.NE.PSTR(IPRES)) THEN                              
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         ELSE                                                             
                                                                          
            IPRES=IPRES+INT(NPRES/4.)                                     
                                                                          
            IF (PEFF.GT.PSTR(IPRES)) THEN                                 
               IF (PEFF.GE.PSTR(NPRES)) THEN                              
                  IPRES=NPRES                                             
               ELSE                                                       
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (PEFF.NE.PSTR(IPRES)) THEN                              
                  IPRES=INT(NPRES/2.)                                     
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
C Transmittance calculation                                               
                                                                          
      TEFF=TEFF-250.0                                                     
                                                                          
      TRANS=DWVL(MOD,IABS,IPRES)+TEFF*(CWVL(MOD,IABS,IPRES)+              
     $      TEFF*(BWVL(MOD,IABS,IPRES)+                                   
     $      TEFF*AWVL(MOD,IABS,IPRES)))                                   
                                                                          
C-------------------------------------------                              
C Simple interpolation/extrapolation scheme                               
C-------------------------------------------                              
                                                                          
      IF (IABS.GT.1) THEN  ! Avoid extrapolation for small abs. amounts   
                                                                          
         IF (UPATH.NE.USTR(IABS)) THEN                                    
                                                                          
            TRU=DWVL(MOD,IABS-1,IPRES)+TEFF*(CWVL(MOD,IABS-1,IPRES)+      
     $          TEFF*(BWVL(MOD,IABS-1,IPRES)+                             
     $          TEFF*AWVL(MOD,IABS-1,IPRES)))                             
            DTU=(TRANS-TRU)*                                              
     $          (UPATH-USTR(IABS))/(USTR(IABS)-USTR(IABS-1))              
            TRU=TRANS+DTU       ! Correction for the abs. amount          
            IF (TRU.LT.0.0) THEN                                          
               TRU=0.0          ! Very high abs. amounts                  
            END IF                                                        
                                                                          
         ELSE                   ! Avoid interpolation                     
            TRU=TRANS                                                     
         END IF                                                           
                                                                          
         IF (PEFF.NE.PSTR(IPRES)) THEN                                    
                                                                          
            TRP=DWVL(MOD,IABS,IPRES-1)+TEFF*(CWVL(MOD,IABS,IPRES-1)+      
     $          TEFF*(BWVL(MOD,IABS,IPRES-1)+                             
     $          TEFF*AWVL(MOD,IABS,IPRES-1)))                             
            DTP=(TRANS-TRP)*                                              
     $          (PEFF-PSTR(IPRES))/(PSTR(IPRES)-PSTR(IPRES-1))            
            TRP=TRU+DTP         ! Correction for the pressure             
            IF (TRP.LT.0.0) THEN                                          
               TRP=0.0          ! Very high pressures                     
            END IF                                                        
            TRANS=TRP                                                     
                                                                          
         ELSE                   ! Avoid interpolation                     
            TRANS=TRU                                                     
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
                                                                          
      END                                                                 
                                                                          
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                         SUBROUTINE CONSEARCH                        *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      SUBROUTINE CONSEARCH(MOD,UPATH,PEFF,TEFF,WVEFF,AWV,BWV,             
     $                     CWV,DWV,TRANS)                                 
                                                                          
C Subroutine CONSEARCH calculates the water vapour transmittane of a      
C homogeneous atmospheric path for both line+continuum absorption over    
C the whole thermal IR spectrum (0-3000 cm-1). The subroutine uses the    
C data of pre-computed tables stored in 'wvcon.f'                         
C                                                                         
C INPUT:  Index for transmittance type (MOD)                              
C         Absorber amount (UPATH)                                         
C         Pressure (PEFF)                                                 
C         Temperature (TEFF)                                              
C         Effective water vapour abs. amount (WVEFF)                      
C         Pre-computed tables (AWV,BWV,CWV,DWV)                           
C                                                                         
C OUTPUT: Transmittance (TRANS)                                           
                                                                          
      IMPLICIT NONE                                                       
                                                                          
C-----------------                                                        
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER MOD                                                         
      REAL UPATH,PEFF,TEFF,WVEFF                                          
      REAL AWV(2,41,29,22),BWV(2,41,29,22)                                
      REAL CWV(2,41,29,22),DWV(2,41,29,22)                                
                                                                          
C-----------------                                                        
C Output Variable                                                         
C-----------------                                                        
                                                                          
      REAL TRANS                                                          
                                                                          
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER IABS,IPRES,IWV,NABS,NPRES,NWV                               
      PARAMETER(NABS=29,NPRES=22,NWV=41)                                  
      REAL USTR(NABS)          ! Stored values of abs. amount,            
      REAL PSTR(NPRES)         ! pressure and                             
      REAL WVSTR(NWV)          ! logarithm of eff. H2O amount             
                                                                          
      REAL TR(0:2),TRU,TRP     ! Variables used                           
      REAL DTU,DTP,DTU1,DTU2   ! for interpolation                        
      REAL SLP, INTCPT                                                    
                                                                          
C Stored values of absorber amount (kg m-2)                               
                                                                          
      DATA(USTR(IABS),IABS=1,NABS)/                                       
     $                               1.0D-06, 5.0D-06, 1.0D-05,           
     $                               5.0D-05, 1.0D-04, 5.0D-04,           
     $                               1.0D-03, 5.0D-03, 1.0D-02,           
     $                               2.0D-02, 4.0D-02, 6.0D-02,           
     $                               8.0D-02, 1.0D-01, 2.0D-01,           
     $                               4.0D-01, 6.0D-01, 8.0D-01,           
     $                               1.0, 2.0, 4.0, 6.0, 8.0,             
     $                               10.0, 20.0, 40.0, 60.0,              
     $                               80.0, 100.0/                         
                                                                          
                                                                          
C Stored values of pressure (Pa)                                          
                                                                          
      DATA(PSTR(IPRES),IPRES=1,NPRES)/ 1.0D+02, 2.0D+02, 4.0D+02,         
     $                                 6.0D+02, 8.0D+02, 1.0D+03,         
     $                                 2.0D+03, 4.0D+03, 6.0D+03,         
     $                                 8.0D+03, 1.0D+04, 2.0D+04,         
     $                                 3.0D+04, 4.0D+04, 5.0D+04,         
     $                                 6.0D+04, 7.0D+04, 8.0D+04,         
     $                                 9.0D+04, 1.0D+05, 1.1D+05,         
     $                                 1.2D+05/                           
                                                                          
C Stored values of the logarithm of the H2O effective amount              
                                                                          
      DATA(WVSTR(IWV),IWV=1,NWV)/ -4.6, -4.4, -4.2, -4.0, -3.8, -3.6,     
     $                            -3.4, -3.2, -3.0, -2.9, -2.8, -2.7,     
     $                            -2.6, -2.5, -2.4, -2.3, -2.2, -2.1,     
     $                            -2.0, -1.9, -1.8, -1.7, -1.6, -1.5,     
     $                            -1.4, -1.3, -1.2, -1.1, -1.0, -0.9,     
     $                            -0.8, -0.7, -0.6, -0.5, -0.4, -0.3,     
     $                            -0.2, -0.1,  0.0,  0.1,  0.2/           
                                                                          
      WVEFF=LOG10(WVEFF)                                                  
                                                                          
C Find the table index for the mixing ratio                               
                                                                          
      IWV=INT(NWV/2.)                                                     
                                                                          
      IF (WVEFF.NE.WVSTR(IWV)) THEN                                       
                                                                          
         IF (WVEFF.LT.WVSTR(IWV)) THEN                                    
                                                                          
            IWV=INT(IWV/2.)                                               
                                                                          
            IF (WVEFF.LT.WVSTR(IWV)) THEN                                 
               IWV=1                                                      
               DO WHILE(WVEFF.GT.WVSTR(IWV))                              
                  IWV=IWV+1                                               
               END DO                                                     
            ELSE                                                          
               IF (WVEFF.NE.WVSTR(IWV)) THEN                              
                  DO WHILE(WVEFF.GT.WVSTR(IWV))                           
                     IWV=IWV+1                                            
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         ELSE                                                             
                                                                          
            IWV=IWV+INT(NWV/4.)                                           
                                                                          
            IF (WVEFF.GT.WVSTR(IWV)) THEN                                 
               IF (WVEFF.GE.WVSTR(NWV)) THEN                              
                  IWV=NWV                                                 
               ELSE                                                       
                  DO WHILE(WVEFF.GT.WVSTR(IWV))                           
                     IWV=IWV+1                                            
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (WVEFF.NE.WVSTR(IWV)) THEN                              
                  IWV=INT(NWV/2.)                                         
                  DO WHILE(WVEFF.GT.WVSTR(IWV))                           
                     IWV=IWV+1                                            
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
C Find the table index for the absorber amount                            
                                                                          
      IABS=INT(NABS/2.)                                                   
                                                                          
      IF (UPATH.NE.USTR(IABS)) THEN                                       
                                                                          
         IF (UPATH.LT.USTR(IABS)) THEN                                    
                                                                          
            IABS=INT(IABS/2.)                                             
                                                                          
            IF (UPATH.LT.USTR(IABS)) THEN                                 
               IABS=1                                                     
               IF (UPATH.GT.USTR(1)) THEN                                 
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (UPATH.NE.USTR(IABS)) THEN                              
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         ELSE                                                             
                                                                          
            IABS=IABS+INT(NABS/4.)                                        
                                                                          
            IF (UPATH.GT.USTR(IABS)) THEN                                 
               IF (UPATH.GE.USTR(NABS)) THEN                              
                  IABS=NABS                                               
               ELSE                                                       
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (UPATH.NE.USTR(IABS)) THEN                              
                  IABS=INT(NABS/2.)                                       
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
C Find the table index for the pressure                                   
                                                                          
      IF(PEFF.LT.PSTR(1))THEN                                             
        IPRES=1                                                           
        WRITE(2,*)'WARNING: Top level pressure < 1mbar. Out               
     :             of range of pre-computed table.'                       
      ELSE                                                                
                                                                          
        IPRES=INT(NPRES/2.)                                               
                                                                          
        IF (PEFF.NE.PSTR(IPRES)) THEN                                     
                                                                          
           IF (PEFF.LT.PSTR(IPRES)) THEN                                  
                                                                          
              IPRES=INT(IPRES/2.)                                         
                                                                          
              IF (PEFF.LE.PSTR(IPRES)) THEN                               
                 IPRES=2                                                  
                 IF (PEFF.GT.PSTR(1)) THEN                                
                    DO WHILE(PEFF.GT.PSTR(IPRES))                         
                       IPRES=IPRES+1                                      
                    END DO                                                
                 END IF                                                   
              ELSE                                                        
                 IF (PEFF.NE.PSTR(IPRES)) THEN                            
                    DO WHILE(PEFF.GT.PSTR(IPRES))                         
                       IPRES=IPRES+1                                      
                    END DO                                                
                 END IF                                                   
              END IF                                                      
                                                                          
           ELSE                                                           
                                                                          
              IPRES=IPRES+INT(NPRES/4.)                                   
                                                                          
              IF (PEFF.GT.PSTR(IPRES)) THEN                               
                 IF (PEFF.GE.PSTR(NPRES)) THEN                            
                    IPRES=NPRES                                           
                 ELSE                                                     
                    DO WHILE(PEFF.GT.PSTR(IPRES))                         
                       IPRES=IPRES+1                                      
                    END DO                                                
                 END IF                                                   
              ELSE                                                        
                 IF (PEFF.NE.PSTR(IPRES)) THEN                            
                    IPRES=INT(NPRES/2.)                                   
                    DO WHILE(PEFF.GT.PSTR(IPRES))                         
                       IPRES=IPRES+1                                      
                    END DO                                                
                 END IF                                                   
              END IF                                                      
                                                                          
           END IF                                                         
                                                                          
        END IF                                                            
      ENDIF                                                               
C-------------------------------------------                              
C Simple interpolation/extrapolation scheme                               
C-------------------------------------------                              
                                                                          
C PIERS' fix                                                              
      IF (IWV.LE.1) IWV=2                                                 
C--- MOD=1                                                                
                                                                          
C An interpolated transmittance (with respect to pressure and absorber    
C amount) is calculated on two successive levels of effective H2O amount  
C and the final transmittance results from the interpolation between      
C the 2 values                                                            
                                                                          
      TEFF=TEFF-250.0                                                     
                                                                          
      IF (MOD.EQ.1) THEN                                                  
                                                                          
         TR(0)=DWV(MOD,IWV,IABS,IPRES)+                                   
     $         TEFF*(CWV(MOD,IWV,IABS,IPRES)+                             
     $         TEFF*(BWV(MOD,IWV,IABS,IPRES)+                             
     $         TEFF*AWV(MOD,IWV,IABS,IPRES)))                             
                                                                          
         TRANS=TR(0)                                                      
                                                                          
         IF (IABS.GT.1) THEN    ! Avoid extrap. for small abs. amounts    
                                                                          
            IF (UPATH.NE.USTR(IABS)) THEN                                 
                                                                          
               TRU=DWV(MOD,IWV,IABS-1,IPRES)+                             
     $             TEFF*(CWV(MOD,IWV,IABS-1,IPRES)+                       
     $             TEFF*(BWV(MOD,IWV,IABS-1,IPRES)+                       
     $             TEFF*AWV(MOD,IWV,IABS-1,IPRES)))                       
               DTU=(TRANS-TRU)*                                           
     $             (UPATH-USTR(IABS))/(USTR(IABS)-USTR(IABS-1))           
               TRU=TR(0)+DTU    ! Correction for the abs. amount          
               IF (TRU.LT.0.0) THEN                                       
                  TRU=0.0       ! Very high abs. amounts                  
               END IF                                                     
                                                                          
            ELSE                ! Avoid interpolation                     
               TRU=TR(0)                                                  
            END IF                                                        
                                                                          
            IF (PEFF.NE.PSTR(IPRES)) THEN                                 
                                                                          
               TRP=DWV(MOD,IWV,IABS,IPRES-1)+                             
     $             TEFF*(CWV(MOD,IWV,IABS,IPRES-1)+                       
     $             TEFF*(BWV(MOD,IWV,IABS,IPRES-1)+                       
     $             TEFF*AWV(MOD,IWV,IABS,IPRES-1)))                       
               DTP=(TRANS-TRP)*                                           
     $              (PEFF-PSTR(IPRES))/(PSTR(IPRES)-PSTR(IPRES-1))        
               TRP=TRU+DTP      ! Correction for the pressure             
               IF (TRP.LT.0.0) THEN                                       
                  TRP=0.0       ! Very high pressures                     
               END IF                                                     
               TR(1)=TRP                                                  
                                                                          
            ELSE                ! Avoid interpolation                     
               TR(1)=TRU                                                  
            END IF                                                        
                                                                          
            TRANS=TR(1)                                                   
                                                                          
         END IF                                                           
                                                                          
C Interpolation between effective absorber amounts                        
                                                                          
         IF (WVEFF.NE.WVSTR(IWV)) THEN                                    
                                                                          
            TR(0)=DWV(MOD,IWV-1,IABS,IPRES)+                              
     $            TEFF*(CWV(MOD,IWV-1,IABS,IPRES)+                        
     $            TEFF*(BWV(MOD,IWV-1,IABS,IPRES)+                        
     $            TEFF*AWV(MOD,IWV-1,IABS,IPRES)))                        
                                                                          
            IF (IABS.GT.1) THEN ! Avoid extrap. for small abs. amounts    
                                                                          
               IF (UPATH.NE.USTR(IABS)) THEN                              
                                                                          
                  TRU=DWV(MOD,IWV-1,IABS-1,IPRES)+                        
     $                TEFF*(CWV(MOD,IWV-1,IABS-1,IPRES)+                  
     $                TEFF*(BWV(MOD,IWV-1,IABS-1,IPRES)+                  
     $                TEFF*AWV(MOD,IWV-1,IABS-1,IPRES)))                  
                  DTU=(TRANS-TRU)*                                        
     $                (UPATH-USTR(IABS))/(USTR(IABS)-USTR(IABS-1))        
                  TRU=TR(0)+DTU   ! Correction for the abs. amount        
                  IF (TRU.LT.0.0) THEN                                    
                     TRU=0.0      ! Very high abs. amounts                
                  END IF                                                  
                                                                          
               ELSE               ! Avoid interpolation                   
                  TRU=TR(0)                                               
               END IF                                                     
                                                                          
               IF (PEFF.NE.PSTR(IPRES)) THEN                              
                                                                          
                  TRP=DWV(MOD,IWV-1,IABS,IPRES-1)+                        
     $                TEFF*(CWV(MOD,IWV-1,IABS,IPRES-1)+                  
     $                TEFF*(BWV(MOD,IWV-1,IABS,IPRES-1)+                  
     $                TEFF*AWV(MOD,IWV-1,IABS,IPRES-1)))                  
                  DTP=(TRANS-TRP)*                                        
     $                (PEFF-PSTR(IPRES))/(PSTR(IPRES)-PSTR(IPRES-1))      
                  TRP=TRU+DTP     ! Correction for the pressure           
                  IF (TRP.LT.0.0) THEN                                    
                     TRP=0.0      ! Very high pressures                   
                  END IF                                                  
                  TR(2)=TRP                                               
                                                                          
               ELSE               ! Avoid interpolation                   
                  TR(2)=TRU                                               
               END IF                                                     
                                                                          
            END IF                                                        
                                                                          
            SLP=(TR(1)-TR(2))/(WVSTR(IWV)-WVSTR(IWV-1))                   
            INTCPT=WVSTR(IWV)*TR(2)-WVSTR(IWV-1)*TR(1)                    
            INTCPT=INTCPT/(WVSTR(IWV)-WVSTR(IWV-1))                       
            TRANS=SLP*WVEFF+INTCPT                                        
                                                                          
         END IF                                                           
                                                                          
C--- MOD=2                                                                
                                                                          
C The transmittances between 2 successive levels of the effective H2O     
C are projected on a value on the actual level and interpolation is       
C performed on the actual level. Interpolation between pressure levels    
C is not considered                                                       
                                                                          
C Transmittance calculation                                               
                                                                          
      ELSE                                                                
                                                                          
         TR(1)=DWV(MOD,IWV,IABS,IPRES)+                                   
     $         TEFF*(CWV(MOD,IWV,IABS,IPRES)+                             
     $         TEFF*(BWV(MOD,IWV,IABS,IPRES)+                             
     $         TEFF*AWV(MOD,IWV,IABS,IPRES)))                             
                                                                          
         TR(2)=DWV(MOD,IWV-1,IABS,IPRES)+                                 
     $         TEFF*(CWV(MOD,IWV-1,IABS,IPRES)+                           
     $         TEFF*(BWV(MOD,IWV-1,IABS,IPRES)+                           
     $         TEFF*AWV(MOD,IWV-1,IABS,IPRES)))                           
                                                                          
C Interpolation between effective absorber amounts                        
                                                                          
         IF (WVEFF.NE.WVSTR(IWV)) THEN                                    
            TR(0)=TR(2)+(TR(1)-TR(2))*(WVEFF-WVSTR(IWV-1))/               
     $            (WVSTR(IWV)-WVSTR(IWV-1))                               
         ELSE                                                             
            TR(0)=TR(1)                                                   
         END IF                                                           
                                                                          
C Interpolation between absorber amounts                                  
                                                                          
         IF (IABS.GT.1) THEN    ! Avoid extrap. for small abs. amounts    
                                                                          
            IF (UPATH.NE.USTR(IABS)) THEN                                 
                                                                          
               TRU=DWV(MOD,IWV,IABS-1,IPRES)+                             
     $             TEFF*(CWV(MOD,IWV,IABS-1,IPRES)+                       
     $             TEFF*(BWV(MOD,IWV,IABS-1,IPRES)+                       
     $             TEFF*AWV(MOD,IWV,IABS-1,IPRES)))                       
               DTU1=(TR(1)-TRU)*                                          
     $              (UPATH-USTR(IABS))/(USTR(IABS)-USTR(IABS-1))          
                                                                          
               TRU=DWV(MOD,IWV-1,IABS-1,IPRES)+                           
     $             TEFF*(CWV(MOD,IWV-1,IABS-1,IPRES)+                     
     $             TEFF*(BWV(MOD,IWV-1,IABS-1,IPRES)+                     
     $             TEFF*AWV(MOD,IWV-1,IABS-1,IPRES)))                     
               DTU2=(TR(2)-TRU)*                                          
     $              (UPATH-USTR(IABS))/(USTR(IABS)-USTR(IABS-1))          
                                                                          
               IF (WVEFF.NE.WVSTR(IWV)) THEN                              
                  DTU=DTU2+(DTU1-DTU2)*(WVEFF-WVSTR(IWV-1))/              
     $               (WVSTR(IWV)-WVSTR(IWV-1))                            
               ELSE                                                       
                  DTU=0.0                                                 
               END IF                                                     
                                                                          
               TRU=TR(0)+DTU    ! Correction for the abs. amount          
               IF (TRU.LT.0.0) THEN                                       
                  TRU=0.0       ! Very high abs. amounts                  
               END IF                                                     
                                                                          
            ELSE                ! Avoid interpolation                     
               TRU=TR(0)                                                  
            END IF                                                        
                                                                          
            TRANS =TRU                                                    
                                                                          
         ELSE                                                             
                                                                          
            TRANS=TR(0)                                                   
                                                                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
      IF (TRANS.LT.0.0) THEN                                              
         TRANS=0.0                                                        
      END IF                                                              
                                                                          
      IF (TRANS.GT.1.0) THEN                                              
         TRANS=1.0                                                        
      END IF                                                              
                                                                          
      END                                                                 
                                                                          
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                         SUBROUTINE GASSEARCH                        *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      SUBROUTINE GASSEARCH(IGAS,IBND,UPATH,PEFF,TEFF,AH2O,BH2O,CH2O,      
     $                     DH2O,AGAS,BGAS,CGAS,DGAS,AN03,BN03,CN03,       
     $                     DN03,AC02,BC02,CC02,DC02,AC06,BC06,CC06,       
     $                     DC06,AC07,BC07,CC07,DC07,TRANS)                

C Subroutine GASSEARCH calculates the transmittane of a homogeneous
C atmospheric path, using the data of pre-computed tables stored
C in 'gastab'
C
C INPUT:  Gas index (IGAS)
C         Band index (IBND)
C         Absorber amount (UPATH)
C         Pressure (PEFF)
C         Temperature (TEFF)
C         Pre-computed tables (AH2O,BH2O,CH2O,DH2O,AGAS,BGAS,CGAS,DGAS,
C                              AN03,BN03,CN03,DN03,AC02,BC02,CC02,DC02,
C                              AC06,BC06,CC06,DC06,AC07,BC07,CC07,DC07)
C
C OUTPUT: Transmittance (TRANS)
                                                                          
      IMPLICIT NONE                                                       
C***********************************************************************  
C*                                                                     *  
C*                         P A R R A Y                                 *  
C*                                                                     *  
C***********************************************************************  
C-----------------------------------------------------------------------  
C In this part of the code the user can set the dimensions of arrays      
C used for the calculations, by setting the values in the paramter        
C statements                                                              
C-----------------------------------------------------------------------  
                      
C      include 'params.i'
                                                    
      INTEGER MXGAS,MXLEV,MXBAND,MXCL                                     
                                                                          
      PARAMETER(MXGAS=8)     ! Maximum number of gases                    
                                                                          
      PARAMETER(MXLEV=5)     ! Maximum number of levels in the            
                             ! atmosphere (not including the surface)     
                                                                          
      PARAMETER(MXBAND=9)    ! Maximum number of spectral bands (not      
                             ! including the whole spectrum, 0-3000cm-1)  
                                                                          
      PARAMETER(MXCL=3)      ! Maximum number of cloud types              
                                                                          
C-----------------------------------------------------------------------  
                                                                          
C-----------------                                                        
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER IGAS,IBND                                                   
      REAL UPATH,PEFF,TEFF                                                
                                                                          
      REAL AH2O(MXBAND,48,28),AGAS(MXBAND-1,48,28)                        
      REAL BH2O(MXBAND,48,28),BGAS(MXBAND-1,48,28)                        
      REAL CH2O(MXBAND,48,28),CGAS(MXBAND-1,48,28)                        
      REAL DH2O(MXBAND,48,28),DGAS(MXBAND-1,48,28)                        
      REAL AN03(48,28),BN03(48,28),CN03(48,28),DN03(48,28)                
      REAL AC02(48,28),BC02(48,28),CC02(48,28),DC02(48,28)                
      REAL AC06(48,28),BC06(48,28),CC06(48,28),DC06(48,28)                
      REAL AC07(48,28),BC07(48,28),CC07(48,28),DC07(48,28)                
                                                                          
C-----------------                                                        
C Output Variable                                                         
C-----------------                                                        
                                                                          
      REAL TRANS                                                          
                                                                          
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER IABS,IPRES,NABS,NPRES                                       
      PARAMETER(NABS=48,NPRES=28)                                         
      REAL USTR(NABS),PSTR(NPRES)  ! Stored values of absorber            
                                     ! amount and pressure                
                                                                          
      REAL FU,FFU,FP,FFP          ! Variables                             
      REAL A0,A1,A2,B0,B1,B2      ! used for                              
      REAL C0,C1,C2,D0,D1,D2      ! interpolation                         
                                                                          
C Stored values of absorber amount (kg m-2)                               
                                                                          
      DATA(USTR(IABS),IABS=1,NABS)/  1.0D-08, 1.0D-07, 5.0D-07,           
     $                               1.0D-06, 2.0D-06, 4.0D-06,           
     $                               6.0D-06, 8.0D-06, 1.0D-05,           
     $                               2.0D-05, 4.0D-05, 6.0D-05,           
     $                               8.0D-05, 1.0D-04, 2.0D-04,           
     $                               4.0D-04, 6.0D-04, 8.0D-04,           
     $                               1.0D-03, 2.0D-03, 4.0D-03,           
     $                               6.0D-03, 8.0D-03, 1.0D-02,           
     $                               2.0D-02, 4.0D-02, 6.0D-02,           
     $                               8.0D-02, 1.0D-01, 2.0D-01,           
     $                               4.0D-01, 6.0D-01, 8.0D-01,           
     $                               1.0, 2.0, 4.0, 6.0, 8.0,             
     $                               1.0D+01, 2.0D+01, 4.0D+01,           
     $                               6.0D+01, 8.0D+01, 1.0D+02,           
     $                               2.0D+02, 4.0D+02, 6.0D+02,           
     $                               8.0D+02/                             
                                                                          
C Stored values of pressure (Pa)                                          
                                                                          
      DATA(PSTR(IPRES),IPRES=1,NPRES)/ 1.0, 1.0D+01, 2.0D+01, 4.0D+01,    
     $                                 6.0D+01, 8.0D+01, 1.0D+02,         
     $                                 2.0D+02, 4.0D+02, 6.0D+02,         
     $                                 8.0D+02, 1.0D+03, 2.0D+03,         
     $                                 4.0D+03, 6.0D+03, 8.0D+03,         
     $                                 1.0D+04, 2.0D+04, 3.0D+04,         
     $                                 4.0D+04, 5.0D+04, 6.0D+04,         
     $                                 7.0D+04, 8.0D+04, 9.0D+04,         
     $                                 1.0D+05, 1.1D+05, 1.2D+05/         
                                                                          
                                                                          
      FU=99.                                                              
      FP=99.                                                              
                                                                          
C Find the table index for the absorber amount                            
                                                                          
      IABS=INT(NABS/2.)                                                   
                                                                          
      IF (UPATH.NE.USTR(IABS)) THEN                                       
                                                                          
         IF (UPATH.LT.USTR(IABS)) THEN                                    
                                                                          
            IABS=INT(IABS/2.)                                             
                                                                          
            IF (UPATH.LT.USTR(IABS)) THEN                                 
               IABS=1                                                     
               IF (UPATH.GT.USTR(1)) THEN                                 
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
               IF (UPATH.LE.USTR(1)) THEN                                 
                  IABS=2                                                  
                  FFU=0.                                                  
                  FU=1.                                                   
               END IF                                                     
            ELSE                                                          
               IF (UPATH.NE.USTR(IABS)) THEN                              
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         ELSE                                                             
                                                                          
            IABS=IABS+INT(NABS/4.)                                        
                                                                          
            IF (UPATH.GT.USTR(IABS)) THEN                                 
               IF (UPATH.GE.USTR(NABS)) THEN                              
                  IABS=NABS                                               
                  FFU=1.                                                  
                  FU=0.                                                   
               ELSE                                                       
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (UPATH.NE.USTR(IABS)) THEN                              
                  IABS=INT(NABS/2.)                                       
                  DO WHILE(UPATH.GT.USTR(IABS))                           
                     IABS=IABS+1                                          
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
                                                                          
C Find the table index for the pressure                                   
                                                                          
      IPRES=INT(NPRES/2.)                                                 
                                                                          
      IF (PEFF.NE.PSTR(IPRES)) THEN                                       
                                                                          
         IF (PEFF.LT.PSTR(IPRES)) THEN                                    
                                                                          
            IPRES=INT(IPRES/2.)                                           
                                                                          
            IF (PEFF.LE.PSTR(IPRES)) THEN                                 
               IPRES=2                                                    
               IF (PEFF.GT.PSTR(1)) THEN                                  
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
               IF (PEFF.LE.PSTR(1)) THEN                                  
                  FFP=0.                                                  
                  FP=1.                                                   
               END IF                                                     
            ELSE                                                          
               IF (PEFF.NE.PSTR(IPRES)) THEN                              
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         ELSE                                                             
                                                                          
            IPRES=IPRES+INT(NPRES/4.)                                     
                                                                          
            IF (PEFF.GT.PSTR(IPRES)) THEN                                 
               IF (PEFF.GE.PSTR(NPRES)) THEN                              
                  IPRES=NPRES                                             
                  FFP=1.                                                  
                  FP=0.                                                   
               ELSE                                                       
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
            ELSE                                                          
               IF (PEFF.NE.PSTR(IPRES)) THEN                              
                  IPRES=INT(NPRES/2.)                                     
                  DO WHILE(PEFF.GT.PSTR(IPRES))                           
                     IPRES=IPRES+1                                        
                  END DO                                                  
               END IF                                                     
            END IF                                                        
                                                                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
C Transmittance calculation                                               
                                                                          
      TEFF=TEFF-250.0                                                     
                                                                          
C----------------------                                                   
C Interpolation scheme                                                    
C----------------------                                                   
                                                                          
      IF ((FU.GT.98.).AND.(FP.GT.98)) THEN                                
         FU=(USTR(IABS)-UPATH)/(USTR(IABS)-USTR(IABS-1))                  
         FP=(PSTR(IPRES)-PEFF)/(PSTR(IPRES)-PSTR(IPRES-1))                
         FFU=1.-FU                                                        
         FFP=1.-FP                                                        
      END IF                                                              
                                                                          
      IF (IGAS.EQ.1) THEN                                                 
                                                                          
         A1=AH2O(IBND,IABS,IPRES)*FFP+                                    
     $      AH2O(IBND,IABS,IPRES-1)*FP                                    
         A2=AH2O(IBND,IABS-1,IPRES)*FFP+                                  
     $      AH2O(IBND,IABS-1,IPRES-1)*FP                                  
         A0=A1*FFU+A2*FU                                                  
                                                                          
         B1=BH2O(IBND,IABS,IPRES)*FFP+                                    
     $      BH2O(IBND,IABS,IPRES-1)*FP                                    
         B2=BH2O(IBND,IABS-1,IPRES)*FFP+                                  
     $      BH2O(IBND,IABS-1,IPRES-1)*FP                                  
         B0=B1*FFU+B2*FU                                                  
                                                                          
         C1=CH2O(IBND,IABS,IPRES)*FFP+                                    
     $      CH2O(IBND,IABS,IPRES-1)*FP                                    
         C2=CH2O(IBND,IABS-1,IPRES)*FFP+                                  
     $      CH2O(IBND,IABS-1,IPRES-1)*FP                                  
         C0=C1*FFU+C2*FU                                                  
                                                                          
         D1=DH2O(IBND,IABS,IPRES)*FFP+                                    
     $      DH2O(IBND,IABS,IPRES-1)*FP                                    
         D2=DH2O(IBND,IABS-1,IPRES)*FFP+                                  
     $      DH2O(IBND,IABS-1,IPRES-1)*FP                                  
         D0=D1*FFU+D2*FU                                                  
                                                                          
         TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                             
         GOTO 1976                                                        
                                                                          
      END IF                                                              
                                                                          
      IF (IGAS.EQ.2) THEN                                                 
                                                                          
         IF (IBND.EQ.1) THEN                                              
            A1=AGAS(IBND,IABS,IPRES)*FFP+                                 
     $         AGAS(IBND,IABS,IPRES-1)*FP                                 
            A2=AGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         AGAS(IBND,IABS-1,IPRES-1)*FP                               
            A0=A1*FFU+A2*FU                                               
                                                                          
            B1=BGAS(IBND,IABS,IPRES)*FFP+                                 
     $         BGAS(IBND,IABS,IPRES-1)*FP                                 
            B2=BGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         BGAS(IBND,IABS-1,IPRES-1)*FP                               
            B0=B1*FFU+B2*FU                                               
                                                                          
            C1=CGAS(IBND,IABS,IPRES)*FFP+                                 
     $         CGAS(IBND,IABS,IPRES-1)*FP                                 
            C2=CGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         CGAS(IBND,IABS-1,IPRES-1)*FP                               
            C0=C1*FFU+C2*FU                                               
                                                                          
            D1=DGAS(IBND,IABS,IPRES)*FFP+                                 
     $         DGAS(IBND,IABS,IPRES-1)*FP                                 
            D2=DGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         DGAS(IBND,IABS-1,IPRES-1)*FP                               
            D0=D1*FFU+D2*FU                                               
                                                                          
            TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                          
            GOTO 1976                                                     
         END IF                                                           
                                                                          
         IF (IBND.EQ.6) THEN                                              
            A1=AC06(IABS,IPRES)*FFP+AC06(IABS,IPRES-1)*FP                 
            A2=AC06(IABS-1,IPRES)*FFP+AC06(IABS-1,IPRES-1)*FP             
            A0=A1*FFU+A2*FU                                               
                                                                          
            B1=BC06(IABS,IPRES)*FFP+BC06(IABS,IPRES-1)*FP                 
            B2=BC06(IABS-1,IPRES)*FFP+BC06(IABS-1,IPRES-1)*FP             
            B0=B1*FFU+B2*FU                                               
                                                                          
            C1=CC06(IABS,IPRES)*FFP+CC06(IABS,IPRES-1)*FP                 
            C2=CC06(IABS-1,IPRES)*FFP+CC06(IABS-1,IPRES-1)*FP             
            C0=C1*FFU+C2*FU                                               
                                                                          
            D1=DC06(IABS,IPRES)*FFP+DC06(IABS,IPRES-1)*FP                 
            D2=DC06(IABS-1,IPRES)*FFP+DC06(IABS-1,IPRES-1)*FP             
            D0=D1*FFU+D2*FU                                               
                                                                          
            TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                          
            GOTO 1976                                                     
         END IF                                                           
                                                                          
         IF (IBND.EQ.7) THEN                                              
            A1=AC07(IABS,IPRES)*FFP+AC07(IABS,IPRES-1)*FP                 
            A2=AC07(IABS-1,IPRES)*FFP+AC07(IABS-1,IPRES-1)*FP             
            A0=A1*FFU+A2*FU                                               
                                                                          
            B1=BC07(IABS,IPRES)*FFP+BC07(IABS,IPRES-1)*FP                 
            B2=BC07(IABS-1,IPRES)*FFP+BC07(IABS-1,IPRES-1)*FP             
            B0=B1*FFU+B2*FU                                               
                                                                          
            C1=CC07(IABS,IPRES)*FFP+CC07(IABS,IPRES-1)*FP                 
            C2=CC07(IABS-1,IPRES)*FFP+CC07(IABS-1,IPRES-1)*FP             
            C0=C1*FFU+C2*FU                                               
                                                                          
            D1=DC07(IABS,IPRES)*FFP+DC07(IABS,IPRES-1)*FP                 
            D2=DC07(IABS-1,IPRES)*FFP+DC07(IABS-1,IPRES-1)*FP             
            D0=D1*FFU+D2*FU                                               
                                                                          
            TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                          
            GOTO 1976                                                     
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
      IF ((IGAS.EQ.3).OR.(IGAS.EQ.4).OR.                                  
     $    (IGAS.EQ.7).OR.(IBND.EQ.8)) THEN                                
                                                                          
         A1=AGAS(IBND,IABS,IPRES)*FFP+                                    
     $      AGAS(IBND,IABS,IPRES-1)*FP                                    
         A2=AGAS(IBND,IABS-1,IPRES)*FFP+                                  
     $      AGAS(IBND,IABS-1,IPRES-1)*FP                                  
         A0=A1*FFU+A2*FU                                                  
                                                                          
         B1=BGAS(IBND,IABS,IPRES)*FFP+                                    
     $      BGAS(IBND,IABS,IPRES-1)*FP                                    
         B2=BGAS(IBND,IABS-1,IPRES)*FFP+                                  
     $      BGAS(IBND,IABS-1,IPRES-1)*FP                                  
         B0=B1*FFU+B2*FU                                                  
                                                                          
         C1=CGAS(IBND,IABS,IPRES)*FFP+                                    
     $      CGAS(IBND,IABS,IPRES-1)*FP                                    
         C2=CGAS(IBND,IABS-1,IPRES)*FFP+                                  
     $      CGAS(IBND,IABS-1,IPRES-1)*FP                                  
         C0=C1*FFU+C2*FU                                                  
                                                                          
         D1=DGAS(IBND,IABS,IPRES)*FFP+                                    
     $      DGAS(IBND,IABS,IPRES-1)*FP                                    
         D2=DGAS(IBND,IABS-1,IPRES)*FFP+                                  
     $      DGAS(IBND,IABS-1,IPRES-1)*FP                                  
         D0=D1*FFU+D2*FU                                                  
                                                                          
         TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                             
         GOTO 1976                                                        
                                                                          
      END IF                                                              
                                                                          
      IF (IGAS.EQ.5) THEN                                                 
                                                                          
         IF ((IBND.EQ.4).OR.(IBND.EQ.7)) THEN                             
            A1=AGAS(IBND,IABS,IPRES)*FFP+                                 
     $         AGAS(IBND,IABS,IPRES-1)*FP                                 
            A2=AGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         AGAS(IBND,IABS-1,IPRES-1)*FP                               
            A0=A1*FFU+A2*FU                                               
                                                                          
            B1=BGAS(IBND,IABS,IPRES)*FFP+                                 
     $         BGAS(IBND,IABS,IPRES-1)*FP                                 
            B2=BGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         BGAS(IBND,IABS-1,IPRES-1)*FP                               
            B0=B1*FFU+B2*FU                                               
                                                                          
            C1=CGAS(IBND,IABS,IPRES)*FFP+                                 
     $         CGAS(IBND,IABS,IPRES-1)*FP                                 
            C2=CGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         CGAS(IBND,IABS-1,IPRES-1)*FP                               
            C0=C1*FFU+C2*FU                                               
                                                                          
            D1=DGAS(IBND,IABS,IPRES)*FFP+                                 
     $         DGAS(IBND,IABS,IPRES-1)*FP                                 
            D2=DGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         DGAS(IBND,IABS-1,IPRES-1)*FP                               
            D0=D1*FFU+D2*FU                                               
                                                                          
            TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                          
            GOTO 1976                                                     
         END IF                                                           
                                                                          
         IF (IBND.EQ.3) THEN                                              
            A1=AN03(IABS,IPRES)*FFP+AN03(IABS,IPRES-1)*FP                 
            A2=AN03(IABS-1,IPRES)*FFP+AN03(IABS-1,IPRES-1)*FP             
            A0=A1*FFU+A2*FU                                               
                                                                          
            B1=BN03(IABS,IPRES)*FFP+BN03(IABS,IPRES-1)*FP                 
            B2=BN03(IABS-1,IPRES)*FFP+BN03(IABS-1,IPRES-1)*FP             
            B0=B1*FFU+B2*FU                                               
                                                                          
            C1=CN03(IABS,IPRES)*FFP+CN03(IABS,IPRES-1)*FP                 
            C2=CN03(IABS-1,IPRES)*FFP+CN03(IABS-1,IPRES-1)*FP             
            C0=C1*FFU+C2*FU                                               
                                                                          
            D1=DN03(IABS,IPRES)*FFP+DN03(IABS,IPRES-1)*FP                 
            D2=DN03(IABS-1,IPRES)*FFP+DN03(IABS-1,IPRES-1)*FP             
            D0=D1*FFU+D2*FU                                               
                                                                          
            TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                          
            GOTO 1976                                                     
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
      IF (IGAS.EQ.6) THEN                                                 
                                                                          
         IF (IBND.EQ.5) THEN                                              
            A1=AGAS(IBND,IABS,IPRES)*FFP+                                 
     $         AGAS(IBND,IABS,IPRES-1)*FP                                 
            A2=AGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         AGAS(IBND,IABS-1,IPRES-1)*FP                               
            A0=A1*FFU+A2*FU                                               
                                                                          
            B1=BGAS(IBND,IABS,IPRES)*FFP+                                 
     $         BGAS(IBND,IABS,IPRES-1)*FP                                 
            B2=BGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         BGAS(IBND,IABS-1,IPRES-1)*FP                               
            B0=B1*FFU+B2*FU                                               
                                                                          
            C1=CGAS(IBND,IABS,IPRES)*FFP+                                 
     $         CGAS(IBND,IABS,IPRES-1)*FP                                 
            C2=CGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         CGAS(IBND,IABS-1,IPRES-1)*FP                               
            C0=C1*FFU+C2*FU                                               
                                                                          
            D1=DGAS(IBND,IABS,IPRES)*FFP+                                 
     $         DGAS(IBND,IABS,IPRES-1)*FP                                 
            D2=DGAS(IBND,IABS-1,IPRES)*FFP+                               
     $         DGAS(IBND,IABS-1,IPRES-1)*FP                               
            D0=D1*FFU+D2*FU                                               
                                                                          
            TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                          
            GOTO 1976                                                     
         END IF                                                           
                                                                          
         IF (IBND.EQ.2) THEN                                              
            A1=AC02(IABS,IPRES)*FFP+AC02(IABS,IPRES-1)*FP                 
            A2=AC02(IABS-1,IPRES)*FFP+AC02(IABS-1,IPRES-1)*FP             
            A0=A1*FFU+A2*FU                                               
                                                                          
            B1=BC02(IABS,IPRES)*FFP+BC02(IABS,IPRES-1)*FP                 
            B2=BC02(IABS-1,IPRES)*FFP+BC02(IABS-1,IPRES-1)*FP             
            B0=B1*FFU+B2*FU                                               
                                                                          
            C1=CC02(IABS,IPRES)*FFP+CC02(IABS,IPRES-1)*FP                 
            C2=CC02(IABS-1,IPRES)*FFP+CC02(IABS-1,IPRES-1)*FP             
            C0=C1*FFU+C2*FU                                               
                                                                          
            D1=DC02(IABS,IPRES)*FFP+DC02(IABS,IPRES-1)*FP                 
            D2=DC02(IABS-1,IPRES)*FFP+DC02(IABS-1,IPRES-1)*FP             
            D0=D1*FFU+D2*FU                                               
                                                                          
            TRANS=D0+TEFF*(C0+TEFF*(B0+TEFF*A0))                          
         END IF                                                           
                                                                          
      END IF                                                              
                                                                          
 1976 CONTINUE                                                            
                                                                          
      TEFF=TEFF+250.0                                                     
                                                                          
                                                                          
      END                                                                 
                                                                          
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                     FUNCTION HALTRAN                                *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      REAL FUNCTION HALTRAN (IBND,UPATH)                                  
                                                                          
C This function returns the halocarbon transmittance, calculated by the   
C simple expression Tr=1-a*u, where u is the column amount per molecular  
C weight The values of 'a' were chosed, so that the model can reproduce   
C the radiative forcing of a typical halocarbon (average of CFC-11,       
C CFC-12, CFC-113, HCFC-141b, HFC-134a, HCFC-22) computed by a            
C narrow-band model with a mid-latitude summer profile, for different     
C overlap cases                                                           
C                                                                         
C INPUT:  Band index (IBND)                                               
C         Absorber amount (UPATH)                                         
C OUTPUT: Halocarbon Transmittance (HALTRAN)                              
                                                                          
      IMPLICIT NONE                                                       
                                                                          
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER IBND                                                        
      REAL UPATH                                                          
                                                                          
                                                                          
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      REAL ALPHA(4)                                                       
                                                                          
      DATA ALPHA(1), ALPHA(2), ALPHA(3), ALPHA(4) /                       
     $     1.300D+05, 1.200D+05, 1.200D+05, 1.133D+05/                    
                                                                          
      IF (IBND.EQ.9) THEN                                                 
         HALTRAN=1.0-ALPHA(1)*UPATH                                       
      END IF                                                              
                                                                          
      IF (IBND.EQ.10) THEN                                                
         HALTRAN=1.0-ALPHA(2)*UPATH                                       
      END IF                                                              
                                                                          
      IF (IBND.EQ.11) THEN                                                
         HALTRAN=1.0-ALPHA(3)*UPATH                                       
      END IF                                                              
                                                                          
      IF (IBND.EQ.12) THEN                                                
         HALTRAN=1.0-ALPHA(4)*UPATH                                       
      END IF                                                              
                                                                          
      END                                                                 
                                                                          
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                         SUBROUTINE GASCNT                           *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      SUBROUTINE GASCNT(IBND,WVEFF,TEFF,TRCNT)                            
                                                                          
C Subroutine GASCNT calculates the self-broadened continuum               
C transmittance in the bands of gases                                     
C                                                                         
C INPUT:  Band index (IBND)                                               
C         Effective H2O amount for continuum absorption (WVEFF)           
C         Effective temperature (TEFF)                                    
C                                                                         
C OUTPUT: Self-broadened continuum transmittance (TRCNT)                  
                                                                          
      IMPLICIT NONE                                                       
                                                                          
C-----------------                                                        
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER IBND                                                        
      REAL WVEFF,TEFF                                                     
                                                                          
C-----------------                                                        
C Output Variable                                                         
C-----------------                                                        
                                                                          
      REAL TRCNT                                                          
                                                                          
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      REAL BND1    ! Abs. Xsection in band 1                              
      REAL BND2    ! Abs. Xection  in band 2                              
      REAL BND3    ! Abs. Xsection in band 3                              
      REAL BND4    ! Abs. Xsection in band 4                              
      REAL BND5    ! Abs. Xsection in band 5                              
      REAL BND6    ! Abs. Xsection in band 6                              
      REAL BND7    ! Abs. Xsection in band 7                              
      REAL BND8    ! Abs. Xsection in band 8                              
      REAL BND9    ! Abs. Xsection in band 9                              
                                                                          
      DATA BND1, BND2 / 5.47, 0.83 /                                      
      DATA BND3, BND4 / 2.03, 0.13 /                                      
      DATA BND5, BND6 / 1.33, 3.99 /                                      
      DATA BND7, BND8, BND9 / 8.61, 0.65, 1.86 /                          
                                                                          
C--- Band 1                                                               
      IF (IBND.EQ.1) THEN                                                 
         TRCNT=EXP(-WVEFF*BND1)                                           
      END IF                                                              
                                                                          
C--- Band 2                                                               
      IF (IBND.EQ.2) THEN                                                 
         TRCNT=EXP(-WVEFF*BND2)                                           
      END IF                                                              
                                                                          
C--- Band 3                                                               
      IF (IBND.EQ.3) THEN                                                 
         TRCNT=EXP(-WVEFF*BND3)                                           
      END IF                                                              
                                                                          
C--- Band 4                                                               
      IF (IBND.EQ.4) THEN                                                 
         TRCNT=EXP(-WVEFF*BND4)                                           
      END IF                                                              
                                                                          
C--- Band 5                                                               
      IF (IBND.EQ.5) THEN                                                 
         TRCNT=EXP(-WVEFF*BND5)                                           
      END IF                                                              
                                                                          
C--- Band 6                                                               
      IF (IBND.EQ.6) THEN                                                 
         TRCNT=EXP(-WVEFF*BND6)                                           
      END IF                                                              
                                                                          
C--- Band 7                                                               
      IF (IBND.EQ.7) THEN                                                 
         TRCNT=EXP(-WVEFF*BND7)                                           
      END IF                                                              
                                                                          
C--- Band 8                                                               
      IF (IBND.EQ.8) THEN                                                 
         TRCNT=EXP(-WVEFF*BND8)                                           
      END IF                                                              
C--- Band 9                                                               
      IF (IBND.EQ.9) THEN                                                 
         TRCNT=EXP(-WVEFF*BND9)                                           
      END IF                                                              
                                                                          
      END                                                                 
                                                                          
                                                                          
C***********************************************************************  
C*                                                                     *  
C*                         SUBROUTINE SMFLUX                           *  
C*                                                                     *  
C***********************************************************************  
                                                                          
      SUBROUTINE SMFLUX(NLEV,PFLUX,FLUX)                                  
                                                                          
C Subroutine SMFLUX smooths out the flux profiles, assuming that          
C the fluxes in an atmospheric layers of a certain thickness above and    
C below an atmospheric level vary as   FLUX = A0 + A1*P + A2*(P**2)       
C                                                                         
C INPUT:  Number of levels in the atmosphere without the surface (NLEV)   
C         Pressure at each level (PFLUX)                                  
C         Irradiance at each mid-level (FUP)                              
C                                                                         
C OUTPUT: Smoothed upward and downward irradiances at                     
C         each level (FUP,FDWN)                                           
                                                                          
      IMPLICIT NONE                                                       
C***********************************************************************  
C*                                                                     *  
C*                         P A R R A Y                                 *  
C*                                                                     *  
C***********************************************************************  
C-----------------------------------------------------------------------  
C In this part of the code the user can set the dimensions of arrays      
C used for the calculations, by setting the values in the paramter        
C statements                                                              
C-----------------------------------------------------------------------  

      INTEGER NN,MM,NHEM,NL,MOCT,MG,JG,NWJ2,NCRAY,JGL,NTRAC,NLEVRF            
      include 'params.i'
                                                    
      INTEGER MXGAS,MXLEV,MXBAND,MXCL                                     
                                                                          
      PARAMETER(MXGAS=8)     ! Maximum number of gases                    
                                                                          
      PARAMETER(MXLEV=NL)     ! Maximum number of levels in the            
                             ! atmosphere (not including the surface)     
                                                                          
      PARAMETER(MXBAND=9)    ! Maximum number of spectral bands (not      
                             ! including the whole spectrum, 0-3000cm-1)  
                                                                          
      PARAMETER(MXCL=3)      ! Maximum number of cloud types              
                                                                          
C-----------------------------------------------------------------------  
                                                                          
C-----------------                                                        
C Input Variables                                                         
C-----------------                                                        
                                                                          
      INTEGER NLEV                                                        
      REAL PFLUX(0:MXLEV),FLUX(0:MXLEV)                                   
                                                                          
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER ILAY,ILEV,KOUNTSM                                           
      REAL THICKN    ! Thickness of a layer above/below a level in mb     
      REAL DET,A0(0:MXLEV),A1(0:MXLEV),A2(0:MXLEV)                        
      INTEGER SN                               ! Sums used for least      
      REAL SX,SX2,SX3,SX4,SX2Y,SXY,SY        ! square fit method          
                                                                          
                                                                          
C--- Set the thickness of the layers in which regression is performed     
                                                                          
      IF (NLEV.GT.30) THEN                                                
         THICKN=75.0                                                      
      ELSE                                                                
         THICKN=100.0                                                     
      ENDIF                                                               
                                                                          
C--- Calculate regression coefficients for the first atmospheric level    
                                                                          
      SN=0                                                                
      SX=0.0                                                              
      SX2=0.0                                                             
      SX3=0.0                                                             
      SX4=0.0                                                             
      SX2Y=0.0                                                            
      SXY=0.0                                                             
      SY=0.0                                                              
                                                                          
      DO ILAY=0,NLEV                                                      
         IF ((PFLUX(0)-PFLUX(ILAY)).LE.100.0) THEN                        
            SN=SN+1                                                       
            SX=SX+PFLUX(ILAY)                                             
            SX2=SX2+(PFLUX(ILAY)**2)                                      
            SX3=SX3+(PFLUX(ILAY)**3)                                      
            SX4=SX4+(PFLUX(ILAY)**4)                                      
            SX2Y=SX2Y+((PFLUX(ILAY)**2)*FLUX(ILAY))                       
            SXY=SXY+(PFLUX(ILAY)*FLUX(ILAY))                              
            SY=SY+FLUX(ILAY)                                              
         END IF                                                           
      END DO                                                              
                                                                          
      DET=(SX2**3)+(SX*SX*SX4)+(SN*SX3*SX3)-                              
     $     SX2*(SN*SX4+2*SX*SX3)                                          
                                                                          
      IF ((DET.EQ.0.0).OR.(SN.LE.2)) THEN                                 
         SN=0                                                             
         SX=0.0                                                           
         SX2=0.0                                                          
         SX3=0.0                                                          
         SX4=0.0                                                          
         SX2Y=0.0                                                         
         SXY=0.0                                                          
         SY=0.0                                                           
         DO ILAY=0,2                                                      
            SN=SN+1                                                       
            SX=SX+PFLUX(ILAY)                                             
            SX2=SX2+(PFLUX(ILAY)**2)                                      
            SX3=SX3+(PFLUX(ILAY)**3)                                      
            SX4=SX4+(PFLUX(ILAY)**4)                                      
            SX2Y=SX2Y+((PFLUX(ILAY)**2)*FLUX(ILAY))                       
            SXY=SXY+(PFLUX(ILAY)*FLUX(ILAY))                              
            SY=SY+FLUX(ILAY)                                              
         END DO                                                           
         DET=(SX2**3)+(SX*SX*SX4)+(SN*SX3*SX3)-                           
     $        SX2*(SN*SX4+2*SX*SX3)                                       
      END IF                                                              
                                                                          
      A0(0)=(SX2*SX2*SX2Y)+SX4*(SX*SXY-SX2*SY)+                           
     $       SX3*(SY*SX3-SXY*SX2-SX2Y*SX)                                 
                                                                          
      A1(0)=(SX*SY*SX4)+SN*(SX3*SX2Y-SX4*SXY)+                            
     $       SX2*(SXY*SX2-SX*SX2Y-SX3*SY)                                 
                                                                          
      A2(0)=(SY*SX2*SX2)+SN*(SXY*SX3-SX2*SX2Y)+                           
     $       SX*(SX*SX2Y-SX3*SY-SX2*SXY)                                  
                                                                          
      A0(0)=A0(0)/DET                                                     
      A1(0)=A1(0)/DET                                                     
      A2(0)=A2(0)/DET                                                     
                                                                          
C--- Calculate regression coefficients for the upper atmospheric level    
                                                                          
      SN=0.0                                                              
      SX=0.0                                                              
      SX2=0.0                                                             
      SX3=0.0                                                             
      SX4=0.0                                                             
      SX2Y=0.0                                                            
      SXY=0.0                                                             
      SY=0.0                                                              
                                                                          
      DO ILAY=1,NLEV                                                      
         IF ((PFLUX(ILAY)-PFLUX(NLEV)).LE.200.0) THEN                     
            SN=SN+1                                                       
            SX=SX+PFLUX(ILAY)                                             
            SX2=SX2+(PFLUX(ILAY)**2)                                      
            SX3=SX3+(PFLUX(ILAY)**3)                                      
            SX4=SX4+(PFLUX(ILAY)**4)                                      
            SX2Y=SX2Y+((PFLUX(ILAY)**2)*FLUX(ILAY))                       
            SXY=SXY+(PFLUX(ILAY)*FLUX(ILAY))                              
            SY=SY+FLUX(ILAY)                                              
         END IF                                                           
      END DO                                                              
                                                                          
      DET=(SX2**3)+(SX*SX*SX4)+(SN*SX3*SX3)-                              
     $     SX2*(SN*SX4+2*SX*SX3)                                          
      IF ((DET.EQ.0.0).OR.(SN.LE.2)) THEN                                 
         SN=0.0                                                           
         SX=0.0                                                           
         SX2=0.0                                                          
         SX3=0.0                                                          
         SX4=0.0                                                          
         SX2Y=0.0                                                         
         SXY=0.0                                                          
         SY=0.0                                                           
         DO ILAY=NLEV-2,NLEV                                              
            SN=SN+1                                                       
            SX=SX+PFLUX(ILAY)                                             
            SX2=SX2+(PFLUX(ILAY)**2)                                      
            SX3=SX3+(PFLUX(ILAY)**3)                                      
            SX4=SX4+(PFLUX(ILAY)**4)                                      
            SX2Y=SX2Y+((PFLUX(ILAY)**2)*FLUX(ILAY))                       
            SXY=SXY+(PFLUX(ILAY)*FLUX(ILAY))                              
            SY=SY+FLUX(ILAY)                                              
         END DO                                                           
         DET=(SX2**3)+(SX*SX*SX4)+(SN*SX3*SX3)-                           
     $        SX2*(SN*SX4+2*SX*SX3)                                       
      END IF                                                              
                                                                          
      A0(NLEV)=(SX2*SX2*SX2Y)+SX4*(SX*SXY-SX2*SY)+                        
     $          SX3*(SY*SX3-SXY*SX2-SX2Y*SX)                              
                                                                          
      A1(NLEV)=(SX*SY*SX4)+SN*(SX3*SX2Y-SX4*SXY)+                         
     $          SX2*(SXY*SX2-SX*SX2Y-SX3*SY)                              
                                                                          
      A2(NLEV)=(SY*SX2*SX2)+SN*(SXY*SX3-SX2*SX2Y)+                        
     $          SX*(SX*SX2Y-SX3*SY-SX2*SXY)                               
                                                                          
      A0(NLEV)=A0(NLEV)/DET                                               
      A1(NLEV)=A1(NLEV)/DET                                               
      A2(NLEV)=A2(NLEV)/DET                                               
                                                                          
C--- Calculate regression coefficients for the rest of the levels         
                                                                          
      DO ILAY=1,NLEV-1                                                    
                                                                          
         KOUNTSM=0.                                                       
                                                                          
 1997    CONTINUE                                                         
                                                                          
         SN=0                                                             
         SX=0.0                                                           
         SX2=0.0                                                          
         SX3=0.0                                                          
         SX4=0.0                                                          
         SX2Y=0.0                                                         
         SXY=0.0                                                          
         SY=0.0                                                           
                                                                          
         DO ILEV=0,NLEV                                                   
            IF (ABS(PFLUX(ILAY)-PFLUX(ILEV)).LE.THICKN) THEN              
               SN=SN+1                                                    
               SX=SX+PFLUX(ILEV)                                          
               SX2=SX2+(PFLUX(ILEV)**2)                                   
               SX3=SX3+(PFLUX(ILEV)**3)                                   
               SX4=SX4+(PFLUX(ILEV)**4)                                   
               SX2Y=SX2Y+((PFLUX(ILEV)**2)*FLUX(ILEV))                    
               SXY=SXY+(PFLUX(ILEV)*FLUX(ILEV))                           
               SY=SY+FLUX(ILEV)                                           
            END IF                                                        
         END DO                                                           
                                                                          
         DET=(SX2**3)+(SX*SX*SX4)+(SN*SX3*SX3)-                           
     $        SX2*(SN*SX4+2*SX*SX3)                                       
                                                                          
         IF ((DET.EQ.0.0).OR.(SN.LE.2)) THEN                              
            KOUNTSM=KOUNTSM+1                                             
            IF (KOUNTSM.GT.5)THEN                                         
               WRITE(2,*)'SUBROUTINE SMFLUX STUCK IN INFINITE             
     $                    LOOP. ABORT.'                                   
               CALL ABORT                                                 
            ENDIF                                                         
            THICKN=THICKN+100.0                                           
            GOTO 1997                                                     
         END IF                                                           
                                                                          
         A0(ILAY)=(SX2*SX2*SX2Y)+SX4*(SX*SXY-SX2*SY)+                     
     $             SX3*(SY*SX3-SXY*SX2-SX2Y*SX)                           
                                                                          
         A1(ILAY)=(SX*SY*SX4)+SN*(SX3*SX2Y-SX4*SXY)+                      
     $             SX2*(SXY*SX2-SX*SX2Y-SX3*SY)                           
                                                                          
         A2(ILAY)=(SY*SX2*SX2)+SN*(SXY*SX3-SX2*SX2Y)+                     
     $             SX*(SX*SX2Y-SX3*SY-SX2*SXY)                            
                                                                          
         A0(ILAY)=A0(ILAY)/DET                                            
         A1(ILAY)=A1(ILAY)/DET                                            
         A2(ILAY)=A2(ILAY)/DET                                            
                                                                          
      END DO                                                              
                                                                          
      DO ILAY=0,NLEV                                                      
         FLUX(ILAY)=A0(ILAY)+PFLUX(ILAY)*(A1(ILAY)+                       
     $              A2(ILAY)*PFLUX(ILAY))                                 
         IF (FLUX(ILAY).LT.0.0) THEN                                      
            FLUX(ILAY)=0.0                                                
         END IF                                                           
      END DO                                                              
                                                                          
                                                                          
      END                                                                 
