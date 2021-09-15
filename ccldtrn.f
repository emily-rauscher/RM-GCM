C***********************************************************************  
C*                         SUBROUTINE CLDTRAN                          *  
C***********************************************************************  
                                                                          
      SUBROUTINE CLDTRAN(NLEV,PFLUX,DP,CL,ACL,LWCCL,BASCL,TOPCL,          
     $                   TRCL,MTRCL)                                      
                                                                          
C Subroutine CLDTRAN is the cloud parametrization of the WBM. It is       
C based on the Chou & Suarez (1994) scheme, which extended the method     
C of Harshvardhan et al. (1987), written for black clouds, for the case   
C of grey clouds. The scheme approximates the cloud by an equivalent      
C effective overcast cloud and derives the cloud transmittance of         
C different atmospheric layers. The cloud emissivity follows the          
C parametrization of Ramanathan et al. (1983). The model can accomodate   
C up to 3 different cloud types, which are considered randomly            
C overlapped. No more than one cloud can be included in the same layer.   
C The cloud base and top must be on pressure levels of the input profile  
C                                                                         
C INPUT:  Number of levels in the atmosphere without the surface (NLEV)   
C         Pressure at flux levels (PFLUX)                                 
C         Thickness of the atmospheric layers (DP)                        
C         Selection of cloud types (CL)                                   
C         CL(1)low level, CL(2)mid level, CL(3)high level,                
C         CL(0) deep convective                                           
C         Cloud amount as a fraction (ACL)                                
C         Liquid water content (LWCCL)                                    
C         Cloud base and top (BASCL,TOPCL)                                
C                                                                         
C OUTPUT: Cloud transmittance (TRCL)                                      
C         Cloud transmittance for the modified transm. case (MTRCL)       
                                                                          
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
                                                                          
      INTEGER NLEV                                                        
      REAL PFLUX(0:MXLEV),DP(MXLEV)                                       
                                                                          
      INTEGER CL(0:MXCL)                                                  
      REAL ACL(0:MXCL),LWCCL(MXCL),BASCL(0:MXCL),TOPCL(0:MXCL)            
                                                                          
C------------------                                                       
C Output Variables                                                        
C------------------                                                       
                                                                          
      REAL TRCL(0:MXLEV-1,MXLEV),MTRCL(0:MXLEV-1,MXLEV)                   
                                                                          
C--------------------                                                     
C Internal Variables                                                      
C--------------------                                                     
                                                                          
      INTEGER ICL,ILEV,ILAY,ICOUNT,IDWN,IUP                               
      REAL DPCL            ! Cloud thickness                              
      REAL LWC(MXCL,MXLEV) ! LWC between 2 succesive levels               
      REAL LWPATH          ! LWC in atmospheric paths                     
      REAL EMSSCL(MXCL)    ! Cloud emissivity                             
      REAL CEMSS           ! Constant used in emissivity computations     
      PARAMETER (CEMSS=1.08942549)                                        
                                                                          
                                                                          
C Calculate the liquid water content between two succesive flux levels    
                                                                          
C--- Initialise LWC                                                       
      DO ILAY=1,NLEV                                                      
         DO ICL=1,MXCL                                                    
            LWC(ICL,ILAY)=0.0                                             
         END DO                                                           
      END DO                                                              
                                                                          
      DO ICL=1,MXCL                                                       
                                                                          
         ILAY=1                                                           
                                                                          
         IF (CL(ICL).EQ.1) THEN                                           
                                                                          
            DO WHILE ((100.*PFLUX(ILAY)).GT.BASCL(ICL))                   
               ILAY=ILAY+1                                                
            END DO                                                        
                                                                          
            DPCL=DP(ILAY)                                                 
            ICOUNT=ILAY                                                   
            DO WHILE ((100.*PFLUX(ICOUNT)).GT.TOPCL(ICL))                 
               ICOUNT=ICOUNT+1                                            
               DPCL=DPCL+DP(ICOUNT)                                       
            END DO                                                        
                                                                          
            LWC(ICL,ILAY)=(DP(ILAY)/DPCL)*0.5*LWCCL(ICL)                  
                                                                          
            DO WHILE ((100.*PFLUX(ILAY)).GT.TOPCL(ICL))                   
               ILAY=ILAY+1                                                
               LWC(ICL,ILAY)=(DP(ILAY)/DPCL)*LWCCL(ICL)                   
            END DO                                                        
                                                                          
         END IF                                                           
                                                                          
      END DO                                                              
                                                                          
                                                                          
C-----------------------------------------------------------------------  
C Calculate the cloud transmittance of the atmospheric paths considered   
C in the irradiance calculations                                          
C-----------------------------------------------------------------------  
                                                                          
C Initialise cloud transmittances                                         
      DO IUP=1,NLEV                                                       
         DO IDWN=0,IUP-1                                                  
            TRCL(IDWN,IUP)=1.0                                            
         END DO                                                           
      END DO                                                              
                                                                          
                                                                          
      DO IUP=1,NLEV                                                       
         DO IDWN=0,IUP-1                                                  
                                                                          
C Calculate LWC in the atmospheric path                                   
            DO ICL=1,MXCL                                                 
               LWPATH=0.                                                  
               IF (CL(ICL).EQ.1) THEN                                     
                  DO ILEV=IDWN+1,IUP                                      
                     LWPATH=LWPATH+LWC(ICL,ILEV)                          
                  END DO                                                  
               END IF                                                     
C Calculate cloud emissivity and transmittance                            
               IF (LWPATH.NE.0.) THEN                                     
                  IF (LWPATH.LE.25E-3) THEN                               
                     EMSSCL(ICL)=CEMSS*(1.-EXP(-100.*LWPATH))             
                     TRCL(IDWN,IUP)=TRCL(IDWN,IUP)*                       
     $                              (1.-ACL(ICL)*EMSSCL(ICL))             
                  ELSE                                                    
                     TRCL(IDWN,IUP)=TRCL(IDWN,IUP)*(1.-ACL(ICL))          
                  END IF                                                  
               END IF                                                     
            END DO                                                        
C Check for deep convective cloud                                         
            IF (CL(0).EQ.1) THEN                                          
               IF ((100.*PFLUX(IUP).GT.BASCL(0)).OR.                      
     $            (100.*PFLUX(IDWN).LT.TOPCL(0))) THEN                    
                  CONTINUE                                                
               ELSE                                                       
                  TRCL(IDWN,IUP)=(1.-ACL(0))*TRCL(IDWN,IUP)               
               END IF                                                     
            END IF                                                        
                                                                          
         END DO                                                           
      END DO                                                              
                                                                          
                                                                          
C-----------------------------------------------------------------------  
C Calculate cloud transmittances for the modified transmittance case      
C-----------------------------------------------------------------------  
                                                                          
C Initialise cloud transmittances                                         
      DO IUP=1,NLEV                                                       
         DO IDWN=0,NLEV-1                                                 
            MTRCL(IDWN,IUP)=1.0                                           
         END DO                                                           
      END DO                                                              
                                                                          
                                                                          
      DO IUP=1,NLEV                                                       
         DO IDWN=0,IUP-1                                                  
                                                                          
C Calculate LWC in the atmospheric path                                   
            DO ICL=1,MXCL                                                 
               LWPATH=0.                                                  
               IF (CL(ICL).EQ.1) THEN                                     
                  DO ILEV=IDWN+1,IUP                                      
                     IF (LWC(ICL,ILEV).NE.0.) THEN                        
                        IF (LWPATH.EQ.0.) THEN                            
                           LWPATH=0.5*LWC(ICL,ILEV)                       
                        ELSE                                              
                           LWPATH=LWPATH+LWC(ICL,ILEV)                    
                        END IF                                            
                     END IF                                               
                  END DO                                                  
               END IF                                                     
C Calculate cloud emissivity and transmittance                            
               IF (LWPATH.NE.0.) THEN                                     
                  IF (LWPATH.LE.25E-3) THEN                               
                     EMSSCL(ICL)=CEMSS*(1.-EXP(-100.*LWPATH))             
                     MTRCL(IDWN,IUP)=MTRCL(IDWN,IUP)*                     
     $                               (1.-ACL(ICL)*EMSSCL(ICL))            
                  ELSE                                                    
                     MTRCL(IDWN,IUP)=MTRCL(IDWN,IUP)*(1.-ACL(ICL))        
                  END IF                                                  
               END IF                                                     
            END DO                                                        
C Check for deep convective cloud                                         
            IF (CL(0).EQ.1) THEN                                          
               IF ((100.*PFLUX(IUP).GT.BASCL(0)).OR.                      
     $            (100.*PFLUX(IDWN).LT.TOPCL(0))) THEN                    
                  CONTINUE                                                
               ELSE                                                       
                  MTRCL(IDWN,IUP)=(1.-ACL(0))*MTRCL(IDWN,IUP)             
               END IF                                                     
            END IF                                                        
                                                                          
         END DO                                                           
      END DO                                                              
                                                                          
                                                                          
      END                                                                 
