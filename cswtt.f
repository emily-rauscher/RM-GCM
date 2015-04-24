C**********************************************************               
C             SUBROUTINE SWTT                                             
C**********************************************************               
      SUBROUTINE SWTT (KNU,KA,APAD,BPAD,D,PU,PTR)                         
C                                                                         
C input everything apart from PTR (called ZR1 in SW and PU ZW)            
C     PURPOSE.                                                            
C     --------                                                            
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE  
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL  
C     INTERVALS.                                                          
C                                                                         
C**   INTERFACE.                                                          
C     ----------                                                          
C          *SWTT* IS CALLED FROM *SW*.                                    
C                                                                         
C        EXPLICIT ARGUMENTS                                               
C        --------------------                                             
C KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL           
C KA     :                     ; INDEX OF THE ABSORBER                    
C                                                                         
C        IMPLICIT ARGUMENTS                                               
C        --------------------                                             
C                                                                         
C     ==== INPUTS ===                                                     
C PU     : (KDLON)             ; ABSORBER AMOUNT                          
C     ==== OUTPUTS ===                                                    
C PTR    : (KDLON)             ; TRANSMISSION FUNCTION                    
C                                                                         
C     METHOD.                                                             
C     -------                                                             
C                                                                         
C          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS     
C     AND HORNER'S ALGORITHM.                                             
C                                                                         
C PMF 23-2-98                                                             
C-----------------------------------------------------------------        
      IMPLICIT NONE                                                       
      INTEGER NLON,JLON                                                   
      PARAMETER (NLON=1)                                                  
      INTEGER KNU,KA                                                      
      REAL PTR(NLON),PU(NLON)                                             
      REAL D(2,3),APAD(2,3,7),BPAD(2,3,7)                                 
      REAL ZR2(NLON)                                                      
C                                                                         
C*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION    
C                                                                         
      DO 101 JLON = 1 , NLON                                              
        PTR(JLON) = APAD(KNU,KA,1) + PU(JLON) *(APAD(KNU,KA,2)+PU(JLON)   
     S         * ( APAD(KNU,KA,3) + PU(JLON) * (APAD(KNU,KA,4)+PU(JLON)   
     S         * ( APAD(KNU,KA,5) + PU(JLON) * (APAD(KNU,KA,6)+PU(JLON)   
     S         * ( APAD(KNU,KA,7) ))))))                                  
C                                                                         
        ZR2(JLON) = BPAD(KNU,KA,1) + PU(JLON)*(BPAD(KNU,KA,2)+PU(JLON)    
     S         * ( BPAD(KNU,KA,3) + PU(JLON) *(BPAD(KNU,KA,4)+PU(JLON)    
     S         * ( BPAD(KNU,KA,5) + PU(JLON) *(BPAD(KNU,KA,6)+PU(JLON)    
     S         * ( BPAD(KNU,KA,7) ))))))                                  
 101  CONTINUE                                                            
C                                                                         
C*         2.      ADD THE BACKGROUND TRANSMISSION                        
C                                                                         
C                                                                         
      DO 201 JLON = 1 , NLON                                              
         PTR(JLON) = (PTR(JLON) / ZR2(JLON))*(1.-D(KNU,KA))+D(KNU,KA)     
 201  CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
