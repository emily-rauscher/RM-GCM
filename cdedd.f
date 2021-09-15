C**********************************************************               
C             SUBROUTINE DEDD                                             
C**********************************************************               
      SUBROUTINE DEDD(PGG,PREF,PRMUZ,PTO1,PW,PRE1,PTR1,PRE2,PTR2)         
C                                                                         
C**** *DEDD* - COMPUTES REFLECTIVITY, TRANSMISSIVITY OF A CLOUDY LAYER    
C                                                                         
C     PURPOSE.                                                            
C     --------                                                            
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY      
C     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.                    
C                                                                         
C**   INTERFACE.                                                          
C     ----------                                                          
C          *DEDD* IS CALLED BY *SW*.                                      
C                                                                         
C                                                                         
C        IMPLICIT ARGUMENTS                                               
C        --------------------                                             
C                                                                         
C     ==== INPUTS ===                                                     
C PGG    : (NDLON)             ; ASSYMETRY FACTOR                         
C PREF   : (NDLON)             ; REFLECTIVITY OF THE UNDERLYING LAYER     
C PRMUZ  : (NDLON)             ; COSINE OF SOLAR ZENITH ANGLE             
C PTO1   : (NDLON)             ; OPTICAL THICKNESS                        
C PW     : (NDLON)             ; SINGLE SCATTERING ALBEDO                 
C     ==== OUTPUTS ===                                                    
C PRE1   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING NO           
C                              ; REFLECTION FROM UNDERLYING LAYER         
C PTR1   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING NO         
C                              ; REFLECTION FROM UNDERLYING LAYER         
C PRE2   : (NDLON)             ; LAYER REFLECTIVITY ASSUMING              
C                              ; REFLECTION FROM UNDERLYING LAYER         
C PTR2   : (NDLON)             ; LAYER TRANSMISSIVITY ASSUMING            
C                              ; REFLECTION FROM UNDERLYING LAYER         
C                                                                         
C     METHOD.                                                             
C     -------                                                             
C                                                                         
C          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.                   
C                                                                         
C PMF 23-2-98                                                             
      parameter(NLON=1)  !assumes 1  longitude only                       
C     -------------------------------------------------------------       
C                                                                         
      REAL PGG(NLON),PREF(NLON),PRE1(NLON),PRE2(NLON)                     
     : ,PRMUZ(NLON),PTO1(NLON)                                            
     :  ,PTR1(NLON),PTR2(NLON),PW(NLON)                                   
C                                                                         
C*         1.      DELTA-EDDINGTON CALCULATIONS                           
C                                                                         
C                                                                         
      DO 131 JLON =  1 , NLON                                             
C                                                                         
C*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS                   
C                                                                         
C                                                                         
         ZFF = PGG(JLON)*PGG(JLON)                                        
         ZGP = PGG(JLON)/(1.+PGG(JLON))                                   
         ZTOP = (1.- PW(JLON) * ZFF) * PTO1(JLON)                         
         ZWCP = (1-ZFF)* PW(JLON) /(1.- PW(JLON) * ZFF)                   
         ZDT = 2./3.                                                      
         ZX1 = 1.-ZWCP*ZGP                                                
         ZWM = 1.-ZWCP                                                    
         ZRM2 =  PRMUZ(JLON) * PRMUZ(JLON)                                
         ZZRK = SQRT(3.*ZWM*ZX1)                                          
         ZX2 = 4.*(1.-ZZRK*ZZRK*ZRM2)                                     
         ZRP = SQRT(3.*ZWM/ZX1)                                           
         ZALPHA = 3.*ZWCP*ZRM2*(1.+ZGP*ZWM)/ZX2                           
         ZBETA = 3.*ZWCP* PRMUZ(JLON) *(1.+3.*ZGP*ZRM2*ZWM)/ZX2           
         ZEXMUO = EXP(-ZTOP/ PRMUZ(JLON) )                                
         ZEXKP = EXP(ZZRK*ZTOP)                                           
         ZEXKM = 1./ZEXKP                                                 
         ZXP2P = 1.+ZDT*ZRP                                               
         ZXM2P = 1.-ZDT*ZRP                                               
         ZAP2B = ZALPHA+ZDT*ZBETA                                         
         ZAM2B = ZALPHA-ZDT*ZBETA                                         
C                                                                         
C*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER           
C                                                                         
 120  CONTINUE                                                            
C                                                                         
         ZA11 = ZXP2P                                                     
         ZA12 = ZXM2P                                                     
         ZA13 = ZAP2B                                                     
         ZA22 = ZXP2P*ZEXKP                                               
         ZA21 = ZXM2P*ZEXKM                                               
         ZA23 = ZAM2B*ZEXMUO                                              
         ZDENA = ZA11 * ZA22 - ZA21 * ZA12                                
         ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA                               
         ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA                               
         ZRI0A = ZC1A+ZC2A-ZALPHA                                         
         ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA                                    
         PRE1(JLON) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JLON)                      
         ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMUO                      
         ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMUO                 
         PTR1(JLON) = ZEXMUO+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JLON)               
C                                                                         
C*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER              
C                                                                         
 130  CONTINUE                                                            
C                                                                         
         ZB21 = ZA21- PREF(JLON) *ZXP2P*ZEXKM                             
         ZB22 = ZA22- PREF(JLON) *ZXM2P*ZEXKP                             
         ZB23 = ZA23- PREF(JLON) *ZEXMUO*(ZAP2B - PRMUZ(JLON) )           
         ZDENB = ZA11 * ZB22 - ZB21 * ZA12                                
         ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB                               
         ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB                               
         ZRI0C = ZC1B+ZC2B-ZALPHA                                         
         ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA                                    
         PRE2(JLON) = (ZRI0C-ZDT*ZRI1C) / PRMUZ(JLON)                     
         ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMUO                  
         ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMUO           
         PTR2(JLON) = ZEXMUO + (ZRI0D + ZDT*ZRI1D) / PRMUZ(JLON)          
C                                                                         
 131  CONTINUE                                                            
      RETURN                                                              
      END                                                                 
