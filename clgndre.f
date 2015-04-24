C**********************************************************               
C             SUBROUTINE LGNDRE                                           
C**********************************************************               
      SUBROUTINE LGNDRE(NN,MM,MOCT,ALP,DALP,MJP,JL,SIJ,CSJ)               
C                                                                         
C     Calculates legendre polynomials (ALP) and their derivatives         
C     (DALP) at the JL'th latitude (JL is in LEGAU) using                 
C     recurrence relationships.                                           
C                                                                         
      REAL ALP(MJP,JL),DALP(MJP,JL)                                       
C                                                                         
      LM=2                                                                
C                                                                         
C     Set P(0,0) and P(0,1)                                               
C                                                                         
      ALP(1,JL)=SQRT(.5)                                                  
      F1M=SQRT(1.5)                                                       
      ALP(2,JL)=F1M*SIJ                                                   
      DALP(1,JL)=0.                                                       
C                                                                         
C     Loop over wavenumbers                                               
C                                                                         
      DO 1 M1=1,MM                                                        
         M=M1-1                                                           
         AM=M                                                             
         A2M=M+M                                                          
         E2=SQRT(A2M+3.)                                                  
         IF (M.GT.0) THEN                                                 
            F2M=-F1M*CSJ/SQRT(A2M)                                        
            F1M=F2M*E2                                                    
            IF (M.NE.MMO) GOTO 1                                          
            LM=LM+1                                                       
            ALP(LM,JL)=F2M                                                
            LM=LM+1                                                       
            ALP(LM,JL)=F1M*SIJ                                            
            DALP(LM-1,JL)=-AM*ALP(LM,JL)/E2                               
         ENDIF                                                            
         M2=M+2                                                           
         MMO=M+MOCT                                                       
         JFM=((NN-M1)/2)*2+M2-1                                           
         IF (JFM.GE.M2) THEN                                              
            K=LM-M2+1                                                     
            AMSQ=AM*AM                                                    
C                                                                         
C           Loop over degree N                                            
C                                                                         
            DO 4 N=M2,JFM                                                 
               AN=N                                                       
               AN2=N*N                                                    
               ANM2=(N-1)*(N-1)                                           
               E1=SQRT((ANM2-AMSQ)/(4.*ANM2-1.))                          
               E2=SQRT((4.*AN2-1.)/(AN2-AMSQ))                            
               ALP(K+N,JL)=E2*(SIJ*ALP(K+N-1,JL)-E1*ALP(K+N-2,JL))        
               DALP(K+N-1,JL)=(1.-AN)*ALP(K+N,JL)/E2+AN*E1*ALP(K+N-2,JL)  
    4       CONTINUE                                                      
            LM=LM+JFM-M2+1                                                
         ENDIF                                                            
         DALP(LM,JL)=-AN*SIJ*ALP(LM,JL)+(AN+AN+1.0)*ALP(LM-1,JL)/E2       
    1 CONTINUE                                                            
      RETURN                                                              
      END                                                                 
