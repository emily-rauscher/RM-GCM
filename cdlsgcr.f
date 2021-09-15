C**********************************************************               
C             SUBROUTINE DLSGCR                                           
C**********************************************************               
      SUBROUTINE DLSGCR(NCOL,F,RGG,DF,IDIM,NL)                            
C                                                                         
C     Black-box routine to take vertical derivative of a 2D array         
C     in the model's vertical co-ordinates using the RGG matrix:          
C                                                                         
C              DF   =   RGG   *   F                                       
C                                                                         
C     gives    DF   =   d(F)/dln(sigma)   =   (sigma)d(F)/d(sigma).       
C                                                                         
C     Calculates DF at NL levels for NCOL columns from the array F.       
C     Uses vector loops over NCOL at each level.  Note that for           
C     NCOL <= 3 inner loops unroll to give vectorisation over levels.     
C                                                                         
C     Uses explicit values of RGG matrix rather than function SDOT        
C     over NL elements of RGG at each level.  Due to sparse nature        
C     of RGG, this gives large speed-up for NL >> 3.                      
C                                                                         
C     Reading model - sigma.                                              
C                                                                         
      REAL F(IDIM,NL),DF(IDIM,NL),RGG(NL,NL)                              
C                                                                         
 6900 FORMAT(/' ***ABORT IN DLSGCR: MORE COLUMNS THAN ARRAY DIMENSION:'   
     :       ,' NCOL IDIM ='2I10)                                         
C                                                                         
      IF (NCOL.GT.IDIM) THEN                                              
         WRITE(2,6900) NCOL,IDIM                                          
         CALL ABORT                                                       
      ENDIF                                                               
C                                                                         
      NLM=NL-1                                                            
      NLMM=NL-2                                                           
C     Top and bottom levels.                                              
      DO 10 I=1,NCOL                                                      
      DF(I,1 )=F(I,1)*RGG(1,1)+F(I,2)*RGG(2,1)+F(I,3)*RGG(3,1)            
      DF(I,NL)=F(I,NLMM)*RGG(NLMM,NL)+F(I,NLM)*RGG(NLM,NL)                
     :        +F(I,NL)*RGG(NL,NL)                                         
   10 CONTINUE                                                            
C     Intermediate levels.                                                
      DO 20 L=2,NLM                                                       
      DO 20 I=1,NCOL                                                      
      DF(I,L)=F(I,L-1)*RGG(L-1,L)+F(I,L+1)*RGG(L+1,L)                     
   20 CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
