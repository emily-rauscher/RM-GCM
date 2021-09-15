C**********************************************************               
C             SUBROUTINE SPDEL2                                           
C**********************************************************               
      SUBROUTINE SPDEL2(Z,FILT,NWJ2,NN,MM,MOCT,NHEM,NL,IPAR,ITYPE)        
C                                                                         
C     Perform del**2 or del**(-2) operation on a spectral field.          
C     The input spectral array is assumed to use the jagged triangular    
C     truncation of the Reading baroclinic spectral models.  This is      
C     also used in the diagnostics program for fields from both Reading   
C     models and the UGCM.                                                
C                                                                         
C     Input arguments:                                                    
C        Z      - Complex array containing input spectral field.          
C                 The truncation is assumed to be jagged triangular.      
C                 i.e. symmetric (even) coefficients are included up      
C                 to total wavenumber (NN-1), while anti-symmetric        
C                 (odd) coefficients are included up to wavenumber NN.    
C                 This gives equal numbers of even and odd coefficients   
C                 in the truncated series.  Ordering is of increasing     
C                 total wavenumber within increasing zonal wavenumber.    
C        FILT   - Real array, unset, to receive filter coefficients.      
C        NWJ2   - First dimension of Z: number of even or odd coeffs      
C                 at a single level in the jagged triangular truncation.  
C        NN     - Highest total wavenumber of input truncation.           
C        MM     - Highest zonal wavenumber of input truncation.           
C        MOCT   - Symmetry in longitude.  Only zonal wavenumbers 0,MOCT,  
C                 2*MOCT,..,MM are included in the input truncation.      
C        NHEM   - Symmetry in latitude: { 1 = hemispheric, only even or   
C                                       {     odd coefficients included.  
C                                       { 2 = global, both even and odd   
C                                       {     coefficients included.      
C        NL     - Number of levels in vertical.                           
C        IPAR   - Parity of field: { IPAR=even for even hem symmetry,     
C                                  { IPAR=odd  for  odd hem symmetry.     
C        ITYPE  - Type of operation required: { +2, del**(+2),            
C                                             { -2, del**(-2).            
C     Output arguments:                                                   
C        Z      - Filtered spectral field.  Ordering of spectral          
C                 coefficients is unchanged.                              
C        FILT   - Real array of filter coefficients.                      
C        Other arguments unchanged.                                       
C                                                                         
C     Method:                                                             
C        This routine can perform the following operations:               
C     ITYPE=+2 : Del**(+2), in which the non-dimensional spectral         
C                coefficient Z(n,m) is multiplied by -[n*(n+1)].          
C     ITYPE=-2 : Del**(-2), in which the non-dimensional spectral         
C                coefficient Z(n,m) is divided by -[n*(n+1)].             
C                                                                         
C     Author:                                                             
C        Original version.                    Mike Blackburn,  25.11.94.  
C        FILT is now a dummy array, (0:NN).   Mike Blackburn,  04.09.96.  
C        FILT(0)=-1 for inverse Laplacian.    Mike Blackburn,  17.04.00.  
C                                                                         
      COMPLEX Z(NWJ2,NHEM,NL)                                             
      REAL FILT(0:NN)                                                     
C                                                                         
C     Set up filter coefficients.                                         
C     Note that FILT(n)=-n*(n+1) (or its inverse).                        
C                                                                         
      IF (ITYPE.EQ.2) THEN                                                
         DO 10 N=0,NN                                                     
            FILT(N)=-REAL(N*(N+1))                                        
   10    CONTINUE                                                         
      ELSE IF (ITYPE.EQ.-2) THEN                                          
         FILT(0)=-1.0                                                     
         DO 20 N=1,NN                                                     
            FILT(N)=-1./REAL(N*(N+1))                                     
   20    CONTINUE                                                         
      ELSE                                                                
         PRINT *,' ***SPDEL2: INVALID VALUE OF ITYPE SUPPLIED = '         
     :          ,ITYPE,' : MUST BE +-2 FOR DEL**(+-2) OPERATION'          
         CALL ABORT                                                       
      ENDIF                                                               
C                                                                         
C     Del**(+-2) operation.                                               
C     Counting for total wavenumber in inner loop is for even coeffs,     
C     NOF increases total wavenumber for odd coeffs.                      
C                                                                         
      DO 30 IHEM=1,NHEM                                                   
         NOF=1-MOD(IHEM+IPAR,2)                                           
         DO 30 L=1,NL                                                     
            I=0                                                           
            DO 30 M=0,MM-1,MOCT                                           
               DO 30 N=M,NN-1,2                                           
                  I=I+1                                                   
                  Z(I,IHEM,L)=Z(I,IHEM,L)*FILT(N+NOF)                     
   30 CONTINUE                                                            
C                                                                         
C     Check that inner loop count is NWJ2.                                
C                                                                         
      IF (I.NE.NWJ2) THEN                                                 
         PRINT *,' ***SPDEL2: ERROR IN WAVENUMBER COUNTING: FINAL I = '   
     :          ,I,' SHOULD EQUAL NWJ2 = ',NWJ2                           
         CALL ABORT                                                       
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
