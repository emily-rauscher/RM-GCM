C**********************************************************               
C             SUBROUTINE MATINV                                           
C**********************************************************               
      SUBROUTINE MATINV(A,N,LDA,IWORK,WORK)                               
C                                                                         
C     This subroutine calculates the inverse of NxN array A.              
C                                                                         
C     Arguements:                                                         
C                                                                         
C       A     - Array of dimension (LDA,N). Contains inverse              
C               of A on exit.                                             
C       N     - Number of rows and columns of A.                          
C       LDA   - Leading dimension of A.                                   
C       IWORK - Integer array contains the pivot indices of A.            
C       WORK  - Real array used as workspace for SGETRI.                  
C                                                                         
C     This subroutine replaces the Cray specific subroutine               
C     MINV with portable LAPACK subroutines SGETRF and SGETRI.            
C                                                                         
      INTEGER N,LDA,IWORK(N),INFO                                         
      REAL A(LDA,N),WORK(N)                                               
C                                                                         
      CALL SGETRF(N,N,A,LDA,IWORK,INFO)                                   
      IF (INFO.NE.0) THEN                                                 
         WRITE(*,*) 'Error: SGETRF returned INFO = ',INFO                 
         CALL ABORT                                                       
      ENDIF                                                               
      CALL SGETRI(N,A,LDA,IWORK,WORK,N,INFO)                              
      IF (INFO.NE.0) THEN                                                 
         WRITE(*,*) 'Error: SGETRI returned INFO = ',INFO                 
         CALL ABORT                                                       
      ENDIF                                                               
C                                                                         
      END                                                                 
