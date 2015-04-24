C**********************************************************               
C             SUBROUTINE PQSAT                                            
C**********************************************************               
      real function PQSAT(T)                                              
       PARAMETER (NQSTAB=50000)                                           
       REAL PQSVAL                                                        
       COMMON /QSTABS/ PQSVAL(NQSTAB)                                     
      real T       
      PQSAT=PQSVAL(min(NQSTAB,(int(T*1e5))))
!      write(*,*) 'PQSVAL',PQSAT
      END                                                                 
