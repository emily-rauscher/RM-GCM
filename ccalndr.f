C*********************************************************************    
C*                     SUBROUTINE CALNDR                             *    
C*********************************************************************    
C                                                                         
C--------------------------------------------------------------------     
C*  Given the julian day, routine returns month and day as                
C*    character*7 variable. Leap years are ignored.                       
C*                                                                        
C---------------------------------------------------------------------    
C                                                                         
C   Declarations.                                                         
C   -------------                                                         
      SUBROUTINE CALNDR(DOY,MDUM,amfrac)                                  
                                                                          
      IMPLICIT NONE                                                       
                                                                          
C Input day number                                                        
                                                                          
      INTEGER JDAY                                                        
      REAL DOY,DAY,GAP,AMFRAC                                             
      INTEGER MDUM                                                        
                                                                          
      CHARACTER*4 MONTH                                                   
      CHARACTER*4 MNTH(12)                                                
                                                                          
      INTEGER MM, MDAY                                                    
      INTEGER NDAY(13)                                                    
                                                                          
C      DATA NDAY/0,31,59,90,120,151,181,212,243,273,304,334,365/          
      DATA NDAY/0,30,60,90,120,150,180,210,240,270,300,330,360/           
     :    ,MNTH/'JAN.','FEB.','MAR.','APR.','MAY ','JUNE','JULY',         
     :           'AUG.','SEPT','OCT.','NOV.','DEC'/                       
                                                                          
C   Perform operation.                                                    
C   ------------------                                                    
                                                                          
       DAY=DOY-1.0                                                        
      DO 100 MM=12,1,-1                                                   
          IF (DAY.GE.(REAL(NDAY(MM+1)+NDAY(MM))/2.0)) THEN                
           MDUM=MM                                                        
           GAP=30.0  ! const month length                                 
           AMFRAC=(DAY-0.5*REAL(NDAY(MM+1)+NDAY(MM)))/GAP                 
           GOTO 200                                                       
          ENDIF                                                           
  100 CONTINUE                                                            
c this means day is between 1-15                                          
      MDUM=12                                                             
      GAP=30.0                                                            
      AMFRAC=(DAY+15.)/GAP                                                
                                                                          
  200 CONTINUE                                                            
                                                                          
      END                                                                 
