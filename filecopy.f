C**********************************************************               
C             SUBROUTINE FILECOPY 
C**********************************************************               
      SUBROUTINE FILECOPY(ITS,IFT,ISF)

      REAL TLON,TLAT,TU,TV,TT,TSP,TF,DAY,SSLON
      INTEGER JG2,MG,NL,LL


      REWIND 26 
      REWIND 50
      REWIND 64

      READ(26,106) JG2,MG,NL
      WRITE(ITS,106) JG2,MG,NL
      READ(50,101) JG2,MG
      WRITE(IFT,101) JG2,MG
      READ(64,106) JG2,MG,NL
      WRITE(ISF,106) JG2,MG,NL
 106  FORMAT(3I5)     
 101  FORMAT(2I5)

      DO 20 L=1,NL
         DO 21 I=1,MG 
            DO 22 J=1,JG2
               READ(26,107) TLON,TLAT,LL,TU,TV,TT
               WRITE(ITS,107) TLON,TLAT,LL,TU,TV,TT
               IF (L.EQ.NL) THEN
                  READ(50,102) TLON,TLAT,TSP
                  WRITE(IFT,102) TLON,TLAT,TSP
                  READ(64,102) TLON,TLAT,TF
                  WRITE(ISF,102) TLON,TLAT,TF
               ENDIF

 22         CONTINUE
 21      CONTINUE
 20   CONTINUE



 107  FORMAT(2E13.5,I4,2E13.5)
 102  FORMAT(3E13.5)

      READ(26,105) DAY,SSLON,SSLAT
      WRITE(ITS,105) DAY,SSLON,SSLAT
      READ(50,105) DAY,SSLON,SSLAT
      WRITE(IFT,105) DAY,SSLON,SSLAT
      READ(64,105) DAY,SSLON,SSLAT
      WRITE(ISF,105) DAY,SSLON,SSLAT

 105  FORMAT(/' OUTPUTS FOR DAY ',F10.4,', SUBSTELLAR LON, LAT:',2F8.3)

      CLOSE(ITS)
      CLOSE(IFT)
      CLOSE(ISF)

      RETURN
      END                                                                 
                                                                          
