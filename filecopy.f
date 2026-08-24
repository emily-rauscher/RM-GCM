C**********************************************************               
C             SUBROUTINE FILECOPY 
C**********************************************************               
      SUBROUTINE FILECOPY(ITS,IFT,ISF)

      REAL TLON,TLAT,TSP,TF,DAY,SSLON
      INTEGER JG2,MG,NL
C     The unit-26 record is relayed verbatim rather than parsed into a
C     fixed variable list, which used to drop any column added in
C     xsect2.f until this list was widened to match.
      CHARACTER*512 TLINE


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
               READ(26,'(A)') TLINE
               WRITE(ITS,'(A)') TLINE(1:LEN_TRIM(TLINE))
               IF (L.EQ.NL) THEN
                  READ(50,102) TLON,TLAT,TSP
                  WRITE(IFT,102) TLON,TLAT,TSP
                  READ(64,102) TLON,TLAT,TF
                  WRITE(ISF,102) TLON,TLAT,TF
               ENDIF

 22         CONTINUE
 21      CONTINUE
 20   CONTINUE



 102  FORMAT(3E13.5)

C     Format 105 (used for WRITE) contains character literals, which are not
C     legal in an INPUT format.  ifort accepts them and skips that many
C     columns; gfortran rejects the READ at runtime.  Format 106 is the same
C     layout with the literals replaced by equal-width X descriptors
C     (17 and 22 columns), so the READs consume exactly what ifort consumed.
      READ(26,107) DAY,SSLON,SSLAT
      WRITE(ITS,105) DAY,SSLON,SSLAT
      READ(50,107) DAY,SSLON,SSLAT
      WRITE(IFT,105) DAY,SSLON,SSLAT
      READ(64,107) DAY,SSLON,SSLAT
      WRITE(ISF,105) DAY,SSLON,SSLAT

 105  FORMAT(/' OUTPUTS FOR DAY ',F10.4,', SUBSTELLAR LON, LAT:',2F8.3)
 107  FORMAT(/17X,F10.4,22X,2F8.3)

      CLOSE(ITS)
      CLOSE(IFT)
      CLOSE(ISF)

      RETURN
      END                                                                 
                                                                          
