C**********************************************************               
C             SUBROUTINE FINALORB
C**********************************************************               
      SUBROUTINE FINALORB(KOUNT,KTOTAL,TSPO,ITSOUT,IFTOUT,ISFOUT,ISWOUT)

      REAL FK,CHECK

      IF (TSPO.LE.90.) THEN
         CALL XSECT2
         CALL FILECOPY(ITSOUT,IFTOUT,ISFOUT,ISWOUT)
         ITSOUT=ITSOUT+1
         IFTOUT=IFTOUT+1
         ISFOUT=ISFOUT+1
         ISWOUT=ISWOUT+1
      ELSE
         FK=KOUNT-(KTOTAL-TSPO+1)
         CHECK=MOD(FK,TSPO/90.)
         IF (INT(CHECK).EQ.0) THEN
            CALL XSECT2
            CALL FILECOPY(ITSOUT,IFTOUT,ISFOUT,ISWOUT)
            ITSOUT=ITSOUT+1
            IFTOUT=IFTOUT+1
            ISFOUT=ISFOUT+1
            ISWOUT=ISWOUT+1
         ENDIF
      ENDIF

      RETURN
      END
