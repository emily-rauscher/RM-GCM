        IF (mod(kount,ntstep).eq.1) then
          ilast=0

        ELSE                   ! ntstep requirement
          DO i=1,mg
            DO LD=1,NL
              im=i+IOFM
              TTRD(im,LD)=(htnet(ihem,jh,i,ld))/CHRF
            ENDDO
          ENDDO
        ENDIF
        IOFM=MGPP

        END DO

