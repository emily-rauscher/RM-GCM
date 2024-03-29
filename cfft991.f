*DECK FFT991F
C     SUBROUTINE 'FFT991' - MULTIPLE FAST REAL PERIODIC TRANSFORM
C
C     Minor modifications - Clive Temperton, January 1991:
C         modified vector-chopping for better performance
C
C     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
C     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
C
C     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
C     WORK IS AN AREA OF SIZE (N+1)*MIN(LOT,64)
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS
C     IN PARALLEL
C
C     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
C
      NFAX=IFAX(1)
      NX=N+1
      IF (MOD(N,2).EQ.1) NX=N
      NBLOX=1+(LOT-1)/64
      LEFT=LOT
      IF (ISIGN.EQ.-1) GO TO 300
C
C     ISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
C     -----------------------------------------
  100 CONTINUE
      ISTART=1
      DO 220 NB=1,NBLOX
      IF (LEFT.LE.64) THEN
         NVEX=LEFT
      ELSE IF (LEFT.LT.128) THEN
         NVEX=LEFT/2
      ELSE
         NVEX=64
      ENDIF
      LEFT=LEFT-NVEX
      IA=ISTART
      I=ISTART
CDIR$ IVDEP
!DEC$ VECTOR ALWAYS
      DO 110 J=1,NVEX
      A(I+INC)=0.5*A(I)
      I=I+JUMP
  110 CONTINUE
      IF (MOD(N,2).EQ.1) GO TO 130
      I=ISTART+N*INC
      DO 120 J=1,NVEX
      A(I)=0.5*A(I)
      I=I+JUMP
  120 CONTINUE
  130 CONTINUE
      IA=ISTART+INC
      LA=1
      IGO=+1
C
      DO 160 K=1,NFAX
      IFAC=IFAX(K+1)
      IERR=-1
      IF (IGO.EQ.-1) GO TO 140
      CALL RPASSM(A(IA),A(IA+LA*INC),WORK(1),WORK(IFAC*LA+1),TRIGS,
     *    INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
      GO TO 150
  140 CONTINUE
      CALL RPASSM(WORK(1),WORK(LA+1),A(IA),A(IA+IFAC*LA*INC),TRIGS,
     *    1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
  150 CONTINUE
      IF (IERR.NE.0) GO TO 500
      LA=IFAC*LA
      IGO=-IGO
      IA=ISTART
  160 CONTINUE
C
C     IF NECESSARY, COPY RESULTS BACK TO A
C     ------------------------------------
      IF (MOD(NFAX,2).EQ.0) GO TO 190
      IBASE=1
      JBASE=IA
      DO 180 JJ=1,NVEX
      I=IBASE
      J=JBASE
      DO 170 II=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
  170 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  180 CONTINUE
  190 CONTINUE
C
C     FILL IN ZEROS AT END
C     --------------------
      IX=ISTART+N*INC
CDIR$ IVDEP
!DEC$ VECTOR ALWAYS
      DO 210 J=1,NVEX
      A(IX)=0.0
      A(IX+INC)=0.0
      IX=IX+JUMP
  210 CONTINUE
C
      ISTART=ISTART+NVEX*JUMP
  220 CONTINUE
      RETURN
C
C     ISIGN=-1, GRIDPOINT TO SPECTRAL TRANSFORM
C     -----------------------------------------
  300 CONTINUE
      ISTART=1
      DO 410 NB=1,NBLOX
      IF (LEFT.LE.64) THEN
         NVEX=LEFT
      ELSE IF (LEFT.LT.128) THEN
         NVEX=LEFT/2
      ELSE
         NVEX=64
      ENDIF
      LEFT=LEFT-NVEX
      IA=ISTART
      LA=N
      IGO=+1
C
      DO 340 K=1,NFAX
      IFAC=IFAX(NFAX+2-K)
      LA=LA/IFAC
      IERR=-1
      IF (IGO.EQ.-1) GO TO 320
      CALL QPASSM(A(IA),A(IA+IFAC*LA*INC),WORK(1),WORK(LA+1),TRIGS,
     *    INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
      GO TO 330
  320 CONTINUE
      CALL QPASSM(WORK(1),WORK(IFAC*LA+1),A(IA),A(IA+LA*INC),TRIGS,
     *    1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
  330 CONTINUE
      IF (IERR.NE.0) GO TO 500
      IGO=-IGO
      IA=ISTART+INC
  340 CONTINUE
C
C     IF NECESSARY, COPY RESULTS BACK TO A
C     ------------------------------------
      IF (MOD(NFAX,2).EQ.0) GO TO 370
      IBASE=1
      JBASE=IA
      DO 360 JJ=1,NVEX
      I=IBASE
      J=JBASE
      DO 350 II=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
  350 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  360 CONTINUE
  370 CONTINUE
C
C     SHIFT A(0) & FILL IN ZERO IMAG PARTS
C     ------------------------------------
      IX=ISTART
CDIR$ IVDEP
!DEC$ VECTOR ALWAYS
      DO 380 J=1,NVEX
      A(IX)=A(IX+INC)
      A(IX+INC)=0.0
      IX=IX+JUMP
  380 CONTINUE
      IF (MOD(N,2).EQ.1) GO TO 400
      IZ=ISTART+(N+1)*INC
      DO 390 J=1,NVEX
      A(IZ)=0.0
      IZ=IZ+JUMP
  390 CONTINUE
  400 CONTINUE
C
      ISTART=ISTART+NVEX*JUMP
  410 CONTINUE
      RETURN
C
C     ERROR MESSAGES
C     --------------
  500 CONTINUE
      GO TO (510,530,550) IERR
  510 CONTINUE
      WRITE(6,520) NVEX
  520 FORMAT(16H1VECTOR LENGTH =,I4,17H, GREATER THAN 64)
      GO TO 570
  530 CONTINUE
      WRITE(6,540) IFAC
  540 FORMAT( 9H1FACTOR =,I3,17H, NOT CATERED FOR)
      GO TO 570
  550 CONTINUE
      WRITE(6,560) IFAC
  560 FORMAT(9H1FACTOR =,I3,31H, ONLY CATERED FOR IF LA*IFAC=N)
  570 CONTINUE
      RETURN
      END
