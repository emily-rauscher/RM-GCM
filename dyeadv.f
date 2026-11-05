C**********************************************************
C             MODULE DYEGRID
C**********************************************************
      MODULE dyegrid
C
C     Grid-point transport of the cloud dyes (QDYE in tracerclds): a
C     positive, limited, flux-form finite-volume scheme driven by the
C     spectral model's winds. The physics in DGRMLT then changes the
C     dyes in place.
C
C     Cells are the Gaussian grid points: MG longitude bins, and rows
C     1..2*JG (north to south) whose width in mu=sin(latitude) is the
C     Gaussian weight, so the poles are edges of zero length. Time is
C     the model's (1/Omega) and winds are UG=u cos(lat)/(a Omega),
C     VG=v cos(lat)/(a Omega). Cell air mass is M=p*/P0*DSIGMA*W (the
C     common dlambda dropped); a flux is the air mass crossing a face in
C     one step, positive toward the higher index (east, south, down).
C
C     DYEADV builds the face fluxes and applies them in HSTEP
C     (horizontal) and SWEEP1D (vertical), Easter (1993) split, with
C     face values from PPMFLUX.
C
      use tracerclds
      SAVE
C
C     Layer thicknesses in sigma, from DYE_START.
      REAL DSIGDY(NL)
C     Row edges in mu: EMUDY(R) separates rows R and R+1, EMUDY(0)=1.
      REAL EMUDY(0:2*JG)
C     Sweep-order switch, and the diagnostics accumulated between lines
C     of dye_adv.txt: largest subcycle counts (horizontal, vertical),
C     and the largest |fixer correction| per dye.
      INTEGER :: KDYSTP=0,IUDYA=0
      INTEGER NSUBMX(2)
      REAL CORRMX(2:NDYE+1)
C     Work arrays of DYEADV, kept here rather than on the stack: cell
C     air masses, face fluxes, and the winds used for the step.
      REAL AM(MG,2*JG,NL)
      REAL FX(MG,2*JG,NL),FY(MG,0:2*JG,NL),FZ(MG,2*JG,0:NL)
      REAL UE(MG,2*JG,NL),VE(MG,2*JG,NL)
C
      CONTAINS
C
C**********************************************************
      SUBROUTINE DYE_STOREROW(IROW,IOFM,IGCD,UG,VG,SPG,AMU,WGT)
C     Called by MGRMLT for one hemisphere of one latitude: keep the
C     row's winds and p*/P0 (=1+SPG) at time t, and its mu and weight.
      INTEGER IROW,IOFM,IGCD
      REAL UG(IGCD,NL),VG(IGCD,NL),SPG(IGCD),AMU,WGT
      INTEGER I,L
      IF (NDYE.LT.1) RETURN
      DO 20 L=1,NL
         DO 10 I=1,MG
            UDY(I,IROW,L)=UG(I+IOFM,L)
            VDY(I,IROW,L)=VG(I+IOFM,L)
   10    CONTINUE
   20 CONTINUE
      DO 30 I=1,MG
         PSDYE(I,IROW)=1.0+SPG(I+IOFM)
   30 CONTINUE
      AMUDY(IROW)=AMU
      WDY(IROW)=WGT
      END SUBROUTINE DYE_STOREROW
C
C**********************************************************
      SUBROUTINE DYE_START(LRST,DAY,SIGMA,DSIGMA,P0)
C     First call, once MGRMLT has stored every row: set up the grid and
C     the initial state, from fort.11.dye on a restart if there is one,
C     else the initial profile.
      LOGICAL LRST
      REAL DAY,SIGMA(NL),DSIGMA(NL),P0
      INTEGER I,R,L,KK
      REAL PR
      LOGICAL LEX
      IF (NDYE.LT.1 .OR. LDYEINI) RETURN
      DO 5 L=1,NL
         DSIGDY(L)=DSIGMA(L)
    5 CONTINUE
      EMUDY(0)=1.0
      DO 10 R=1,2*JG
         EMUDY(R)=EMUDY(R-1)-WDY(R)
   10 CONTINUE
      IF (ABS(EMUDY(2*JG)+1.0).GT.1.0E-10) THEN
         WRITE(*,*) 'DYE_START: Gaussian weights sum to ',
     &              1.0-EMUDY(2*JG),', not 2'
         STOP 'DYE_START: bad Gaussian weights'
      ENDIF
      DO 12 R=1,2*JG
         IF (AMUDY(R).GE.EMUDY(R-1) .OR. AMUDY(R).LE.EMUDY(R)) THEN
            WRITE(*,*) 'DYE_START: row ',R,' at mu=',AMUDY(R),
     &                 ' is outside its cell ',EMUDY(R),EMUDY(R-1)
            STOP 'DYE_START: rows out of order'
         ENDIF
   12 CONTINUE
      DAYDY=DAY
      DTDY=0.0
      LDYEPRV=.FALSE.
      LDYMTGT=.FALSE.
      DO 15 KK=2,NDYE+1
         DMPHYS(KK)=0.0
         CORRMX(KK)=0.0
   15 CONTINUE
      NSUBMX=0
      IF (LRST) THEN
         INQUIRE(FILE='fort.11.dye',EXIST=LEX)
         IF (LEX) THEN
            CALL DYE_RSREAD(DAY)
            LDYEINI=.TRUE.
            RETURN
         ENDIF
         WRITE(*,*) 'DYE_START: restart without fort.11.dye, so the ',
     &              'dyes start from the initial profile'
         WRITE(2,*) 'DYE_START: restart without fort.11.dye, so the ',
     &              'dyes start from the initial profile'
      ENDIF
      DO 40 KK=2,NDYE+1
         DO 30 L=1,NL
            DO 20 R=1,2*JG
               DO 20 I=1,MG
                  PR=SIGMA(L)*PSDYE(I,R)*P0
                  IF (KDYEMODE.EQ.2 .AND. MOD(KK,2).EQ.1) THEN
                     QDYE(I,R,L,KK)=0.0
                  ELSE IF (PR.GT.PDYEFIX .OR. PR.LT.PDYEUPPER) THEN
                     QDYE(I,R,L,KK)=1.0
                  ELSE
                     QDYE(I,R,L,KK)=0.0
                  ENDIF
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
      LDYEINI=.TRUE.
      END SUBROUTINE DYE_START
C
C**********************************************************
      SUBROUTINE DYE_RSWRITE(IUNIT)
C     Dye restart record fort.IUNIT.dye, written beside unit IUNIT (11 or
C     12) and read back as fort.11.dye; it holds everything the next
C     steps depend on, so a restart reproduces a straight run exactly.
      INTEGER IUNIT,IU
      CHARACTER*16 FNAME
      IF (NDYE.LT.1) RETURN
      WRITE(FNAME,'(A5,I2,A4)') 'fort.',IUNIT,'.dye'
      OPEN(NEWUNIT=IU,FILE=FNAME,FORM='UNFORMATTED',STATUS='REPLACE')
      WRITE(IU) NDYE,MG,JG,NL
      WRITE(IU) DAYDY,DTDYM,LDYEPRV,LDYMTGT,KDYSTP
      WRITE(IU) DYMTGT,DMPHYS
      WRITE(IU) UDYM,VDYM
      WRITE(IU) QDYE
      CLOSE(IU)
      END SUBROUTINE DYE_RSWRITE
C
C**********************************************************
      SUBROUTINE DYE_RSREAD(DAY)
      REAL DAY
      INTEGER IU,ND,MGR,JGR,NLR
      OPEN(NEWUNIT=IU,FILE='fort.11.dye',FORM='UNFORMATTED',
     &     STATUS='OLD')
      READ(IU) ND,MGR,JGR,NLR
      IF (ND.NE.NDYE .OR. MGR.NE.MG .OR. JGR.NE.JG .OR. NLR.NE.NL)
     &   THEN
         WRITE(*,*) 'DYE_RSREAD: fort.11.dye is NDYE,MG,JG,NL=',
     &              ND,MGR,JGR,NLR
         STOP 'DYE_RSREAD: fort.11.dye does not match params.i'
      ENDIF
      READ(IU) DAYDY,DTDYM,LDYEPRV,LDYMTGT,KDYSTP
      READ(IU) DYMTGT,DMPHYS
      READ(IU) UDYM,VDYM
      READ(IU) QDYE
      CLOSE(IU)
      IF (ABS(DAYDY-DAY).GT.1.0E-6*MAX(1.0,ABS(DAY))) THEN
         WRITE(*,*) 'DYE_RSREAD: fort.11.dye is for day ',DAYDY,
     &              ' but the restart is for day ',DAY
         STOP 'DYE_RSREAD: fort.11.dye is for another day'
      ENDIF
      DAYDY=DAY
      END SUBROUTINE DYE_RSREAD
C
C**********************************************************
      SUBROUTINE DYEADV(DAY,PI2,KOUNT,NDIAG)
C     Move the dyes from DAYDY to DAY. NDIAG: steps between lines of
C     dye_adv.txt.
      REAL DAY,PI2
      INTEGER KOUNT,NDIAG
      REAL DT,DLAM,PE,TI,WIN,SC,ACT,TGT,QMN,CR0,CR1,CRE
      REAL CDIV(NL),CSUM
      INTEGER I,IP,R,L,KK,ID,NSUB
      REAL ML(NL),QL(NL,NDYE),FL(0:NL)
C
      IF (NDYE.LT.1) RETURN
      DT=PI2*(DAY-DAYDY)
      IF (DT.LE.0.0) RETURN
      DLAM=PI2/REAL(MG)
C
C     Cell air masses from the spectral p*, then
C     the fixer (LDYEFIX), which rescales each dye to the mass the last
C     advection left plus the physics since (DMPHYS, from DGRMLT).
!$OMP PARALLEL DO PRIVATE(I,R)
      DO 10 L=1,NL
         DO 10 R=1,2*JG
            DO 10 I=1,MG
               AM(I,R,L)=PSDYE(I,R)*DSIGDY(L)*WDY(R)
   10 CONTINUE
!$OMP END PARALLEL DO
C
C     Serial sums, so the result does not depend on the thread count.
      DO 20 KK=2,NDYE+1
         ACT=0.0
         DO 15 L=1,NL
            DO 15 R=1,2*JG
               DO 15 I=1,MG
                  ACT=ACT+AM(I,R,L)*QDYE(I,R,L,KK)
   15    CONTINUE
         IF (LDYEFIX .AND. LDYMTGT .AND. ACT.GT.0.0) THEN
            TGT=DYMTGT(KK)+DMPHYS(KK)
            SC=MAX(TGT,0.0)/ACT
            CORRMX(KK)=MAX(CORRMX(KK),ABS(SC-1.0))
!$OMP PARALLEL DO PRIVATE(I,R)
            DO 18 L=1,NL
               DO 18 R=1,2*JG
                  DO 18 I=1,MG
                     QDYE(I,R,L,KK)=QDYE(I,R,L,KK)*SC
   18       CONTINUE
!$OMP END PARALLEL DO
         ENDIF
   20 CONTINUE
C
C     Winds extrapolated to the middle of the step (1.5*u(t)-0.5*u(t-1))
C     when the previous step had the same length.
      IF (LDYEPRV .AND. ABS(DT-DTDYM).LE.1.0E-9*DT) THEN
         WIN=1.5
      ELSE
         WIN=1.0
      ENDIF
!$OMP PARALLEL DO PRIVATE(I,R)
      DO 30 L=1,NL
         DO 30 R=1,2*JG
            DO 30 I=1,MG
               UE(I,R,L)=WIN*UDY(I,R,L)+(1.0-WIN)*UDYM(I,R,L)
               VE(I,R,L)=WIN*VDY(I,R,L)+(1.0-WIN)*VDYM(I,R,L)
   30 CONTINUE
!$OMP END PARALLEL DO
C
C     Horizontal face fluxes. FX(I,R,L): face between I and I+1 (MG
C     and 1 for I=MG). FY(I,R,L): face between rows R and R+1, with v
C     (not VG=v cos(lat), which is far from linear near the poles) and
C     p* linear in mu between the row centres; the polar faces are 0.
!$OMP PARALLEL DO PRIVATE(I,IP,R,PE,TI,CR0,CR1,CRE)
      DO 50 L=1,NL
         DO 40 R=1,2*JG
            DO 40 I=1,MG
               IP=MOD(I,MG)+1
               PE=0.5*(PSDYE(I,R)+PSDYE(IP,R))
               FX(I,R,L)=PE*DSIGDY(L)*WDY(R)*0.5*(UE(I,R,L)+UE(IP,R,L))
     &                   /(1.0-AMUDY(R)**2)*DT/DLAM
   40    CONTINUE
         DO 45 I=1,MG
            FY(I,0,L)=0.0
            FY(I,2*JG,L)=0.0
   45    CONTINUE
         DO 48 R=1,2*JG-1
            TI=(AMUDY(R)-EMUDY(R))/(AMUDY(R)-AMUDY(R+1))
            CR0=SQRT(1.0-AMUDY(R)**2)
            CR1=SQRT(1.0-AMUDY(R+1)**2)
            CRE=SQRT(1.0-EMUDY(R)**2)
            DO 48 I=1,MG
               PE=PSDYE(I,R)+TI*(PSDYE(I,R+1)-PSDYE(I,R))
               FY(I,R,L)=-PE*DSIGDY(L)*CRE*DT
     &                   *(VE(I,R,L)/CR0+TI*(VE(I,R+1,L)/CR1
     &                   -VE(I,R,L)/CR0))
   48    CONTINUE
   50 CONTINUE
!$OMP END PARALLEL DO
C
C     Vertical fluxes from the horizontal divergence: layer L keeps the
C     share DSIGMA(L) of the column's change, the rest passes on.
!$OMP PARALLEL DO PRIVATE(I,L,CDIV,CSUM)
      DO 70 R=1,2*JG
         DO 70 I=1,MG
            CSUM=0.0
            DO 60 L=1,NL
               CDIV(L)=FX(I,R,L)-FX(MOD(I+MG-2,MG)+1,R,L)
     &                +FY(I,R,L)-FY(I,R-1,L)
               CSUM=CSUM+CDIV(L)
   60       CONTINUE
            FZ(I,R,0)=0.0
            DO 65 L=1,NL-1
               FZ(I,R,L)=FZ(I,R,L-1)-CDIV(L)+DSIGDY(L)*CSUM
   65       CONTINUE
            FZ(I,R,NL)=0.0
   70 CONTINUE
!$OMP END PARALLEL DO
C
C     Horizontal and vertical parts, H then Z on odd steps and Z then
C     H on even ones.
      KDYSTP=KDYSTP+1
      DO 200 ID=1,2
         IF ((ID.EQ.1) .EQV. (MOD(KDYSTP,2).EQ.1)) THEN
!$OMP PARALLEL DO PRIVATE(NSUB) REDUCTION(MAX:NSUBMX)
            DO 110 L=1,NL
               CALL HSTEP(L,NSUB)
               NSUBMX(1)=MAX(NSUBMX(1),NSUB)
  110       CONTINUE
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO PRIVATE(I,L,KK,ML,QL,FL,NSUB) COLLAPSE(2)
!$OMP& REDUCTION(MAX:NSUBMX)
            DO 130 R=1,2*JG
               DO 130 I=1,MG
                  FL(0)=FZ(I,R,0)
                  DO 122 L=1,NL
                     ML(L)=AM(I,R,L)
                     FL(L)=FZ(I,R,L)
                     DO 122 KK=1,NDYE
                        QL(L,KK)=QDYE(I,R,L,KK+1)
  122             CONTINUE
                  CALL SWEEP1D(NL,ML,QL,FL,.FALSE.,NSUB)
                  NSUBMX(2)=MAX(NSUBMX(2),NSUB)
                  DO 124 L=1,NL
                     AM(I,R,L)=ML(L)
                     DO 124 KK=1,NDYE
                        QDYE(I,R,L,KK+1)=QL(L,KK)
  124             CONTINUE
  130       CONTINUE
!$OMP END PARALLEL DO
         ENDIF
  200 CONTINUE
C
C     Mass the advection leaves, the fixer's next target.
      DO 220 KK=2,NDYE+1
         ACT=0.0
         DO 215 L=1,NL
            DO 215 R=1,2*JG
               DO 215 I=1,MG
                  ACT=ACT+AM(I,R,L)*QDYE(I,R,L,KK)
  215    CONTINUE
         DYMTGT(KK)=ACT
         DMPHYS(KK)=0.0
  220 CONTINUE
      LDYMTGT=.TRUE.
C
      DO 230 L=1,NL
         DO 230 R=1,2*JG
            DO 230 I=1,MG
               UDYM(I,R,L)=UDY(I,R,L)
               VDYM(I,R,L)=VDY(I,R,L)
  230 CONTINUE
      LDYEPRV=.TRUE.
      DTDYM=DT
      DTDY=DT
      DAYDY=DAY
C
C     Diagnostics: per dye the lowest q now and the largest fixer
C     correction, and the largest subcycle counts, since the last line.
      IF (MOD(KOUNT,MAX(NDIAG,1)).EQ.0) THEN
         IF (IUDYA.EQ.0) THEN
            OPEN(NEWUNIT=IUDYA,FILE='dye_adv.txt',POSITION='APPEND')
            WRITE(IUDYA,'(A)') '# day, subcycles h z, then per dye'//
     &         ' 2..NDYE+1: min q, max |fixer correction|'
         ENDIF
         WRITE(IUDYA,'(F12.5,2I4)',ADVANCE='NO') DAY,NSUBMX
         DO 240 KK=2,NDYE+1
            QMN=MINVAL(QDYE(:,:,:,KK))
            WRITE(IUDYA,'(2(1X,ES11.4))',ADVANCE='NO') QMN,CORRMX(KK)
            CORRMX(KK)=0.0
  240    CONTINUE
         WRITE(IUDYA,*)
         FLUSH(IUDYA)
         NSUBMX=0
      ENDIF
      END SUBROUTINE DYEADV
C
C**********************************************************
      SUBROUTINE PPMFLUX(N,AM,Q,FS,LPER,G)
C     Dye mass G(0:N) carried through the faces of a line of N cells by
C     air-mass fluxes FS(0:N) (FS(K) between cells K and K+1, positive
C     toward K+1; periodic lines have FS(0)=FS(N) on the same face,
C     closed lines FS(0)=FS(N)=0), for cell air masses AM and mixing
C     ratios Q.
C
C     Each cell holds a PPM parabola in its index coordinate (Colella
C     and Woodward 1984) with the extremum-preserving limiter of Colella
C     and Sekora (2008, J. Comput. Phys. 227, 7069), so smooth peaks are
C     not clipped and jumps do not overshoot. A face carries the mean of
C     the upwind parabola over the slab that crosses it. A parabola that
C     dips below zero is replaced by its cell mean (Lin 2004), which
C     keeps q>=0 while no cell loses more than its air mass. Closed lines
C     extend their end values into two ghost cells.
      INTEGER N
      REAL AM(N),Q(N),FS(0:N),G(0:N)
      LOGICAL LPER
      REAL QE(-1:N+2),E(0:N),QL(N),QR(N),Q6(N)
      REAL C,QF,D,XM,PMN,D2,D2L,D2R,D2C,D2LM,SG,FAC
      REAL, PARAMETER :: CLIM=1.25
      INTEGER K
C     Line with two ghost cells each side.
      DO 10 K=1,N
         QE(K)=Q(K)
   10 CONTINUE
      IF (LPER) THEN
         QE(-1)=Q(N-1)
         QE(0)=Q(N)
         QE(N+1)=Q(1)
         QE(N+2)=Q(2)
      ELSE
         QE(-1)=Q(1)
         QE(0)=Q(1)
         QE(N+1)=Q(N)
         QE(N+2)=Q(N)
      ENDIF
C     Edge values, fourth order, limited where they leave the range of
C     the two cells they separate unless the curvature around them is
C     consistent. E(K) is between cells K and K+1.
      DO 20 K=1,N
         IF (K.EQ.N .AND. .NOT.LPER) GOTO 20
         E(K)=7.0/12.0*(QE(K)+QE(K+1))-(QE(K-1)+QE(K+2))/12.0
         IF ((E(K)-QE(K))*(QE(K+1)-E(K)).LT.0.0) THEN
            D2=3.0*(QE(K)-2.0*E(K)+QE(K+1))
            D2L=QE(K-1)-2.0*QE(K)+QE(K+1)
            D2R=QE(K)-2.0*QE(K+1)+QE(K+2)
            SG=SIGN(1.0,D2)
            IF (SG*D2L.GT.0.0 .AND. SG*D2R.GT.0.0) THEN
               D2LM=SG*MIN(CLIM*ABS(D2L),CLIM*ABS(D2R),ABS(D2))
            ELSE
               D2LM=0.0
            ENDIF
            E(K)=0.5*(QE(K)+QE(K+1))-D2LM/3.0
         ENDIF
   20 CONTINUE
      IF (LPER) THEN
         E(0)=E(N)
      ELSE
         E(0)=Q(1)
         E(N)=Q(N)
      ENDIF
C     Cell parabolas q(x)=QL+x*(D+Q6*(1-x)), x in [0,1], D=QR-QL. At an
C     extremum the parabola's curvature is limited by that of its
C     neighbours (a smooth peak keeps it, a spike loses it); elsewhere
C     the edges are pulled in so the parabola stays monotone.
      DO 30 K=1,N
         QL(K)=E(K-1)
         QR(K)=E(K)
         IF ((QR(K)-Q(K))*(Q(K)-QL(K)).LE.0.0 .OR.
     &       (QE(K-1)-Q(K))*(Q(K)-QE(K+1)).LE.0.0) THEN
            D2=-2.0*(6.0*Q(K)-3.0*(QL(K)+QR(K)))
            D2C=QE(K-1)-2.0*Q(K)+QE(K+1)
            D2L=QE(K-2)-2.0*QE(K-1)+Q(K)
            D2R=Q(K)-2.0*QE(K+1)+QE(K+2)
            SG=SIGN(1.0,D2)
            IF (SG*D2C.GT.0.0 .AND. SG*D2L.GT.0.0 .AND.
     &          SG*D2R.GT.0.0) THEN
               D2LM=SG*MIN(ABS(D2),CLIM*ABS(D2C),CLIM*ABS(D2L),
     &                     CLIM*ABS(D2R))
            ELSE
               D2LM=0.0
            ENDIF
            IF (D2.NE.0.0) THEN
               FAC=D2LM/D2
            ELSE
               FAC=0.0
            ENDIF
            QL(K)=Q(K)+(QL(K)-Q(K))*FAC
            QR(K)=Q(K)+(QR(K)-Q(K))*FAC
         ELSE
            IF (ABS(QR(K)-Q(K)).GE.2.0*ABS(QL(K)-Q(K)))
     &         QR(K)=Q(K)-2.0*(QL(K)-Q(K))
            IF (ABS(QL(K)-Q(K)).GE.2.0*ABS(QR(K)-Q(K)))
     &         QL(K)=Q(K)-2.0*(QR(K)-Q(K))
         ENDIF
         D=QR(K)-QL(K)
         Q6(K)=6.0*Q(K)-3.0*(QL(K)+QR(K))
         PMN=MIN(QL(K),QR(K))
         IF (Q6(K).NE.0.0) THEN
            XM=0.5+D/(2.0*Q6(K))
            IF (XM.GT.0.0 .AND. XM.LT.1.0)
     &         PMN=MIN(PMN,QL(K)+XM*(D+Q6(K)*(1.0-XM)))
         ENDIF
         IF (PMN.LT.0.0) THEN
            QL(K)=Q(K)
            QR(K)=Q(K)
            Q6(K)=0.0
         ENDIF
   30 CONTINUE
      DO 40 K=0,N
         IF ((K.EQ.0 .OR. K.EQ.N) .AND. .NOT.LPER) THEN
            G(K)=0.0
         ELSE IF (FS(K).GE.0.0) THEN
            KM=K
            IF (K.EQ.0) KM=N
            C=FS(K)/AM(KM)
            QF=QR(KM)-0.5*C*(QR(KM)-QL(KM)-(1.0-2.0*C/3.0)*Q6(KM))
            G(K)=FS(K)*QF
         ELSE
            KP=K+1
            IF (K.EQ.N) KP=1
            C=-FS(K)/AM(KP)
            QF=QL(KP)+0.5*C*(QR(KP)-QL(KP)+(1.0-2.0*C/3.0)*Q6(KP))
            G(K)=FS(K)*QF
         ENDIF
   40 CONTINUE
      END SUBROUTINE PPMFLUX
C
C**********************************************************
      SUBROUTINE SWEEP1D(N,AM,Q,F,LPER,NSUB)
C     One flux-form 1-D sweep along a line of N cells (the vertical
C     part): air masses AM, mixing ratios Q(N,NDYE), face fluxes F(0:N)
C     as in PPMFLUX. Goes in fractions of F, each sized from the
C     current masses so no cell loses more than half its air mass;
C     NSUB returns how many it took.
      INTEGER N,NSUB
      REAL AM(N),Q(N,NDYE),F(0:N)
      LOGICAL LPER
      REAL G(0:N),FS(0:N),QK(N),CMX,OUTF,AMN,REM,FRAC
      INTEGER K,KK
C
      NSUB=0
      REM=1.0
   5  CONTINUE
      CMX=0.0
      DO 10 K=1,N
         OUTF=MAX(F(K),0.0)+MAX(-F(K-1),0.0)
         CMX=MAX(CMX,OUTF/AM(K))
   10 CONTINUE
      CALL NEXTFRAC(CMX,REM,FRAC,NSUB)
      DO 15 K=0,N
         FS(K)=F(K)*FRAC
   15 CONTINUE
      DO 80 KK=1,NDYE
         DO 20 K=1,N
            QK(K)=Q(K,KK)
   20    CONTINUE
         CALL PPMFLUX(N,AM,QK,FS,LPER,G)
         DO 50 K=1,N
            AMN=AM(K)-(FS(K)-FS(K-1))
            Q(K,KK)=(AM(K)*Q(K,KK)-(G(K)-G(K-1)))/AMN
   50    CONTINUE
   80 CONTINUE
      DO 90 K=1,N
         AM(K)=AM(K)-(FS(K)-FS(K-1))
   90 CONTINUE
      IF (REM.GT.0.0) GOTO 5
      END SUBROUTINE SWEEP1D
C
C**********************************************************
      SUBROUTINE NEXTFRAC(CMX,REM,FRAC,NSUB)
C     Next fraction of a sweep: all that remains (REM) if no cell then
C     loses more than half its air mass (CMX: largest outgoing share
C     for the whole sweep), otherwise as much as keeps to that.
      REAL CMX,REM,FRAC
      INTEGER NSUB
      IF (CMX*REM.LE.0.5) THEN
         FRAC=REM
      ELSE
         FRAC=0.5/CMX
      ENDIF
      REM=REM-FRAC
      IF (REM.LT.1.0E-12) THEN
         FRAC=FRAC+REM
         REM=0.0
      ENDIF
      NSUB=NSUB+1
      IF (NSUB.GT.100000) STOP 'dyeadv: sweep does not converge'
      END SUBROUTINE NEXTFRAC
C
C**********************************************************
      SUBROUTINE HSTEP(L,NSUB)
C     Horizontal transport of every dye in layer L: the Lin and Rood
C     (1996, Mon. Wea. Rev. 124, 2046) combination of the x and y
C     operators (x fluxes from q half-advanced in y and vice versa, the
C     inner half steps in advective form), outer update in flux form.
C     Outgoing dye is scaled down where it would exceed what a cell
C     holds, keeping q>=0. Subcycled like SWEEP1D.
C
C     Polar rows. Near the poles the cells are so narrow that the zonal
C     flux alone can take more than half a cell's air mass in a step,
C     and subcycling the whole layer for those few rows is costly. So
C     first, in any row where a cell's outflow would exceed XLIM of its
C     air mass, the share of the zonal flux beyond that is moved by a
C     plain periodic 1-D sweep (SWEEP1D, subcycled on its own), and the
C     Lin-Rood step gets the rest. Rows that do not need it are
C     untouched, so away from the poles nothing changes.
      INTEGER L,NSUB
C     Layer-sized work arrays on the heap, not the (small) stack of an
C     OpenMP thread.
      REAL, ALLOCATABLE :: FXS(:,:),FYS(:,:),Q0(:,:),QXA(:,:),QYA(:,:)
      REAL, ALLOCATABLE :: GX(:,:),GY(:,:),SCL(:,:),FXH(:,:)
      REAL LN(MG),LNM(MG),LF(0:MG),LG(0:MG)
      REAL CN(2*JG),CNM(2*JG),CF(0:2*JG),CG(0:2*JG)
      REAL CMX,OUTF,REM,FRAC,AMN,AVL
      REAL ALF,XO,YO,APS(MG),QPS(MG,NDYE),FPS(0:MG)
      INTEGER I,R,KK,IM,NSX
C     Outflow share the Lin-Rood step is left with in a pre-swept row,
C     a little under the 1/2 that would make it subcycle.
      REAL, PARAMETER :: XLIM=0.49
C
      ALLOCATE(FXS(MG,2*JG),FYS(MG,0:2*JG),Q0(MG,2*JG),QXA(MG,2*JG),
     &         QYA(MG,2*JG),GX(0:MG,2*JG),GY(MG,0:2*JG),SCL(MG,2*JG),
     &         FXH(MG,2*JG))
C
C     Polar rows: the fraction ALF of the zonal flux that keeps every
C     cell's outflow within XLIM goes to the Lin-Rood step (FXH); the
C     rest is swept first, here.
      DO 4 R=1,2*JG
         ALF=1.0
         DO 2 I=1,MG
            IM=MOD(I+MG-2,MG)+1
            XO=MAX(FX(I,R,L),0.0)+MAX(-FX(IM,R,L),0.0)
            YO=MAX(FY(I,R,L),0.0)+MAX(-FY(I,R-1,L),0.0)
            IF (XO+YO.GT.XLIM*AM(I,R,L) .AND. XO.GT.0.0)
     &         ALF=MIN(ALF,MAX(XLIM*AM(I,R,L)-YO,0.0)/XO)
    2    CONTINUE
         IF (ALF.LT.1.0) THEN
            DO 3 I=1,MG
               FPS(I)=(1.0-ALF)*FX(I,R,L)
               FXH(I,R)=ALF*FX(I,R,L)
               APS(I)=AM(I,R,L)
               DO 3 KK=1,NDYE
                  QPS(I,KK)=QDYE(I,R,L,KK+1)
    3       CONTINUE
            FPS(0)=FPS(MG)
            CALL SWEEP1D(MG,APS,QPS,FPS,.TRUE.,NSX)
            DO 1 I=1,MG
               AM(I,R,L)=APS(I)
               DO 1 KK=1,NDYE
                  QDYE(I,R,L,KK+1)=QPS(I,KK)
    1       CONTINUE
         ELSE
            DO 6 I=1,MG
               FXH(I,R)=FX(I,R,L)
    6       CONTINUE
         ENDIF
    4 CONTINUE
C
      NSUB=0
      REM=1.0
   5  CONTINUE
      CMX=0.0
      DO 10 R=1,2*JG
         DO 10 I=1,MG
            IM=MOD(I+MG-2,MG)+1
            OUTF=MAX(FXH(I,R),0.0)+MAX(-FXH(IM,R),0.0)
     &          +MAX(FY(I,R,L),0.0)+MAX(-FY(I,R-1,L),0.0)
            CMX=MAX(CMX,OUTF/AM(I,R,L))
   10 CONTINUE
      CALL NEXTFRAC(CMX,REM,FRAC,NSUB)
      DO 15 R=1,2*JG
         DO 15 I=1,MG
            FXS(I,R)=FXH(I,R)*FRAC
   15 CONTINUE
      DO 16 R=0,2*JG
         DO 16 I=1,MG
            FYS(I,R)=FY(I,R,L)*FRAC
   16 CONTINUE
C
      DO 100 KK=2,NDYE+1
         DO 18 R=1,2*JG
            DO 18 I=1,MG
               Q0(I,R)=QDYE(I,R,L,KK)
   18    CONTINUE
C        Inner advective half steps: QXA = q + x-update/2, QYA likewise
C        in y, each without the 1-D divergence term.
         DO 30 R=1,2*JG
            DO 22 I=1,MG
               LN(I)=Q0(I,R)
               LNM(I)=AM(I,R,L)
               LF(I)=FXS(I,R)
   22       CONTINUE
            LF(0)=LF(MG)
            CALL PPMFLUX(MG,LNM,LN,LF,.TRUE.,LG)
            DO 24 I=1,MG
               QXA(I,R)=Q0(I,R)+0.5*(-(LG(I)-LG(I-1))
     &                  +Q0(I,R)*(LF(I)-LF(I-1)))/LNM(I)
   24       CONTINUE
   30    CONTINUE
         DO 40 I=1,MG
            DO 32 R=1,2*JG
               CN(R)=Q0(I,R)
               CNM(R)=AM(I,R,L)
               CF(R)=FYS(I,R)
   32       CONTINUE
            CF(0)=FYS(I,0)
            CALL PPMFLUX(2*JG,CNM,CN,CF,.FALSE.,CG)
            DO 34 R=1,2*JG
               QYA(I,R)=Q0(I,R)+0.5*(-(CG(R)-CG(R-1))
     &                  +Q0(I,R)*(CF(R)-CF(R-1)))/CNM(R)
   34       CONTINUE
   40    CONTINUE
C        Outer fluxes: x from q half-advanced in y, y from q half-
C        advanced in x (floored at 0; they only set face values).
         DO 50 R=1,2*JG
            DO 42 I=1,MG
               LN(I)=MAX(QYA(I,R),0.0)
               LNM(I)=AM(I,R,L)
               LF(I)=FXS(I,R)
   42       CONTINUE
            LF(0)=LF(MG)
            CALL PPMFLUX(MG,LNM,LN,LF,.TRUE.,LG)
            DO 44 I=0,MG
               GX(I,R)=LG(I)
   44       CONTINUE
   50    CONTINUE
         DO 60 I=1,MG
            DO 52 R=1,2*JG
               CN(R)=MAX(QXA(I,R),0.0)
               CNM(R)=AM(I,R,L)
               CF(R)=FYS(I,R)
   52       CONTINUE
            CF(0)=FYS(I,0)
            CALL PPMFLUX(2*JG,CNM,CN,CF,.FALSE.,CG)
            DO 54 R=0,2*JG
               GY(I,R)=CG(R)
   54       CONTINUE
   60    CONTINUE
C        Positivity: scale each cell's outgoing dye to what it holds.
         DO 70 R=1,2*JG
            DO 70 I=1,MG
               OUTF=MAX(GX(I,R),0.0)+MAX(-GX(I-1,R),0.0)
     &             +MAX(GY(I,R),0.0)+MAX(-GY(I,R-1),0.0)
C              Rounding can leave q a hair below zero; such a cell has
C              nothing to give (and with no outflow, AVL/OUTF would be
C              -Infinity, and NaN once it meets a zero flux).
               AVL=MAX(AM(I,R,L)*Q0(I,R),0.0)
               IF (OUTF.GT.AVL) THEN
                  SCL(I,R)=AVL/OUTF
               ELSE
                  SCL(I,R)=1.0
               ENDIF
   70    CONTINUE
         DO 75 R=1,2*JG
            DO 72 I=1,MG
               IF (GX(I,R).GT.0.0) THEN
                  GX(I,R)=GX(I,R)*SCL(I,R)
               ELSE
                  GX(I,R)=GX(I,R)*SCL(MOD(I,MG)+1,R)
               ENDIF
   72       CONTINUE
            GX(0,R)=GX(MG,R)
   75    CONTINUE
         DO 78 R=1,2*JG-1
            DO 78 I=1,MG
               IF (GY(I,R).GT.0.0) THEN
                  GY(I,R)=GY(I,R)*SCL(I,R)
               ELSE
                  GY(I,R)=GY(I,R)*SCL(I,R+1)
               ENDIF
   78    CONTINUE
         DO 80 R=1,2*JG
            DO 80 I=1,MG
               IM=MOD(I+MG-2,MG)+1
               AMN=AM(I,R,L)-(FXS(I,R)-FXS(IM,R))-(FYS(I,R)-FYS(I,R-1))
               QDYE(I,R,L,KK)=(AM(I,R,L)*Q0(I,R)-(GX(I,R)-GX(I-1,R))
     &                        -(GY(I,R)-GY(I,R-1)))/AMN
   80    CONTINUE
  100 CONTINUE
      DO 110 R=1,2*JG
         DO 110 I=1,MG
            IM=MOD(I+MG-2,MG)+1
            AM(I,R,L)=AM(I,R,L)-(FXS(I,R)-FXS(IM,R))
     &                -(FYS(I,R)-FYS(I,R-1))
  110 CONTINUE
      IF (REM.GT.0.0) GOTO 5
      END SUBROUTINE HSTEP
C
      END MODULE dyegrid
