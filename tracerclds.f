C**********************************************************
C             MODULE TRACERCLDS
C**********************************************************
      MODULE tracerclds
C
C     State for the settling dye tracers and their optional coupling
C     into the radiative transfer. This replaces the DYEPAR, DYERAD,
C     DYEDEN, DYEPRV and DYECOL common blocks, which between them had
C     to be repeated in seventeen places: a module is checked by the
C     compiler, whereas a mismatched COMMON is silent corruption.
C
C     No IMPLICIT NONE: params.i is a bare PARAMETER statement relying
C     on implicit integer typing, as everywhere else in this model.
C
      include 'params.i'
      PARAMETER(MGPP=MG+2,IGC=MGPP*NHEM)
      SAVE
C
C     --- namelist configuration (INVARPARAM in fort.7) -------------
C     ADYE is indexed by TRACER number, so ADYE(1) is the water vapour
C     slot and is never used; set ADYE(2:NTRAC), one radius per dye.
      REAL ADYE(NTRAC)
C     Grain density. Read from the namelist, then overwritten with
C     DENSITY(KDYESPEC) when the cloud tables load, so the settling and
C     the OPPRMULTI optics cannot describe different particles.
      REAL RHODYE
C     Deep reservoir below PDYEFIX, upper-atmosphere fill above
C     PDYEUPPER, relaxation timescale in orbits.
      REAL PDYEFIX,PDYEUPPER,TRELAXORB
C
C     --- radiatively active dye ------------------------------------
C     LDYERAD off leaves the RT on the prescribed condensation-curve
C     clouds. On, the single tracer KDYERAD places the cloud and
C     ADYE(KDYERAD) sets its optical grain size, with KDYESPEC
C     selecting which cloud species supplies the optics (indexed as
C     MOLEF is: 7 = Mg2SiO4).
      INTEGER KDYERAD,KDYESPEC
      LOGICAL LDYERAD
C
C     --- cloud scheme ----------------------------------------------
C     KDYEMODE=1: each dye is total species material, and is all-vapor 
C     or all-solid depending on T from FCONDYE. KDYEMODE=2: dyes come in
C     vapour/solid pairs, vapour in slot 2P and solid in 2P+1, so
C     NTRAC=1+2*(number of pairs). The pair exchanges mass by relaxing
C     to the QSATDYE equilibrium on TCONDDYE (supersaturated) or
C     TEVAPDYE (subsaturated), both in seconds. Only the solid falls
C     and only the solid is opaque; its radius is ADYE(2P+1).
      INTEGER KDYEMODE
      REAL TCONDDYE,TEVAPDYE
C
C     --- work arrays ----------------------------------------------
C     Snapshot of the radiatively active dye for RADIATION_ALLLATS,
C     which loops every latitude while TRAG holds only the current one.
      REAL TRAG_forrad(IGC,NL,JG)
C     One column's dye profile, handed to OPPRMULTI. Threadprivate
C     because the RT runs its columns in parallel.
      REAL QDYECOL(NL+1)
!$OMP THREADPRIVATE(QDYECOL)
C
C     --- condensation curve for the dye species --------------------
C     TCONDS(MET_INDEX,:,KDYESPEC) and its pressure grid in Pa, copied
C     out in cloud_properties_set_up.f once MET_INDEX is final. The dye
C     is TOTAL species material, gas plus solid, so both the settling
C     in DGRMLT and the opacity in OPPRMULTI have to take the condensed
C     fraction from FCONDYE below; only solids fall, and only solids
C     are opaque. Held here, and read through one function, so the two
C     cannot drift apart.
      REAL TCDYE(80),PCDYE(80)
C     PCDYE is uniform in ln P (to the 4 digits it is tabulated with),
C     so IPCDYE finds a pressure's bracket from its log rather than by
C     searching: PCLN1=LOG(PCDYE(1)), RDLNPC=79/(LOG(PCDYE(80))-PCLN1),
C     both set beside PCDYE in cloud_properties_set_up.f.
      REAL PCLN1,RDLNPC
C
      CONTAINS
C
      INTEGER FUNCTION IPCDYE(PPA)
C     Index of the PCDYE entry nearest PPA (Pa).
      REAL PPA
      INTEGER IP
      IF (PPA.LE.PCDYE(1)) THEN
         IPCDYE=1
         RETURN
      ENDIF
      IF (PPA.GE.PCDYE(80)) THEN
         IPCDYE=80
         RETURN
      ENDIF
      IP=MIN(MAX(1+INT((LOG(PPA)-PCLN1)*RDLNPC),1),79)
   10 IF (IP.GT.1 .AND. PPA.LT.PCDYE(IP)) THEN
         IP=IP-1
         GOTO 10
      ENDIF
   20 IF (IP.LT.79 .AND. PPA.GT.PCDYE(IP+1)) THEN
         IP=IP+1
         GOTO 20
      ENDIF
      IF (ABS(PCDYE(IP+1)-PPA).LT.ABS(PCDYE(IP)-PPA)) IP=IP+1
      IPCDYE=IP
      END FUNCTION IPCDYE
C
      REAL FUNCTION FCONDYE(PPA,TK)
C     Condensed fraction of the dye species at pressure PPA (Pa) and
C     temperature TK (K): 0 all vapour, 1 all solid, over the same 10 K
C     ramp OPPRMULTI uses for its own species.
      REAL PPA,TK
      FCONDYE=MIN(MAX((TCDYE(IPCDYE(PPA))-TK)/10.0,0.0),1.0)
      END FUNCTION FCONDYE
C
      REAL FUNCTION QSATDYE(PPA,TK)
C     Saturation abundance of the dye species, in the dyes' units
C     (1 = deep abundance), at pressure PPA (Pa) and temperature TK (K).
C     The condensation curve is where an abundance-1 gas saturates, so
C     the saturation vapour pressure at TK is that gas's partial
C     pressure at P*, the pressure where TCDYE(P*)=TK; hence
C     QSATDYE=P*/PPA. Between table points ln P is interpolated
C     linearly in 1/T (Clausius-Clapeyron); beyond the table the
C     end-to-end slope is used. Needs TCDYE monotonic in P, which
C     cloud_properties_set_up.f checks when KDYEMODE=2.
      REAL PPA,TK,RT,SLOPE,PLNST
      INTEGER ILO,IHI,IMID,IEND
      LOGICAL LUP
      RT=1.0/TK
      LUP=TCDYE(80).GT.TCDYE(1)
      IF ((TK.LE.TCDYE(1)).EQV.LUP .OR. TK.EQ.TCDYE(1)) THEN
         IEND=1
      ELSE IF ((TK.GE.TCDYE(80)).EQV.LUP .OR. TK.EQ.TCDYE(80)) THEN
         IEND=80
      ELSE
         IEND=0
      ENDIF
      IF (IEND.NE.0) THEN
         SLOPE=(LOG(PCDYE(80))-LOG(PCDYE(1)))
     &        /(1.0/TCDYE(80)-1.0/TCDYE(1))
         PLNST=LOG(PCDYE(IEND))+(RT-1.0/TCDYE(IEND))*SLOPE
      ELSE
C        Bisect for TCDYE(ILO) and TCDYE(IHI) bracketing TK.
         ILO=1
         IHI=80
   10    IF (IHI-ILO.GT.1) THEN
            IMID=(ILO+IHI)/2
            IF ((TCDYE(IMID).LT.TK).EQV.LUP) THEN
               ILO=IMID
            ELSE
               IHI=IMID
            ENDIF
            GOTO 10
         ENDIF
         SLOPE=(LOG(PCDYE(IHI))-LOG(PCDYE(ILO)))
     &        /(1.0/TCDYE(IHI)-1.0/TCDYE(ILO))
         PLNST=LOG(PCDYE(ILO))+(RT-1.0/TCDYE(ILO))*SLOPE
      ENDIF
      QSATDYE=EXP(PLNST)/PPA
      END FUNCTION QSATDYE
C
      SUBROUTINE SETTLECOL(QIN,PCOL,TCOL,SIG,DSIG,AD,GA,GASCON,WW,
     &                     DELT2,FSCL,GLAM,GRHO,GETA,QOUT)
C     One column's gravitational settling over DELT2: backward-
C     Euler donor-cell sweep,
C        q(L) = [ qin(L) + GSED*q(L-1) ] / (1 + CSED),  q(0)=0,
C     bidiagonal because grains only fall, so one top-down pass.
C     PCOL in Pa and TCOL in K are the layer pressures and
C     temperatures, AD the grain radius in m. The caller supplies the
C     per-level quantities that do not depend on the grain, so a column
C     computes them once for all its dyes: GLAM, GRHO, GETA are the gas
C     mean free path, density and viscosity from SED_GAS, and FSCL
C     scales VFALL - FCONDYE for a total-material dye (KDYEMODE=1), of
C     which only the condensed part falls, and 1 for a solid tracer.
      REAL QIN(NL),PCOL(NL),TCOL(NL),SIG(NL),DSIG(NL),QOUT(NL)
      REAL FSCL(NL),GLAM(NL),GRHO(NL),GETA(NL)
      REAL AD,GA,GASCON,WW,DELT2
      REAL UPSED,PSED,TSED,VFALL,DZSED,CSED,GSED,PSEDU,TSEDU,VFUP
      INTEGER L
      UPSED=0.0
      DO 300 L=1,NL
         PSED=PCOL(L)
         TSED=TCOL(L)
         CALL V_FALL_GRAIN(AD,GA,RHODYE,GLAM(L),GRHO(L),GETA(L),VFALL)
         DZSED=(GASCON*TSED/GA)*(DSIG(L)/SIG(L))
         VFALL=VFALL*FSCL(L)
         CSED=VFALL/DZSED/WW*DELT2
         IF (L.EQ.1) THEN
            GSED=0.0
         ELSE
            GSED=(TSED*PSEDU)/(PSED*TSEDU)
     +           *VFUP/DZSED/WW*DELT2
         ENDIF
         QOUT(L)=(QIN(L)+GSED*UPSED)/(1.0+CSED)
         UPSED=QOUT(L)
         PSEDU=PSED
         TSEDU=TSED
         VFUP=VFALL
  300 CONTINUE
      END SUBROUTINE SETTLECOL
C
      END MODULE tracerclds
