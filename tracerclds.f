C**********************************************************
C             MODULE TRACERCLDS
C**********************************************************
      MODULE tracerclds
C
C     State for the settling dye tracers and their optional coupling
C     into the radiative transfer.
C
      include 'params.i'
      PARAMETER(MGPP=MG+2,IGC=MGPP*NHEM)
C     Most vapour/solid pairs (KDYEMODE=2), hence most radiatively
C     active dyes, NDYE allows: one cloud species per pair.
      PARAMETER(NPAIRMX=MAX(NDYE/2,1))
      SAVE
C
C     Grain radius per dye, indexed 2..NDYE+1; ADYE(1) is unused.
      REAL ADYE(NDYE+1)
C     Grain size distribution: ADYE is the median radius of a number-
C     weighted lognormal of width SIGDYE, as in the cloud optics, whose
C     Mie tables (CLOUD_KEXT etc.) and exp_92_lnsig2_pi correction were
C     built for it; SIGDYE must match the sigma those tables used.
C     LDYELOGN (INVARPARAM, default T): settle that distribution, with
C     its mass-weighted mean fall speed; F settles single grains of
C     radius ADYE, as before (and bit-identical to it).
      REAL, PARAMETER :: SIGDYE=1.65
      LOGICAL LDYELOGN
C     8-point Gauss-Hermite rule (weight exp(-x**2)) for that mean.
      INTEGER, PARAMETER :: NQDYE=8
      REAL, PARAMETER :: XQDYE(NQDYE)=(/-2.9306374202572441,
     &   -1.981656756695843,-1.1571937124467802,-0.38118699020732211,
     &   0.38118699020732211,1.1571937124467802,1.981656756695843,
     &   2.9306374202572441/)
      REAL, PARAMETER :: WQDYE(NQDYE)=(/0.00019960407221136783,
     &   0.017077983007413467,0.20780232581489183,0.66114701255824149,
     &   0.66114701255824149,0.20780232581489183,0.017077983007413467,
     &   0.00019960407221136783/)
C     Grain density of each pair's species (index as KDYESPEC), set to
C     DENSITY(KDYESPEC(P)) when the cloud tables load, so the settling
C     and the OPPRMULTI optics cannot describe different particles.
      REAL RHODYE(NPAIRMX)
C     Below PDYEFIX the dyes relax toward 1 on TRELAXORB (orbits). Above
C     PDYEUPPER they start at 1 (initial state only).
      REAL PDYEFIX,PDYEUPPER,TRELAXORB
C
C     --- radiatively active dye ------------------------------------
C     KDYESPEC(P): the cloud species (indexed as MOLEF is: 7 = Mg2SiO4,
C     13 = Al2O3) of vapour/solid pair P in KDYEMODE=2, one species per
C     pair; KDYEMODE=1 uses KDYESPEC(1) for every dye. LDYERAD off
C     leaves the RT on the prescribed condensation-curve clouds. On,
C     the dyes place the clouds instead: in KDYEMODE=2 every pair's
C     solid, in KDYEMODE=1 the single tracer KDYERAD; each with its own
C     ADYE as its optical grain size.
      INTEGER KDYERAD,KDYESPEC(NPAIRMX)
      LOGICAL LDYERAD
C     The radiatively active dyes, set in INIVARPARAM: NDYERAD of them,
C     number N in dye slot KSLOTRAD(N) and of species KDYESPEC(N);
C     ISPRAD(I) is the one of species I (0 = none).
      INTEGER NDYERAD,KSLOTRAD(NPAIRMX),ISPRAD(NCLOUDS)
C
C     --- cloud scheme ----------------------------------------------
C     KDYEMODE=1: each dye is total species material, and is all-vapor 
C     or all-solid depending on T from FCONDYE. KDYEMODE=2: dyes come in
C     vapour/solid pairs, vapour in slot 2P and solid in 2P+1, so
C     NDYE=2*(number of pairs). The pair exchanges mass by relaxing
C     to the QSATDYE equilibrium on TCONDDYE (supersaturated) or
C     TEVAPDYE (subsaturated), both in seconds. Only the solid falls
C     and only the solid is opaque; its radius is ADYE(2P+1).
      INTEGER KDYEMODE
      REAL TCONDDYE,TEVAPDYE
C
C     T: rescale each dye once a step (dyeadv.f) so its mass changes
C     only through the physics, not through the finite-volume vs
C     spectral p* mismatch; with KDYEMODE=2 each vapour/solid pair is
C     rescaled together, by one factor, so it is the pair's total that
C     is held.
      LOGICAL LDYEFIX
C
C     Budget of the dye mass above PDYEBUD (Pa; <=0 turns it off), per
C     process, written to dye_budget.txt every NDYEBUD steps (0: with
C     the lines of dye_adv.txt). Diagnostic only.
      REAL PDYEBUD
      INTEGER NDYEBUD
C
C     Benchmark mode (INVARPARAM; not dye-specific, kept here as the
C     module the namelist already fills): T skips the snapshot series
C     written in the main loop, fort.26xx/50xx/64xx over the last orbit
C     and fort.27xx/70xx/90xx every NSNAP steps, so timings are not
C     dominated by output. Default F.
      LOGICAL LBENCH
C
C     Single-column RT mode (INVARPARAM; not dye-specific either): T runs
C     the RT once, for one column with the initial-condition profile
C     (RESTTT of fort.7, via TRS in SETZT) at p*=P0 and stellar zenith
C     angle ZEN1COL (degrees), writes rt1col.txt and rt1col_bands.txt,
C     and stops. Default F. T1COL: that profile in K, top to bottom.
      LOGICAL L1COLRT
      REAL ZEN1COL
      REAL T1COL(NL)
C
C     Cloud ramp (INCLOUDY; not dye-specific): every cloud and haze
C     optical depth, prescribed or dye, is scaled by model day/CLDRAMP
C     until day CLDRAMP. Counted from the start of the run (KOUNT is
C     restored on restart). 0 (default) = no ramp.
      REAL CLDRAMP
C
C     --- grid-point dye state ----------------------------------------
C     Mixing ratios q on the Gaussian grid, rows 1..2*JG north to south:
C     latitude JH of the model's hemisphere loop is row JH in the north
C     (IOFM=0) and 2*JG+1-JH in the south (IOFM=MGPP).
      REAL QDYE(MG,2*JG,NL,2:NDYE+1)
C     Winds and surface pressure stored by MGRMLT at the current time
C     level t, as the model carries them: UDY=u cos(lat)/(a*Omega),
C     VDY=v cos(lat)/(a*Omega), PSDYE=p*/P0. UDYM/VDYM hold the winds of
C     the previous step, for extrapolating to the middle of the step.
      REAL UDY(MG,2*JG,NL),VDY(MG,2*JG,NL)
      REAL UDYM(MG,2*JG,NL),VDYM(MG,2*JG,NL)
      REAL PSDYE(MG,2*JG)
C     Per row: sin(latitude) and the Gaussian weight, which is the
C     row's width in sin(latitude) (the weights over all rows sum to 2).
      REAL AMUDY(2*JG),WDY(2*JG)
C     LDYEINI: QDYE holds a state. LDYEPRV: UDYM/VDYM hold the winds of
C     a previous step of length DTDYM. DAYDY: model day of QDYE. DTDY:
C     length of the latest dye step in model time (1/Omega), which the
C     physics in DGRMLT also integrates over.
      LOGICAL :: LDYEINI=.FALSE.,LDYEPRV=.FALSE.
      REAL DAYDY,DTDY,DTDYM
C     Fixer bookkeeping, per dye, in units of sum(p*/P0*DSIGMA*WDY*q):
C     DYMTGT is the mass the last advection left, DMPHYS what the
C     physics has added since. LDYMTGT: DYMTGT is set.
      REAL DYMTGT(2:NDYE+1),DMPHYS(2:NDYE+1)
      LOGICAL :: LDYMTGT=.FALSE.
C     Wall-clock seconds since the last line of dye_adv.txt spent in
C     the dye physics (DGRMLT), for the timing columns written there.
      REAL :: TDYPHY=0.0
C
C     --- work arrays ----------------------------------------------
C     The radiatively active dyes in the model's (IGC,NL,JG) grid
C     layout, copied from QDYE for RADIATION_ALLLATS.
      REAL TRAG_forrad(IGC,NL,JG,NPAIRMX)
C     One column's opaque amount of each active dye per RT layer, handed
C     to OPPRMULTI (filled, aligned to the RT layers, in
C     RADIATION_ALLLATS). Threadprivate because the RT runs its columns
C     in parallel.
      REAL QDYECOL(NL+1,NPAIRMX)
!$OMP THREADPRIVATE(QDYECOL)
C
C     --- condensation curves of the dye species --------------------
C     TCDYE(:,P) is TCONDS(MET_INDEX,:,KDYESPEC(P)), on the pressure
C     grid PCDYE in Pa, copied out in cloud_properties_set_up.f once
C     MET_INDEX is final. A KDYEMODE=1 dye is TOTAL species material,
C     gas plus solid, so both its settling in DGRMLT and its opacity in
C     OPPRMULTI take the condensed fraction from FCONDYE; a KDYEMODE=2
C     pair takes its equilibrium from QSATDYE. Only solids fall, and
C     only solids are opaque.
      REAL TCDYE(80,NPAIRMX),PCDYE(80)
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
      REAL FUNCTION FCONDYE(PPA,TK,IP)
C     Condensed fraction of species KDYESPEC(IP) at pressure PPA (Pa)
C     and temperature TK (K): 0 all vapour, 1 all solid, over the same
C     10 K ramp OPPRMULTI uses for its own species.
      REAL PPA,TK
      INTEGER IP
      FCONDYE=MIN(MAX((TCDYE(IPCDYE(PPA),IP)-TK)/10.0,0.0),1.0)
      END FUNCTION FCONDYE
C
      REAL FUNCTION QSATDYE(PPA,TK,IP)
C     Saturation abundance of species KDYESPEC(IP), in the dyes' units
C     (1 = deep abundance), at pressure PPA (Pa) and temperature TK (K).
C     The condensation curve is where an abundance-1 gas saturates, so
C     the saturation vapour pressure at TK is that gas's partial
C     pressure at P*, the pressure where TCDYE(P*)=TK; hence
C     QSATDYE=P*/PPA. Between table points ln P is interpolated
C     linearly in 1/T (Clausius-Clapeyron); beyond the table the
C     end-to-end slope is used. Needs TCDYE monotonic in P, which
C     cloud_properties_set_up.f checks when KDYEMODE=2.
      REAL PPA,TK,RT,SLOPE,PLNST
      INTEGER IP,ILO,IHI,IMID,IEND
      LOGICAL LUP
      RT=1.0/TK
      LUP=TCDYE(80,IP).GT.TCDYE(1,IP)
      IF ((TK.LE.TCDYE(1,IP)).EQV.LUP .OR. TK.EQ.TCDYE(1,IP)) THEN
         IEND=1
      ELSE IF ((TK.GE.TCDYE(80,IP)).EQV.LUP .OR. TK.EQ.TCDYE(80,IP)) THEN
         IEND=80
      ELSE
         IEND=0
      ENDIF
      IF (IEND.NE.0) THEN
         SLOPE=(LOG(PCDYE(80))-LOG(PCDYE(1)))
     &        /(1.0/TCDYE(80,IP)-1.0/TCDYE(1,IP))
         PLNST=LOG(PCDYE(IEND))+(RT-1.0/TCDYE(IEND,IP))*SLOPE
      ELSE
C        Bisect for TCDYE(ILO,IP) and TCDYE(IHI,IP) bracketing TK.
         ILO=1
         IHI=80
   10    IF (IHI-ILO.GT.1) THEN
            IMID=(ILO+IHI)/2
            IF ((TCDYE(IMID,IP).LT.TK).EQV.LUP) THEN
               ILO=IMID
            ELSE
               IHI=IMID
            ENDIF
            GOTO 10
         ENDIF
         SLOPE=(LOG(PCDYE(IHI))-LOG(PCDYE(ILO)))
     &        /(1.0/TCDYE(IHI,IP)-1.0/TCDYE(ILO,IP))
         PLNST=LOG(PCDYE(ILO))+(RT-1.0/TCDYE(ILO,IP))*SLOPE
      ENDIF
      QSATDYE=EXP(PLNST)/PPA
      END FUNCTION QSATDYE
C
      SUBROUTINE SETTLECOL(QIN,SIG,SIGH,DSIG,THCOL,AD,RHOD,GA,GASCON,
     &                     WW,DELT2,FSCL,GLAM,GRHO,GETA,QOUT)
C     One column's gravitational settling over DELT2: backward-Euler
C     upwind sweep in sigma-mass form,
C        DSIG(L)*q(L) = DSIG(L)*qin(L) + F(L-1)*qh(L-1) - F(L)*qh(L),
C     F(0)=0. F(L) is the flux coefficient at layer L's lower interface,
C     rho*VFALL*DELT2 over the column mass per unit sigma (ps/g),
C        F(L) = GA*SIGH(L)*VFALL(L)/(GASCON*THCOL(L)) * DELT2/WW,
C     every factor at the interface: SIGH(L) its sigma (1 at the
C     surface), THCOL(L) its temperature. qh(L) is q reconstructed at
C     that interface from upstream, linear in ln(sigma):
C        qh(L) = q(L) + FH*PHI*(q(L)-q(L-1)),
C     FH the interface's distance from SIG(L) in units of the spacing
C     from SIG(L-1) (1/2 on the log grid). PHI is a minmod limiter on
C     the slopes either side of L, taken from QIN so that qh stays
C     linear in the unknowns and grains, which only fall, still allow
C     one top-down pass. A cell whose qh would come out negative is
C     re-solved with PHI=0, so q stays non-negative.
C     First-order donor cell (PHI=0) on this grid, 1.38 per layer in
C     sigma, errs by sqrt(1.38) either way: upwinding q, a uniform q
C     settles at the right rate but a cloud layer descends 17% too
C     fast; upwinding q*SIG instead (the original scheme), the layer
C     moves at the right speed but a uniform q settles 15% too slowly.
C     The reconstruction gets both within a few per cent. Each F*qh
C     appears once with each sign, so column mass changes only by
C     what leaves the bottom.
C     AD is the grain radius in m and RHOD its density in kg/m3 (the
C     species' RHODYE). With LDYELOGN, AD is the median radius of the
C     lognormal number distribution of width SIGDYE, and VFALL is its
C     mass-weighted mean: weighted by mass, ln r is again normal, with
C     the same width s=ln(SIGDYE) and median AD*exp(3*s**2), so
C        VFALL = sum_k WQDYE(k)*v(AD*exp(3*s**2+sqrt(2)*s*XQDYE(k)))
C                / sqrt(pi),
C     v the single-grain speed of V_FALL_GRAIN (to <3e-4 of the exact
C     mean, except near its 100 m/s cap). The caller supplies the grain-
C     independent quantities, so a column computes them once for all
C     its dyes: GLAM, GRHO, GETA are the gas mean free path, density
C     and viscosity from SED_GAS at each interface. FSCL belongs to the
C     donor layer above the interface and scales VFALL: FCONDYE for a
C     total-material dye (KDYEMODE=1), of which only the condensed part
C     falls, and 1 for a solid tracer.
      REAL QIN(NL),SIG(NL),SIGH(NL),DSIG(NL),THCOL(NL),QOUT(NL)
      REAL FSCL(NL),GLAM(NL),GRHO(NL),GETA(NL)
      REAL AD,RHOD,GA,GASCON,WW,DELT2
      REAL VFALL,FSED,FLXUP,FH,PHI,SUP,SDN,C1,QH
      REAL RQ(NQDYE),WQ(NQDYE),SL,VK
      INTEGER L,K
C     The quadrature radii and weights, the same at every level.
      IF (LDYELOGN) THEN
         SL=LOG(SIGDYE)
         DO 10 K=1,NQDYE
            RQ(K)=AD*EXP(3.0*SL*SL+SQRT(2.0)*SL*XQDYE(K))
            WQ(K)=WQDYE(K)/SQRT(4.0*ATAN(1.0))
   10    CONTINUE
      ENDIF
      FLXUP=0.0
      DO 300 L=1,NL
         IF (LDYELOGN) THEN
            VFALL=0.0
            DO 20 K=1,NQDYE
               CALL V_FALL_GRAIN(RQ(K),GA,RHOD,GLAM(L),GRHO(L),GETA(L),
     &                           VK)
               VFALL=VFALL+WQ(K)*VK
   20       CONTINUE
         ELSE
            CALL V_FALL_GRAIN(AD,GA,RHOD,GLAM(L),GRHO(L),GETA(L),VFALL)
         ENDIF
         FSED=GA*SIGH(L)*VFALL*FSCL(L)/(GASCON*THCOL(L))*DELT2/WW
C        Limited upstream slope; none in the top layer, and the bottom
C        layer, with nothing below, takes its upstream slope unlimited.
         IF (L.EQ.1) THEN
            QOUT(L)=DSIG(L)*QIN(L)/(DSIG(L)+FSED)
            QH=QOUT(L)
         ELSE
            FH=LOG(SIGH(L)/SIG(L))/LOG(SIG(L)/SIG(L-1))
            SUP=(QIN(L)-QIN(L-1))/LOG(SIG(L)/SIG(L-1))
            IF (L.LT.NL) THEN
               SDN=(QIN(L+1)-QIN(L))/LOG(SIG(L+1)/SIG(L))
            ELSE
               SDN=SUP
            ENDIF
            IF (SUP*SDN.GT.0.0) THEN
               PHI=MIN(1.0,SDN/SUP)
            ELSE
               PHI=0.0
            ENDIF
            C1=FH*PHI
            QOUT(L)=(DSIG(L)*QIN(L)+FLXUP+FSED*C1*QOUT(L-1))
     &             /(DSIG(L)+FSED*(1.0+C1))
            QH=QOUT(L)+C1*(QOUT(L)-QOUT(L-1))
            IF (QH.LT.0.0) THEN
               C1=0.0
               QOUT(L)=(DSIG(L)*QIN(L)+FLXUP)/(DSIG(L)+FSED)
               QH=QOUT(L)
            ENDIF
         ENDIF
         FLXUP=FSED*QH
  300 CONTINUE
      END SUBROUTINE SETTLECOL
C
      END MODULE tracerclds
