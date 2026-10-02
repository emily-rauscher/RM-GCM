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
C     --- work arrays ----------------------------------------------
C     Previous timestep's dye field, for the implicit settling sweep.
      REAL TRAPRV(IGC,NL,JG,NTRAC-1)
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
C
      CONTAINS
C
      REAL FUNCTION FCONDYE(PPA,TK)
C     Condensed fraction of the dye species at pressure PPA (Pa) and
C     temperature TK (K): 0 all vapour, 1 all solid, over the same 10 K
C     ramp OPPRMULTI uses for its own species.
      REAL PPA,TK
      INTEGER IP
      IP=MINLOC(ABS(PCDYE-PPA),1)
      FCONDYE=MIN(MAX((TCDYE(IP)-TK)/10.0,0.0),1.0)
      END FUNCTION FCONDYE
C
      END MODULE tracerclds
