# Dependencies are coarse but complete (see the rules after OBJS): any
# header or module change recompiles everything that could use it, so
# "make clean" is never needed for correctness and "make -j" is safe.
# Model configuration.  The array sizes (params.i) and the build record
# (buildinfo.inc) are generated from these settings by mkparams.sh, and
# the settings are kept in build.cfg, so a plain "make" rebuilds the
# same configuration and a changed one recompiles everything.
#   make RES=T2 NL=100 RT=picket   T2, 100 levels, picket fence/double-gray
#   make RES=T31 NL=50 RT=11       T31, 50 levels, 11-band corr-k
#   RES     horizontal resolution: T2 or T31
#   NL      number of vertical levels
#   RT      picket, or a corr-k band count (11, 30, 196, ...)
#   NDYE    number of cloud dye tracers
#   NTGRID, NPGRID  size of the corr-k k-table T-P grid
# A fresh checkout (no build.cfg) builds T31 L50 picket fence.
-include build.cfg
RES ?= T31
NL ?= 50
RT ?= 11
NDYE ?= 2
NTGRID ?= 73
NPGRID ?= 20

# Compiler selection.  ifort is the default because it gives the faster
# executable.  It needs the oneAPI environment loaded first
# ("source /opt/intel/oneapi/setvars.sh").  It is an x86_64 binary, so
# on Apple Silicon it runs under Rosetta and stops working once Rosetta
# is gone ("Bad CPU type in executable"); gfortran is the fallback there.
#   make                     -> ifort (default)
#   make COMPILER=gfortran   -> gfortran
COMPILER ?= ifort

ifeq ($(COMPILER),gfortran)

FC = gfortran
# Homebrew GCC derives -mmacosx-version-min from the running OS, and on
# macOS 26+ it emits a version clang rejects ("invalid version number in
# '-mmacosx-version-min=18.0'").  Pin a target clang accepts; override on
# the command line if you need a different one.
MACOSX_DEPLOYMENT_TARGET ?= 15.0
export MACOSX_DEPLOYMENT_TARGET
OPT=-O3 -march=native
PROF=-g -p
REPORT=
# -frecursive is the gfortran spelling of ifort's -recursive; there is no
# -heap-arrays equivalent (gfortran heap-allocates large temporaries anyway).
PARALLEL=-fopenmp -frecursive
# -fdec-static enables the DEC "AUTOMATIC" attribute used in corrk.f and
# ropprmulti_corrk.f.  -r8 -> -fdefault-real-8 -fdefault-double-8,
# -132 -> -ffixed-line-length-132, -traceback -> -fbacktrace.
# -fno-align-commons matches ifort's COMMON layout: ifort packs COMMON blocks
# with no padding, gfortran pads by default.  COMMON/CLOUDY/ leads with
# CHARACTER(30) AEROSOLMODEL followed by REAL*8, so the default would shift
# every subsequent member.  (ifort's -warn noalign silences the same warning.)
DEBUG_OPTS=
# -std=legacy accepts the F66/F77 spellings this model is written in
# (arithmetic IF, shared DO termination labels, DO ending on a statement
# other than CONTINUE) without reporting each one.  They compile
# correctly either way, and silencing them is what makes real warnings
# visible: without it 162 notices buried a genuine out-of-bounds read.
FFLAGS = -g $(OPT) $(REPORT) $(PARALLEL) $(DEBUG_OPTS) -fdec-static \
         -fno-align-commons -std=legacy \
         -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fbacktrace

else

FC = ifort
OPT=-O3 -xHost
#OPT='-Ofast -xHost'
PROF=-g -p
REPORT=-warn noalign
PARALLEL=-fopenmp -recursive -heap-arrays 0
# DEBUG=-debug extended
DEBUG_OPTS=
FFLAGS = -g $(OPT) $(REPORT) $(PARALLEL) $(DEBUG_OPTS) -r8 -132 -traceback

endif

# Unreachable from the program entry point and dropped from the build:
#   the pre-rradiation.f radiation chain hanging off NIKOSRAD
#   (cnikos, cradsw, csw, cswtt, cirrad, cdedd, ccldtrn, ccolamt,
#   co3interp -- its callers cmorc.f/cmorc_old.f are not built),
#   plus double-gray-ropprmulti, rinterpol, chanal1, cinisurf and
#   cinterp.  Their sources are in legacy/ (see legacy/README.md).
OBJS = datapaths.o cloudmods.o corrk.o tracerclds.o dyeadv.o ropprmulti_corrk.o radiative_transfer_picket_fence.o cmltri_nopg.o cbalanc.o cblayer.o cblsurf.o ccalndr.o ccbadj.o ccbcon.o cconvec.o ccubm.o ccudif.o cdanalv.o cdgrmlt.o cdifuse.o cdlsgcr.o cdryadj.o cdstep.o cenergy.o chanal.o chanalv.o chexp.o chexp1.o chexpv.o cictrac.o cinibal.o cinigau.o ciniphys.o ciniqs.o cinires.o ciniresij.o ciniset.o cinisi.o cinisp.o cinistr.o inisimprad.o inivarparam.o cinital.o clgndre.o clscrn.o cltd.o cltddia.o clti.o cltidt.o cmascor.o cmatinv.o cmgrmlt.o cnoise.o cpqsat.o cpvcr.o csetres.o csettee.o csetzt.o csfct.o cspdel2.o cspop.o csurfm.o ctbal.o ctstep.o cvdiff.o cwrsps.o cxsect.o cfft991.o cssum.o csdot.o cicamax.o cgwtlt.o cset99.o cqreig.o csgetrf.o csgetri.o csgemm.o crpassm.o cqpassm.o chessen.o cqrt.o cxerbla.o cilaenv.o csgetf2.o cslaswp.o cstrsm.o cstrtri.o csgemv.o csswap.o clsame.o cisamax.o csscal.o csger.o cstrmm.o cstrti2.o cstrmv.o cranf.o xsect2.o xsect3.o filecopy.o finalorb.o Binary.o rradiation.o rcalc_radheat.o rradsub.o rsetuprad_simple.o rradtran.o roppr1.o rtwostr.o radd.o rnewflux1.o rmakeclouds.o ropprmulti.o cloud_properties_set_up.o radiative_transfer_corrk.o corrk_setup.o roppr1_corrk.o calc_v_fall.o

igcm3_nopg: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<

# params.i sizes the COMMON blocks and arrays (NTRAC, NL, MG, ...), so
# every object must be rebuilt when it changes; a partial rebuild mixes
# layouts and corrupts the COMMONs silently.  mkparams.sh runs on every
# make but rewrites params.i only when the configuration changes, so an
# unchanged one recompiles nothing.
$(OBJS): params.i nclouds.inc $(wildcard *.h)

# Fortran modules: a user needs the module's .mod built first and must
# be recompiled when it changes.  Rather than track which file uses
# which module, every other object depends on all of them.
MODOBJS = datapaths.o cloudmods.o corrk.o tracerclds.o dyeadv.o
corrk.o: datapaths.o
dyeadv.o: tracerclds.o
$(filter-out $(MODOBJS),$(OBJS)): $(MODOBJS)

params.i: FORCE
	./mkparams.sh $(RES) $(NL) $(RT) $(NDYE) $(NTGRID) $(NPGRID)

buildinfo.inc: params.i ;

cinital.o: buildinfo.inc

.PHONY: FORCE clean distclean
FORCE:

clean:
	rm -f *.o *.mod igcm3_nopg

# Also forget the configuration: the next make builds the defaults.
distclean: clean
	rm -f params.i buildinfo.inc build.cfg
