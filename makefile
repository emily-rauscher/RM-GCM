# Compiler selection.  Plain "make" builds with ifort
# using exactly the flags this model has always used.
#   make                     -> ifort   (default, unchanged)
#   make COMPILER=gfortran   -> gfortran
COMPILER ?= ifort

ifeq ($(COMPILER),gfortran)

FC = gfortran
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
FFLAGS = -g $(OPT) $(REPORT) $(PARALLEL) $(DEBUG_OPTS) -fdec-static \
         -fno-align-commons \
         -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fbacktrace

else

FC = ifort
OPT=-O3 -xHost
#OPT='-Ofast -xHost'
PROF=-g -p
REPORT=-warn noalign
PARALLEL=-fopenmp -recursive -heap-arrays 0
# DEBUG=-debug extended
DEBUG_OPTS=''
FFLAGS = -g $(OPT) $(REPORT) $(PARALLEL) $(DEBUG_OPTS) -r8 -132 -traceback

endif

OBJS = corrk.o ropprmulti_corrk.o double-gray-ropprmulti.o radiative_transfer_picket_fence.o cmltri_nopg.o cbalanc.o cblayer.o cblsurf.o ccalndr.o ccbadj.o ccbcon.o ccldtrn.o ccolamt.o cconvec.o ccubm.o ccudif.o cdanalv.o cdedd.o cdgrmlt.o cdifuse.o cdlsgcr.o cdryadj.o cdstep.o cenergy.o chanal.o chanal1.o chanalv.o chexp.o chexp1.o chexpv.o cictrac.o cinibal.o cinigau.o ciniphys.o ciniqs.o cinires.o ciniresij.o ciniset.o cinisi.o cinisp.o cinistr.o cinisurf.o inisimprad.o inivarparam.o cinital.o cinterp.o cirrad.o clgndre.o clscrn.o cltd.o cltddia.o clti.o cltidt.o cmascor.o cmatinv.o cmgrmlt.o cnikos.o cnoise.o co3interp.o cpqsat.o cpvcr.o cradsw.o csetres.o csettee.o csetzt.o csfct.o cspdel2.o cspop.o csurfm.o csw.o cswtt.o ctbal.o ctstep.o cvdiff.o cwrsps.o cxsect.o cfft991.o cssum.o csdot.o cicamax.o cgwtlt.o cset99.o cqreig.o csgetrf.o csgetri.o csgemm.o crpassm.o cqpassm.o chessen.o cqrt.o cxerbla.o cilaenv.o csgetf2.o cslaswp.o cstrsm.o cstrtri.o csgemv.o csswap.o clsame.o cisamax.o csscal.o csger.o cstrmm.o cstrti2.o cstrmv.o cranf.o xsect2.o xsect3.o filecopy.o finalorb.o Binary.o rradiation.o rcalc_radheat.o rradsub.o rsetuprad_simple.o rradtran.o rinterpol.o roppr1.o rtwostr.o radd.o rnewflux1.o rmakeclouds.o ropprmulti.o cloud_properties_set_up.o radiative_transfer_corrk.o corrk_setup.o roppr1_corrk.o calc_v_fall.o

igcm3_nopg: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod igcm3_nopg
