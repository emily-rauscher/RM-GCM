FC = ifort
OPT=-O3 -xHost
#OPT='-Ofast -xHost'
PROF=-g -p
REPORT=-warn noalign
#PARALLEL='-fopenmp'
OTHER_OPTS=-debug extended

FFLAGS = -g $(OPT) $(REPORT) $(PARALLEL) $(OTHER_OPTS) -r8 -132 -traceback

OBJS = double-gray-ropprmulti.o radiative_transfer_picket_fence.o cmltri_nopg.o cbalanc.o cblayer.o cblsurf.o ccalndr.o ccbadj.o ccbcon.o ccldtrn.o ccolamt.o cconvec.o ccubm.o ccudif.o cdanalv.o cdedd.o cdgrmlt.o cdifuse.o cdlsgcr.o cdryadj.o cdstep.o cenergy.o chanal.o chanal1.o chanalv.o chexp.o chexp1.o chexpv.o cictrac.o cinibal.o cinigau.o ciniphys.o ciniqs.o cinires.o ciniresij.o ciniset.o cinisi.o cinisp.o cinistr.o cinisurf.o inisimprad.o inivarparam.o cinital.o cinterp.o cirrad.o clgndre.o clscrn.o cltd.o cltddia.o clti.o cltidt.o cmascor.o cmatinv.o cmgrmlt.o cnikos.o cnoise.o co3interp.o cpqsat.o cpvcr.o cradsw.o csetres.o csettee.o csetzt.o csfct.o cspdel2.o cspop.o csurfm.o csw.o cswtt.o ctbal.o ctstep.o cvdiff.o cwrsps.o cxsect.o cfft991.o cssum.o csdot.o cicamax.o cgwtlt.o cset99.o cqreig.o csgetrf.o csgetri.o csgemm.o crpassm.o cqpassm.o chessen.o cqrt.o cxerbla.o cilaenv.o csgetf2.o cslaswp.o cstrsm.o cstrtri.o csgemv.o csswap.o clsame.o cisamax.o csscal.o csger.o cstrmm.o cstrti2.o cstrmv.o cranf.o xsect2.o xsect3.o filecopy.o finalorb.o Binary.o rradiation.o rcalc_radheat.o rradsub.o rsetuprad_simple.o rradtran.o rinterpol.o roppr1.o rtwostr.o radd.o rnewflux1.o rmakeclouds.o ropprmulti.o cloud_properties_set_up.o

igcm3_nopg: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod igcm3_nopg