This README file is intended to serve as a user's guide to the iGCM code. 
It is a compendium of notes that should be ammended and updated by users as they see fit.

-----------------------------------------------
0. BACKGROUND:
-----------------------------------------------
University of Reading's website with some old documentation of the original code: 
http://www.met.reading.ac.uk/~mike/dyn_models/igcm/

This Intermediate General Circulation Model (igcm) code was modified by Menou and Rauscher,
as described in Menou & Rauscher (2009). It is currently a psuedo-spectral model that solves
the primitive equations and impliments a double-gray radiation scheme.  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-----------------------------------------------
1. HOW TO COMPILE & RUN THE CODE
-----------------------------------------------
~~~~~~~~~~~~~
1.1 Params.i
~~~~~~~~~~~~~ 
To run in parallel, first you must compile it with the proper input file-- params.i

params.i includes information on the model resolution.  It essentially provides the compiler
with the information required to prepare arrays of the correct sizes. 
An example below shows a case for a grid resolution referred to as T42 L30:
C     T42 L30  full sphere
      PARAMETER(NN=42,MM=42,NHEM=2 ,NL=30,MOCT=1,MG=128,JG=32,NWJ2=462
     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)

Note the first line is commented out.
The T42 refers to the spectral resolution.  More information on these resolutions and how
the relate to physical resolutions can be found at:
https://climatedataguide.ucar.edu/climate-model-evaluation/common-spectral-model-grid-resolutions

L30 Refers to the number of vertical layers in the model, corresponding to the parameter NL=30.
MG refers to the number of longitudes, and JG refers to th enumbe rof latitudes per hemisphere.
NHEM is the number of hemispheres (...two is usually good).  

So you need a file named params.i that has a readable PARAMETER().  
This input information can be supplied by uncommenting-out the desired lines of a much larger
file that includes a commented-out list of many possible resolutions.  
Alternatively you can rename the full list params.i_full and copy the desired lines to a seperate
file called params.full, or write a script to generate it as needed.

~~~~~~~~~~~~~~~~~
1.2  Compiling
~~~~~~~~~~~~~~~~~
The details of the compilation are set in the the compile_nopg file.
The file is a list of commands or flags specific to the compiler being used.  For example:

#!/bin/bash -x

rm -f igcm3_nopg 

OPT='-O3 -xHost'
#PROF='-g -p'
REPORT='-warn noalign -openmp-report2'
PARALLEL='-openmp'
#OTHER_OPTS='-ftrapuv'
ifort -g $OPT $REPORT $PARALLEL $OTHER_OPTS -r8 -132 -traceback -o igcm3_nopg cmltri_nopg.f cbalanc.f cblayer.f cblsurf.f ccalndr.f ccbadj.f ccbcon.f ccldtrn.f ccolamt.f cconvec.f ccubm.f ccudif.f cdanalv.f cdedd.f cdgrmlt.f cdifuse.f cdlsgcr.f cdryadj.f cdstep.f cenergy.f chanal.f chanal1.f chanalv.f chexp.f chexp1.f chexpv.f cictrac.f cinibal.f cinigau.f ciniphys.f ciniqs.f cinires.f ciniresij.f ciniset.f cinisi.f cinisp.f cinistr.f cinisurf.f inisimprad.f inivarparam.f cinital.f cinterp.f cirrad.f clgndre.f clscrn.f cltd.f cltddia.f clti.f cltidt.f cmascor.f cmatinv.f cmgrmlt.f cmorc.f  cnikos.f cnoise.f co3interp.f cpqsat.f cpvcr.f cradsw.f csetres.f csettee.f csetzt.f csfct.f cspdel2.f cspop.f csurfm.f csw.f cswtt.f ctbal.f ctstep.f cvdiff.f cwrsps.f cxsect.f cfft991.f cssum.f csdot.f cicamax.f  cgwtlt.f cset99.f cqreig.f csgetrf.f csgetri.f csgemm.f crpassm.f cqpassm.f chessen.f cqrt.f cxerbla.f cilaenv.f csgetf2.f cslaswp.f cstrsm.f cstrtri.f csgemv.f csswap.f clsame.f cisamax.f csscal.f csger.f cstrmm.f cstrti2.f cstrmv.f cranf.f xsect2.f xsect3.f filecopy.f finalorb.f Binary.f


The first line prints this flags to the screen so the user knows what is happening. The second 
removes the old executable (compilers turn code into executables).  And so on... the only one 
you need to worry about is the Parallel flag. If you want the code to run in parallel, this flag
must be set.  If not, comment it out by placing a # in the first column. More instructions on 
running the code in parallel are below.

TO COMPILE the code, just use the following command:
./compile_nopg

A new executable igcm3_nopg should be produced.

~~~~~~~~~~~~~~~~~~~~~~~
1.3  RUNNING THE CODE
~~~~~~~~~~~~~~~~~~~~~~~
Once the executable has been created, you can copy it to a nice clean directory for executing.
In that directory, you will need another input file named fort.7.
================
1.3.1.  fort.7
================
The fort.7 file contains all the other information specific to your simulation beyond resolution.
The complete list is long, but some entries of possible interest include:
 &INPPL   ;
 GA = 17.7, ; Gravitational accration m/s^2
 GASCON  =   3197.7,    ; specific gas constant  J/kg/ K
 RADEA   =   2.94e7,    ;  planetary radius  , meters
 AKAP    =   0.2857,   ; gas constant divided by specific  R/Cp   2/7
 WW      =   1.08e-4,  ; rotation rate  radians/sec
 P0      =   3.0e6,   ; bottom boundary pressure, pascals, though it can vary
 RV      =   461.510000000000     ,; g /kg of condensible  don't care
 CPD     =   12619.5,   ;    Gas constant at constant pressure for dry air (unused)

 &INPRN
 KRUN    =           5.76e6,   ; # of total times steps in the run
 BEGDAY  =  0.000000000000000E+000,   ;  zero for new run, non-zero for restart run to be explained
 TSPD    =   800.,                       ; time steps per day (1 day is one full sideral  2 !pi/ ww)
 KITS    =           3,                  ;3 for initiation, but set this to zero for restart (check)
 PNU     =  2.000000000000000E-002,      ; Robert time filter, leave as is
 TDISS   =  0.40306,                    ;  Very important.  Time scale for hyper dissipation on small scale.  Units of planet day.
 NDEL    =           8,                    ; exponent on the hypderdissipation
 T0      = 30*400.00,        ;Number of layers * initial temperatures

 RNTAPE  =  0.000000000000000E+000,      ; in practice, not used
 KOUNTH  =           50000,            ; timesteps between outputting or saving data: History, some multiple of timestep per day, depends on rotation rate, usually 100 days
 KOUNTR  =           50000,    ; timesteps between outputting or saving data: Restart, some multiple of timestep per day, depends on rotation rate, usually 100 days
 KOUNTP  =           800,      ; set equal to time steps per day
 LGPO    = 30*F,    ; integer multiple of number vertical levels
 LSPO    = 30*F,    ;integer multiple of number vertical levels
 OOM_IN        = 5.0,     ; orders of magnitude in pressure we are covering, defines the top of model ... if 100 bar and 5 oom then top is ~1mbar
NTSTEP_IN    = 25,     ; how often you run the radiative transfer to redifine the tendencies,
FBASEFLUX    = 0.354,      ; Bottom boundary flux condition, (internal flux upwards) W/m^2
OPACIR_POWERLAW    = 0,       ; height dependence of opacity ... integer values are quick.  Non-integer values are make things MUCH, MUCH, SLOWER and should be avoided if at all possible!
OPACIR_REFPRES    = 3.e5,     ; reference pressure for IR absorption coefficient in radiative transfer  K=___(p/pref) ^ exponent (above)
SOLC_IN        = 2051.25,      ; Top boundary condition, W/m^2
TOAALB        = 0.3,            ; top of the atmosphere albedo
PORB        = 450.,         ; ratio between the orbital period and rotational period, ~ for tidally locked, >1 for other cases
OBLIQ        = 0.0,             ; Obliquity in degrees
ECCEN        = 0.0            ; eccentricity but needs work and validation
NLCR    =           29,   ; controls which layers are allowed to be convective, set to number of vertical levels -1
TFRC    = 30*0.,    ; Rayleigh friction at each pressure level  time scale (nlayers * 0)
 RESTIM  = 30*1e30,   ; Timescale for newtonian forcing ;;so large, we are not using it
RESTTT = 341.122, 352.315, 366.800, 384.957, 406.947, 432.636,   ;important; Temperatures initial conditions from top down to bottom, analytically determined
      461.524, 492.674, 524.632, 555.375, 582.425, 603.309,
      616.545, 622.770, 624.683, 625.102, 625.334, 625.652,
      626.119, 626.801, 627.798, 629.253, 631.372, 634.442,
      638.870, 645.208, 654.185, 666.727, 683.945, 707.078,
 REDTEP = 30*0.,   ;array nlayers * # (amplitude at each level if newtonian forcing is used)
 LDIUR                = T,   ; if true , diurnally averaged flux pattern at the top of the atmopshere, if false, explicity day night pattern
ABSSW1  =  8.14e-4,            ; absorption coefficient for optical cm^2/gram
ABSLW1  =  1.49e-2,         ; IR absorption coefficient at the reference pressure, cm^2/gram

...
Some things to note:

a) Several parameters are dependent on the number of layers.  In the above case, 30 layers are expected.
The user must ensure that the number in fort.7 matches the number specified in the compiled params.i file,
or else the code will crash upon running.  A sneaky one is NLCR, which should be set to NL-1.

b) Some parameters are factors of the number of time steps per day.  

c) For numerical stability, the time step, spatial/spectral resolution, and hyper-dissipation must be balanced. This takes some adjusting.

d) Choosing a integer value for the opacity power law will greatly reduce the run time.  Non-integer values should be avoided if possible.

  
=============================
1.3.2. Executing (in serial)
=============================
With the executable produced and the fort.7 file supplied, the code can be run using the following command:
./igcm_nopg

Now that's fine, but it's the bear minimum.  It would be better to include 'nohup' to ensure that the
run is not interrupted when you logout, and 'nice' to share the computer with other processes.  Also, it 
is often convenient to create a file that records how long the run took, and to make the entire process
run in the background by ending everything with an ampersand. The resulting command would be:

'nohup nice /usr/bin/time -p -o timing ./igcm3_nopg &'

Since that's a bit much, I recommend editing your .cshrc file to include an alias similar to:

alias runigcm 'nohup nice /usr/bin/time -p -o timing ./igcm3_nopg &'

This permits you to simply type 'runigcm' as a shortcut for the full command.

=============================
1.3.3. RUNNING IN PARALLEL
=============================
Running the code in parallel requires a few extra steps. 

1) make sure your compile_nopg code has the -openmp flag set.
2) Then if you running tshell, set the number of threads you wish to use with the following command:

setenv OMP_NUM_THREADS #
where # is the number of threads desired

For bash:
export OMP_NUM_THREADS=#

And you can see that it worked by using the 'printenv' command.

3) If this resolution is high and/or the number of threads is great, the stack memory is overwhelmed.
This can be rememdied by increasing the default 10mb limit on stack memory as follows:

limit stacksize unlimited

4) Exectue as usual
e.g. nohup nice /usr/bin/time -p -o timing ./igcm3_nopg &

5) To confirm that the calculations are using multiple threads, you can use the 'top' command. You'll note that the %CPU for the process >100% for a parallel process. 'top -H' will display each parallel process separately, and you can simply count them to see how many threads are currently being used.

=============================
1.3.3. USING RESTART TO RUN
=============================
There are situations where you might want to restart a run, without starting it again from the very beginning (perhaps the atmosphere hasn't reached a statistical equilibrium, perhaps the computer was turned off).  You can do this by using the restart mechanism:

1) No need to recompile.  Copy fort.7, fort.11, and fort.13 from the run you want to restart.

2) In fort.7 change the following parameters, as described:
  a) LRSTRT=T [this is the logical switch to perform a restart run]
 [b) LSHORT=F (this is the logical to use initial short timesteps and should already be =F, but double-check)]
 [c) LRESTIJ=T (related to initialization, should already be =T, but double-check)]
  d) KITS=0 [this is the number of short initial timesteps to use]
  e) BEGDAY: set this to the last day a restart record was written
     (You can figure this out by checking fort.26 for the day on which the run ended, and then calculating the most recent record from the greatest integer multiple of KOUNTR/TSPD, where KOUNTR gives the number of timesteps between each write of the restart record.)
  f) KRUN: set to the number of additional timesteps you want the simulation to run, starting from BEGDAY