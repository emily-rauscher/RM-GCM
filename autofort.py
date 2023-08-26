#!/usr/bin/env python

### This turns the data in eplanetpars.xlsx into sets of fort.7s for each planet
### It now makes one fort.7 depending on your specific parameters
### By changing "destpath", you can make a separate set of fort.7s in their own
### directory (i.e. one directory with 200 day runs, another with 1000 day runs,
### and a third with TSPD=6000 for UHJs)

### NOTE: currently set up to ignore eccentricity unless it's greater than 0.2

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil

df = pd.read_excel("./eplanetpars.xlsx", header=0)
destpath = "./" # IMPORTANT! Where to put all these fort.7 files, needs to end in a forward slash

### Most common choices:
names = ['WASP-76b']
num_forts = 4

oom       = [7,9, 11, 13]  #number of order of magnitudes spanned from surface pressure
nlay      = [50] #number of vertical layers in the model (set above)
tspd      = [4800.] # 4800 for HJs, 6000 for UHJs
tdiss     = [0.005] # 0.005 +- a factor of 5
aerlayers = [10., 20., 30., 40.] # how thick can clouds get
numdays   = [200] # how many orbits
ntstepin  = [8] # 8 works for cooler planets, but some will need 4 instead
sponge = [0.005] # sponge layer timescale (in units of planet day). The top layer of the atmosphere gets
                # sponge, then the strength decreases by a factor of ten for each of the next
                # two layers (only the top three layers are affected). 0.005 is a good starting point
toaalb = [0.05]
CLOUDS = ['T']
hazetype = ['soot']
rayscatlam = [['0.2, 0.5, 0.8,']]
metallicity = [0.0]
ABSLW = [1.0e-2]
ABSSW = [4.0e-3]
PICKET_FENCE_CLOUDS = ['T']
lbdrag = ['T']
lrstrt = ['F']

def dimcheck(variable, num_forts):
    if len(variable) != num_forts:
        if len(variable) != 1:
            raise Exception('Number of input variables not matching, check ', variable)
        elif len(variable) > num_forts:
            raise Exception('Number of input variables not matching, check ', variable)
        assign = variable[0]
        for x in range(1, num_forts):
            variable.append(assign)

def remove_folder(path):
    # check if folder exists
    if os.path.exists(path):
         # remove if exists
         shutil.rmtree(path)
    else:
         pass



### Creating folders & sub-folders to hold the fort.7s

dimcheck(oom, num_forts)
dimcheck(nlay, num_forts)
dimcheck(tspd, num_forts)
dimcheck(tdiss, num_forts)
dimcheck(aerlayers, num_forts)
dimcheck(numdays, num_forts)
dimcheck(ntstepin, num_forts)
dimcheck(sponge, num_forts)
dimcheck(toaalb, num_forts)
dimcheck(CLOUDS, num_forts)
dimcheck(hazetype, num_forts)
dimcheck(rayscatlam, num_forts)
dimcheck(metallicity, num_forts)
dimcheck(ABSLW, num_forts)
dimcheck(ABSSW, num_forts)
dimcheck(PICKET_FENCE_CLOUDS, num_forts)
dimcheck(lbdrag, num_forts)
dimcheck(lrstrt, num_forts)


for name in names:
    remove_folder(destpath + name)
    os.mkdir(destpath + name)
    for x in range(0, num_forts):
        # Loop through all planets, creating fort.7s using numbers in eplanetpars.xslx
        radjup= df["R_p (R_J)"][df['Name']==name].values[0]
        massjup= df["M_p (M_J)"][df['Name']==name].values[0]
        pdays= df["P (days)"][df['Name']==name].values[0]
        rstar= df["R* (R_sun)"][df['Name']==name].values[0]
        Dau= df["a (au)"][df['Name']==name].values[0]
        tstar= df["T* (K)"][df['Name']==name].values[0]
        eccen = df["e"][df['Name']==name].values[0]
        if type(eccen)==float:
            if eccen < 0.2:
                eccen = 0
        #metal = df["[Fe/H] (dex)"][df['Name']==pname].values[0]

        # calculate values necessary for fort.7
        radmet=radjup* (7.14e7)
        masskg=massjup*1.899e27
        bigG=6.67e-11 #m^3 kg^-1 s^-2
        surfgrav=(bigG*masskg)/ ((radmet)**2)
        ww=(2*np.pi)/(pdays*24*3600)
        sigmasb=5.671e-5
        dcm=Dau*1.496e13
        rstarcm=rstar* 6.96e10
        #note Tirr doesn't go in fort.7 but is needed for setting up 1D profile
        Tirrcalc=((rstarcm/dcm)**.5)*(tstar)
        solc_in=(Tirrcalc**4)*sigmasb #erg/cm^2
        solc_in=solc_in*100*100*(1e-7)
        rstarm=rstar*6.96e8
        orbsep=Dau*1.496e11
        #Teq also doesn't go into fort.7, but a useful number
        Teq = (tstar*((rstarm/(2*orbsep))**(1/2)))

        #setting up your tp profile
        f=.375 # is redistribution number, often .375 but can go from 0-1
        # f=.23
        #change the opacities below to force a thermal inversion
        kth=(ABSLW[x]) *1 #ABSLW
        kvis=(ABSSW[x]) *1 #ABSSW
        gamma=(kvis/kth)

        grav=surfgrav #m/s
        surfp=100 #surface pressure in bars, usually 100
        tirr=Tirrcalc #calculated above
        tint=550. #internal temp, sets the behavior of the bottom part of the model

        longpost=False #set True to file.write off temperatures in easy format for pasting into fort.7

        def avgsurface(tau,tirr,tint): #EQ 29. CORRECT ONE!
            first=.75*(tint**4)*((2./3.)+ tau)
            second=.75*(tirr**4)*f*( (2./3.)+ (1/(gamma * np.sqrt(3))) + ((gamma/np.sqrt(3)) - (1/(gamma*np.sqrt(3))))*np.exp(-gamma*tau*np.sqrt(3)))
            final=first+second
            return final**.25
        def makeinitialtprofs(nlay,oom,surfp,tirr,tint):
            sigma=np.empty([nlay])*0.0  #calculating array of pressure values based on # of layers
            if oom>0: #setting up pressure values
                stp=-1.0*oom/nlay
                sigma[nlay-1]=10.**(stp/2.)
                for n in range(nlay-2,-1,-1):
                    sigma[n]=sigma[n+1]*10.**(stp)
            p_BAR=sigma*surfp
            taus=np.copy(p_BAR)
            taus=(10*p_BAR* 1000)/(943.4) #converting to optical depths
            initialts=np.copy(taus)
            for i in range(len(taus)):
                initialts[i]=avgsurface(taus[i],tirr,tint)
            return initialts, p_BAR

        temps,press=makeinitialtprofs(nlay[x],oom[x],surfp,tirr,tint)

        if longpost==True:
            for i in range(len(temps)):
                file.write (round(temps[i],2))
            #file.write (temps)
        temparray=np.linspace(min(temps), max(temps),100)
        vissphere=((6.67e-5)*grav)/kvis
        irsphere=((6.67e-5)*grav)/kth


        if CLOUDS[x] == 'T':
            molef = [1.23e-7,0.0,0.0,0.0,4.4e-7,3.26e-5,1.745e-5,9.56e-9,0.0,0.0,1.99e-6,7.83e-8,1.385e-6]
            for y in range(0, len(molef)):
                molef[x] *= (10 ** metallicity[x])
            cloudsyn = 'nucleation-limited clouds'
            case = 'cloud_'
        elif CLOUDS[x] == 'F':
            molef = [0,0,0,0,0,0,0,0,0,0,0,0,0]
            cloudsyn = 'no clouds'
            case = 'clear_'
        else:
            raise Exception('CLOUDS not called correctly, must be "T" or "F"\n')
        if lbdrag[x] == 'T':
            magyn = 'with magnetic drag'
            case += 'yamag'
        elif lbdrag[x] == 'F':
            magyn = 'no magnetic drag'
            case += 'nomag'
        else:
            raise Exception('Magnetic drag (lbdrag) not called correctly, must be "T" or "F"\n')


        file = open(destpath+name+ "/version" + str(x) + "fort.7", 'x')
        file.write("&COMMENT\n")
        file.write(" THECOMMENT= '{}, {}, {}'\n".format(name, cloudsyn, magyn))
        file.write(" /\n")
        file.write(" &INPPL      \n")
        file.write(" GA      = {},\n".format(np.round(surfgrav, decimals=2))) # Gravitational acceleration, m/s^2
        file.write(" GASCON  = 3523.,\n") # specific gas constant  J/kg/ K
        file.write(" RADEA   = {},\n".format(radmet)) # planetary radius  , meters
        file.write(" AKAP    = 0.2860000000000,\n") # gas constant divided by specific  R/Cp   2/7
        file.write(" WW      = {},\n".format(ww)) # rotation rate  radians/sec
        file.write(" P0      = 1.0E+7     ,\n") # Bottom boundary pressure (pascals), this is 100 bars
        file.write(" RV      = 461.510000000000     ,\n") # g/kg of condensible air (what does that mean?)
        file.write(" CLATNT  = 2500000.00000000     \n")
        file.write(" /\n")
        file.write(" &INPRN\n")
        file.write(" KRUN     = {},\n".format(tspd[x]*numdays[x])) # total number of steps
        file.write(" BEGDAY   = 0.000000000000000E+000,\n") # zero for new run, non-zero for restart run
        file.write(" TSPD     = {},\n".format(tspd[x])) # time steps per day (1 day is one full sideral  2 !pi/ ww)
        file.write(" KITS     = 3,\n") # 3 for initiation, but set this to zero for restart (check)
        file.write(" PNU      = 2.000000000000000E-002,\n")
        file.write(" TDISS    = {},\n".format(tdiss[x])) # Very important.  Time scale for hyper dissipation on small scale.  Units of planet day.
        file.write(" LFLUX    = T,\n")
        file.write(" BEGDOY   = 0.000000000000000E+000,\n")
        file.write(" NDEL     = 8,\n") # exponent on the hypderdissipation
        file.write(" T0       = {}*1700.,\n".format(nlay[x])) # Number of layers * initial temperatures
        file.write(" LRSTRT   = {},\n".format(lrstrt[x]))
        file.write(" LSTRETCH = F,\n")
        file.write(" LSHORT   = F,\n")
        file.write(" LTVEC    = F,\n")
        file.write(" LBALAN   = F,\n")
        file.write(" LRESTIJ  = T,\n")
        file.write(" LCLIM    = F,\n")
        file.write(" LPERPET  = T,\n")
        file.write(" L22L     = F,\n")
        file.write(" LOROG    = F,\n")
        file.write(" LCSFCT   = F,\n")
        file.write(" KOLOUR   = 0,\n")
        file.write(" LNOISE   = T,\n")
        file.write(" LMASCOR  = T,\n")
        file.write(" LMASOLD  = T,\n")
        file.write(" LMASPRT  = F\n")
        file.write(" /\n")
        file.write(" &INPOP\n")
        file.write(" RNTAPE  = 0.000000000000000E+000,\n") # in practice, not used
        file.write(" KOUNTH  = {},\n".format(tspd[x]*10)) # timesteps between outputting or saving data: History, some multiple of timestep per day, depends on rotation rate, usually 100 days
        file.write(" KOUNTR  = {},\n".format(tspd[x]*10)) # timesteps between outputting or saving data: Restart, some multiple of timestep per day, depends on rotation rate, usually 100 days
        file.write(" KOUNTP  = {},\n".format(tspd[x])) # set equal to time steps per day
        file.write(" KOUNTE  = 0,\n")
        file.write(" NCOEFF  = 0,\n")
        file.write(" NLAT    = 4,\n")
        file.write(" LGPO    = {}*F,\n".format(nlay[x])) # integer multiple of number vertical levels
        file.write(" LSPO    = {}*F,\n".format(nlay[x])) # integer multiple of number vertical levels
        file.write(" RNTAPO  = 0.000000000000000E+000,\n")
        file.write(" NTRACO  = 0,\n")
        file.write(" LSHIST  = T,\n")
        file.write(" LMINIH  = T\n")
        file.write(" /\n")
        file.write("  &INVARPARAM\n")
        file.write(" OOM_IN      = {},\n".format(oom[x])) # Pressure OOM
        file.write(" LPLOTMAP 	= F,\n")
        file.write(" NLPLOTMAP_IN	= 10,\n")
        file.write(" RFCOEFF_IN	= 0.0e-0, \n")
        file.write(" NTSTEP_IN	= {},\n".format(ntstepin[x])) # how often you run the radiative transfer to redifine the tendencies
        file.write(" NSKIP_IN	= 0, \n")
        file.write(" BOTRELAXTIME	= 1.0, \n")
        file.write(" FBASEFLUX	= 3543.,\n") # Bottom boundary flux condition, (internal flux upwards) W/m^2
        file.write(" FORCE1DDAYS	= 0.0e0, \n")
        file.write(" OPACIR_POWERLAW	= 0.0,\n")
        file.write(" OPACIR_REFPRES	= 1.0e5, \n") # reference pressure for IR absorption coefficient in radiative transfer  K=___(p/pref) ^ exponent (above)
        file.write(" SOLC_IN		= {},\n".format(solc_in)) # Top boundary condition, W/m^2
        file.write(" TOAALB		= {}, \n".format(toaalb[x])) # top of the atmosphere albedo
        file.write(" PORB		= 1.0, \n") # ratio between the orbital period and rotational period, ~1 for tidally locked, >1 for other cases
        file.write(" OBLIQ		= 0.0, \n") # Obliquity in degrees
        file.write(" ECCEN		= {}\n".format(eccen)) # eccentricity
        file.write(" /\n")
        file.write(" &INMAG\n")
        file.write("  LBDRAG = {},\n".format(lbdrag[x]))
        file.write("  BFIELD = 3.,\n")
        file.write("  TDRAG_MIN = 0.005,\n")
        file.write("  RAMPUP = 25.\n")
        file.write("  /\n")
        file.write(" &INBINVAL\n")
        file.write(" LBIN = F,\n")
        file.write(" PORBST=10.8,\n")
        file.write(" ECCPL=0.400,\n")
        file.write(" ECCST=0.0234,\n")
        file.write(" SMAPL=0.989,\n")
        file.write(" SMAST=0.0836,\n")
        file.write(" STMASS1=1.043,\n")
        file.write(" STMASS2=0.362,\n")
        file.write(" STRAD1=0.964\n")
        file.write(" STRAD2=0.3506,\n")
        file.write(" STTEMP1=5636\n")
        file.write(" STTEMP2=3357\n")
        file.write(" /\n")
        file.write("  &INPHYS\n")
        file.write("  LBL     = F,\n")
        file.write("  LVD     = T,\n")
        file.write("  LCR     = T,\n")
        file.write("  LLR     = F,\n")
        file.write("  LRD     = T,\n")
        file.write("  LCUBM   = F,\n")
        file.write("  LCBADJ  = F,\n")
        file.write("  CD      =  1.000000000000000E-003,\n")
        file.write("  BLRH    =   100.000000000000     ,\n")
        file.write("  AKVV    =   1.0E-20,\n")
        file.write("  AKTV    =   1.0E-20,\n")
        file.write("  AKQV    =   1.0E-20,\n")
        file.write("  AKTC    =   10.0000000000000     ,\n")
        file.write("  AKQC    =   10.0000000000000     ,\n")
        file.write("  NLCR    =           {},\n".format(nlay[x]-1)) # controls which layers are allowed to be convective, set to number of vertical levels -1
        file.write("  NCUTOP  =           5,\n")
        file.write("  FW      =  1.000000000000000E-006,\n")
        file.write("  FRADC   =   1.25000000000000     ,\n")
        file.write("  CURHM   =   200.000000000000     ,\n")
        file.write("  CUBMT   =   3.00000000000000     ,\n")
        file.write("  CBADJT  =   3.00000000000000     ,\n")
        file.write("  CBADJP  =  -30.0000000000000     ,\n")
        file.write("  LSL     = F,\n")
        file.write("  LOC     = F,\n")
        file.write("  LNOICE  = F,\n")
        file.write("  LOLDBL  = F,\n")
        file.write("  LCOND   = T,\n")
        file.write("  LNNSK   = T,\n")
        file.write("  ITSLL   =           6,\n")
        file.write("  ITSLO   =           6\n")
        file.write("  /\n")
        file.write("  &INPRSIJ\n")
        file.write("  TFRC    = {},{},{},{}*0.\n".format(sponge[x], sponge[x]*10,sponge[x]*100,nlay[x]-3)) # Controls sponge layers
        # Rayleigh friction at each pressure level  time scale (nlayers * 0)
        file.write(" RESTIM  = {}*1e30,\n".format(nlay[x])) # timescale for newtonian forcing ;;so large, we are not using it
        file.write(" RESTTT =  {} ,\n".format(str(temps)[1:-1])) # important; Temperatures initial conditions from top down to bottom, analytically determined
        file.write(" REDTEP  = {}*0.,\n".format(nlay[x])) # array nlayers * # (amplitude at each level if newtonian forcing is used)
        file.write(" DTNS    =  0.0000000E+00,\n")
        file.write(" DTEP    =   300.00000    ,\n")
        file.write(" ALR     =  2.0E-06,\n")
        file.write(" DTTRP   =   10.000000    ,\n")
        file.write(" ZTROP   =   2.0e7    ,\n")
        file.write(" TGR     =   2100.0    ,\n")
        file.write(" YRLEN   =  0.0000000E+00\n")
        file.write(" /\n")
        file.write(" &INSIMPRAD\n")
        file.write(" LLOGPLEV        = T,\n")
        file.write(" LFLUXDIAG       = T,\n")
        file.write(" L1DZENITH       = F,\n")
        file.write(" LDIUR           = F,\n")
        file.write(" JSKIPLON        = 1,\n")
        file.write(" JSKIPLAT        = 1,\n")
        file.write(" DOSWRAD         = T,\n")
        file.write(" DOLWRAD         = T,\n")
        file.write(" LWSCAT          = T,\n")
        file.write(" FLXLIMDIF       = F,\n")
        file.write(" SURFEMIS        = 1.,\n")
        file.write(" RAYSCAT         = F,\n")
        file.write(" RAYSCATLAM      = {},\n".format(rayscatlam[x]))
        file.write(" AEROSOLS        = T,\n")
        file.write(" ABSSW           = {}, \n".format(ABSSW[x])) # absorption coefficient for optical cm^2/gram
        file.write(" ABSLW           = {},\n".format(ABSLW[x])) # IR absorption coefficient at the reference pressure, cm^2/gram
        file.write(" ALBSW           = 0.0,\n")
        file.write(" NEWTB           = 0,\n")
        file.write(" NEWTE           = 0,\n")
        file.write(" with_TiO_and_VO = 1,\n")
        file.write(" picket_fence_optical_depths = T\n")
        file.write(" /\n")
        file.write(" &INCLOUDY\n")
        file.write(" AEROSOLMODEL = 'Global',\n")
        file.write(" AEROSOLCOMP  = 'standard',\n")
        file.write(" HAZETYPE     = {},\n".format(hazetype[x])) #hazetype, so soot, tholin, or clear
        file.write(" MTLX         = .1,\n")
        file.write(" METALLICITY  = {},\n".format(metallicity[x]))
        file.write(" HAZES = F,\n")
        file.write(" PICKET_FENCE_CLOUDS = {},\n".format(PICKET_FENCE_CLOUDS[x]))
        file.write(" MOLEF        = {},\n".format(str(molef)[1:-1])) # Clouds
        file.write(" AERLAYERS    = {},\n".format(aerlayers[x])) # how many layers a cloud can extend through
        file.write(" AERTOTTAU    = 100.,\n")
        file.write(" CLOUDBASE    = 10.,\n")
        file.write(" CLOUDTOP     = 0.0001,\n")
        file.write(" AERHFRAC     = 1.,\n")
        file.write(" PI0AERSW     = 1.0,\n")
        file.write(" ASYMSW       = 0.376,\n")
        file.write(" EXTFACTLW    = 0.01,\n")
        file.write(" PI0AERLW     = 0.0452,\n")
        file.write(" ASYMLW       = 0.005896,\n")
        file.write(" DELTASCALE   = T,\n")
        file.write(" SIG_AREA     = 25.,\n")
        file.write(" PHI_LON      = 65.\n")
        file.write(" /\n")
        file.write(" &INQ\n")
        file.write(" LRH     = T,\n")
        file.write(" RH      = {}*1e-5    ,\n".format(nlay[x]))
        file.write(" LNSURF  = F\n")
        file.write(" /")
        file.close()
