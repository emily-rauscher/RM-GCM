C**********************************************************               
C             SUBROUTINE NIKOSRAD                                         
C**********************************************************               
       SUBROUTINE NIKOSRAD(PR,T,H2O,O3,alat,HTLW,HTSW,DOY,CF,IC,          
     $     RFLUXES,SWALB,ALON,KOUNT,ITSPD)

       include 'params.i'
C                                                                       
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 
      
       LOGICAL LPLOTMAP

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,ABSMMRLW,
     & JSKIPLON,JSKIPLAT,NSWMODEL,NLWMODEL,ABSSW1,ABSSTRAT,PRMIN,ALBSW1,
     & ABSSW2,SCATSW2,ASYMSW2,ABSLW1,NEWTB,NEWTE, with_TiO_and_VO, opacity_method
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR
       CHARACTER(len=6) opacity_method

      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC                        
     +              ,RNTAPO                                               
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR                          
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY                          
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG                 
     +              ,LSPO(NL),LGPO(NL)                                    
     $              ,LSHIST,LMINIH

      LOGICAL LSHIST,LMINIH                                               
      LOGICAL LSPO,LGPO    

! need this for the logical switch for if binary or not.
      COMMON/BINVAL/PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2,STRAD1,STRAD2,STTEMP1,STTEMP2,LBIN

      LOGICAL LBIN

C                 -----------                                             
C  Driving program for the extended atmospheric radiation scheme of       
C              Lacis and Hansen and Nikos Christidis.                     
C                                                                         
C                                                                         
C***********************************************************************  
C*                                                                     *  
C*                         P A R R A Y                                 *  
C*                                                                     *  
C***********************************************************************  
C-----------------------------------------------------------------------  
C In this part of the code the user can set the dimensions of arrays      
C used for the calculations, by setting the values in the paramter        
C statements                                                              
C-----------------------------------------------------------------------  
                                                                          
      INTEGER MXGAS,MXLEV,MXBAND,MXCL                                     
                                                                          
      PARAMETER(MXGAS=8)     ! Maximum number of gases                    
                                                                          
      PARAMETER(MXLEV=NL)     ! Maximum number of levels in the            
                             ! atmosphere (not including the surface)     
                                                                          
      PARAMETER(MXBAND=9)    ! Maximum number of spectral bands (not      
                             ! including the whole spectrum, 0-3000cm-1)  
                                                                          
      PARAMETER(MXCL=3)      ! Maximum number of cloud types              
                                                                          
C-----------------------------------------------------------------------  
C                                                                         
C  Switches for long wave radiation scheme                                
C                                                                         
      COMMON/RADLW/VMRCO2,VMRCH4,VMRN2O,VMRHALO,GAS(MXGAS),LLBLM          
      LOGICAL LLBLM                                                       
      INTEGER GAS                                                         
      REAL VMRCO2,VMRCH4,VMRN2O,VMRHALO                                   
C                                                                         
C     Sets basic constants, especially those needed for array dimensions  
C                                                                         
      PARAMETER(MH=2,PI=3.14159265359,PI2=2.0*PI                          
     +,NNP=NN+1,MGPP=MG+2,JGP=JG+1,JGG=JG*NHEM,JGGP=JGG+1,MJP=NWJ2+NWJ2   
     +,NLM=NL-1,NLP=NL+1,NLPP=NL+2,NLA=NL+3,NLB=NL+4,NL2=NL*NL            
     +,IDA=(MG+MG+MG)/2+1,IDB=NWJ2*NL,IDC=IDB+IDB,IDD=MGPP*NL             
     +,IDE=NL2*NN,IDF=NCRAY*(MG+1),IDG=JG*NL,IDH=JG*MG                    
     +,IDI=NNP/2,IDJ=IDI*IDI,IDK=NL*IDI,IDL=MGPP/2,IDM=NNP/2,IDN=IDM*NL   
     +,NWW=1+(MM-1)/MOCT)                                                 
      PARAMETER(IGA=NWJ2*NHEM,IGB=IDB*NHEM,IGC=MGPP*NHEM,IGD=IDD*NHEM     
     +,IGG=IDG*NHEM,IGL=IDL*NHEM,IGM=IDM*NHEM,IGN=IDN*NHEM                
     +,IGO=IGA+IGA,IGP=IGB+IGB,NFTWG=(5+NTRAC)*NL+3                       
     +,NFTGW=(6+3*NTRAC)*NL+2,NFTGD=(3+NTRAC)*NL,NLTR=NL*NTRAC)           
C     Number of 2D (surface) output fields. This value is                 
C     Doubled due to averaged and instantaneous fields.                   
      PARAMETER (N2DFLD=21,NGRPAD=N2DFLD*2*IGC)                           
C                                                                         
C                                                                         
C     Basic planetary parameters for run plus information about           
C     vertical grid structure                                             
C                                                                         
C     Note that RD and GASCON are identical and CPD is set from RD,AKAP.  
      COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)                
     +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)     
     +              ,TKP(NL),C(NL2),SQH(NNP)                              
     +              ,MF,MFP,JZF,NF                                    
     +              ,AKAP,GA,GASCON,RADEA,WW,PFAC,EZ,AIOCT             
     +              ,RD,RV,CPD,CLATNT                                     
     +              ,P0,LRSTRT,LSHORT,LTVEC,LSTRETCH                         
     +              ,LFLUX                                                
     +              ,LBALAN,LRESTIJ                                       
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
     +              ,LNOISE,NFP                                               
      COMPLEX EZ,AIOCT                                                    
      LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ                 
     +       ,LFLUX,LNOISE                                                
     +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                           
C                                                                         
      PARAMETER (NRLEV=MXLEV+1)                                           
      REAL PR(NRLEV),T(NRLEV),H2O(NRLEV),O3(NRLEV),alat                   
      real htlw(nrlev),htsw(nrlev),swalb,DOY,alon                         
      REAL CF(4,2)     !cloud fractions,lwps                              
c clouds in order of convective,high,mid,low                              
c ie convective has index of 1 with low cloud index of 4                  
      INTEGER IC(4,2)   ! cloud positions

      REAL RFLUXES(2,2,2)

C      LOGICAL LDIUR    
********************************************************************      
C Morcrette                                                               
      REAL FLS(nrlev+1)                                                   
      REAL SWDP(NRLEV),pclfr(nrlev),pqlwp(nrlev)                          
C Nikos                                                                   
      INTEGER CKD, IFIRST                                                 
      DATA IFIRST/1/                                                      
      REAL PRES(0:MXLEV),TEMP(0:MXLEV),MMR(MXGAS,0:MXLEV)                 
      INTEGER CL(0:MXCL)                                                  
      REAL ACL(0:MXCL),LWCCL(MXCL),BASCL(0:MXCL),TOPCL(0:MXCL)            
C Pressure at flux levels                                                 
      REAL PFLUX(0:MXLEV)                                                 
C ER Modif
      REAL TLEV(0:MXLEV)                                                 
      REAL SSLON,SSLAT          ! ER: previously delta_zenith angle, now substellar location
      REAL DLENGTH  ! ER: half-length of solar day

C------------------                                                       
C Output Variables                                                        
C------------------                                                       
                                                                          
      REAL FUP(0:MXLEV),FDWN(0:MXLEV),FNET(0:MXLEV)                       
      REAL PDP(MXLEV)                                                     
C********************************************************************     
C                                                                         
C long wave radiation namelist                                            
      NAMELIST/INRADLW/LLBLM,GAS,VMRCO2,VMRCH4,VMRN2O,                    
     :             VMRHALO                                                
      LLBLM=.TRUE.


      write(*,*) 'Isaac doesnt think this radiation scheme should be used'

C ER modif for non-synchronous orbit
      IF (PORB.NE.0) THEN 
         SSLON=(1./PORB-1.)*KOUNT*360./ITSPD
C         SSLON=ALON-SSLON
         SSLON=MOD(SSLON,360.)
      ELSE
         SSLON=0.  ! substellar longitude
C         SSLON=ALON  !local longitudinal angle to star, in degrees
      ENDIF

C ER modif for non-zero obliquity
      IF (OBLIQ.EQ.0) THEN
         SSLAT=0.  ! substellar latitude
C        SSLAT=ALAT  !local latitudinal angle to star, in degrees
         DLENGTH=PI/2.
      ELSE
         SSLAT=ASIN(SIN(OBLIQ*PI/180.)
     +        *SIN(PI2*KOUNT/ITSPD/PORB))*180./PI
C        SSLAT=ALAT-SSLAT
         IF (SSLAT.GT.0) THEN
            IF (ALAT.GT.90.-SSLAT) THEN
               DLENGTH=PI
            ELSEIF (ALAT.LT.-90.+SSLAT) THEN
               DLENGTH=0.
            ELSE
               DLENGTH=ACOS(-1.*TAN(ALAT/360.*PI2)*TAN(SSLAT/360.*PI2))
            ENDIF
         ELSEIF (ALAT.LT.-90.-SSLAT) THEN
            DLENGTH=PI
         ELSEIF (ALAT.GT.90+SSLAT) THEN
            DLENGTH=0.
         ELSE
            DLENGTH=ACOS(-1.*TAN(ALAT/360.*PI2)*TAN(SSLAT/360.*PI2))
         ENDIF
      ENDIF
C                                                                         
C    Setup for LW scheme                                                  
      NLEV=NRLEV-1                                                        
      CKD=1           ! water vapour continium                            
      IF (CKD.EQ.1.AND.IFIRST.EQ.1)THEN                                   
          WRITE (2,*)'Water Vapour Continuum On'                          
      ELSEIF(CKD.EQ.0.AND.IFIRST.EQ.1)THEN                                
          WRITE (2,*)'Water Vapour Continuum Off'                         
      ENDIF                                                               
      IF (NL.NE.MXLEV)THEN                                                
        WRITE(*,*)'JOB ABORTED! NL= ',NL,' MXLEV= ',MXLEV                 
        CALL ABORT                                                        
      ENDIF                                                               
C set up default values.                                                  
      IF (IFIRST.EQ.1) THEN                                               
       IFIRST=0                                                           

C                                                                         
       GAS(1)=1                                                           
       GAS(2)=1                                                           
       GAS(3)=1                                                           
       GAS(4)=0                                                           
       GAS(5)=0                                                           
       GAS(6)=0                                                           
       GAS(7)=0                                                           
       GAS(8)=0                                                           
C run time is 7% quicker with gasses 4-8 switched off (over 12 days)      
       VMRCO2=358.0E-6                                                    
       VMRCH4=1.72E-6                                                     
       VMRN2O=312.0E-9                                                    
       VMRHALO=7.97E-12                                                   
C READ namelists, overwrite defaults and write them out                   
!!       READ(7,INRADLW)                                                    
       WRITE(2,INRADLW)                                                   
      ENDIF                                                               
      DO LHT=1,NRLEV                                                      
C MMR is mass mixing ratio of gases                                       
c       PCLD(LHT)=0.0      ! set SW cloud opotical depth                  
        PCLFR(LHT)=0.0      ! set SW cloud fraction                       
        PQLWP(LHT)=0.0      ! set SW cloud lwp                            
        MMR(1,LHT-1)= ABSMMRLW !H2O(LHT)                                             
        MMR(2,LHT-1)=VMRCO2*(44.011/28.964)                               
        MMR(3,LHT-1)=O3(LHT)                                              
        MMR(4,LHT-1)=VMRCH4*(16.043/28.964)                               
        MMR(5,LHT-1)=VMRN2O*(44.014/28.964)                               
        MMR(6,LHT-1)=MMR(2,LHT-1)                                         
        MMR(7,LHT-1)=MMR(3,LHT-1)                                         
        MMR(8,LHT-1)=VMRHALO*(125.5/28.964)                               
        PRES(LHT-1)=PR(LHT)                                               
        TEMP(LHT-1)=T(LHT)                                                
!        write(*,*) LHT, H2O(LHT)
      ENDDO                                                               
C Calculate pressures at flux levels                                      
      IF(LLOGPLEV) THEN  
!     KM: Refine layer thicknesses for log P grid
C     Redefine Delta-Ps of layers as averages of log(DP), not averages of DP
C     PFLUX unchanged since it is not used by H2OFLUX routine
         PFLUX(0)=PRES(0)/100.0
         DO ILAY=1,NLEV-1
            PFLUX(ILAY)=(PRES(ILAY)+PRES(ILAY+1))/200.0
         END DO
         PFLUX(NLEV)=PRES(NLEV)/200.0
         PDP(1)=EXP(0.5*(LOG(PRES(0)-PRES(1))+LOG(PRES(0)-PRES(2))))
         DO ILAY=2,NLEV-1
            PDP(ILAY)=EXP(0.5*(LOG(PRES(ILAY-1)-PRES(ILAY)) 
     &           +LOG(PRES(ILAY)-PRES(ILAY+1))))
         END DO
         PDP(NLEV)=EXP(0.5*(LOG(PRES(NLEV-1))+LOG(PRES(NLEV))))
      ELSE
         PFLUX(0)=PRES(0)/100.0                                              
         DO ILAY=1,NLEV-1                                                    
            PFLUX(ILAY)=(PRES(ILAY)+PRES(ILAY+1))/200.0                      
         END DO                                                              
         PFLUX(NLEV)=PRES(NLEV)/200.0                                        
         PDP(1)=0.5*(2*PRES(0)-PRES(1)-PRES(2))                              
         DO ILAY=2,NLEV-1                                                    
            PDP(ILAY)=0.5*(PRES(ILAY-1)-PRES(ILAY+1))                        
         END DO                                                               
         PDP(NLEV)=0.5*PRES(NLEV-1)                                          
      ENDIF

      GO TO 456                 ! skip clouds                                        
C clouds                                                                  
      DO ICL=1,3     !low=1,3=high                                        
         CL(ICL)=0                                                        
         ICLIC=5-ICL !for ic 4=low,2=high,1=convective                    
                                                                          
C         IF (IC(ICLIC,1).GT.1) THEN                                      
         IF (CF(ICLIC,1).GT.0.0) THEN                                     
C Cloud present                                                           
           CL(ICL)=1                                                      
C single level cloud                                                      
           LEV=IC(ICLIC,1)                                                
           BASCL(ICL)=PR(IC(ICLIC,1))                                     
           TOPCL(ICL)=PR(IC(ICLIC,1))                                     
           ACL(ICL)=CF(ICLIC,1)                                           
C cloud LWP is 1% super saturation in g/m^2                               
           CF(ICLIC,2)=0.01*SVP(T(LEV))/PR(LEV)*                          
     $          622*PDP(LEV-1)/GA                                         
c           CF(ICLIC,2)=CF(ICLIC,2)*REAL(ICL)                             
           LWCCL(ICL)=CF(ICLIC,2)/1000.0  ! g- Kg convert                 
c sw cloud optical depth                                                  
c           PCLD(IC(ICLIC,1))=1.5*CF(ICLIC,2)/10.0*ACL(ICL)               
           PCLFR(IC(ICLIC,1))=CF(ICLIC,1)                                 
           PQLWP(IC(ICLIC,1))=CF(ICLIC,2)                                 
        ENDIF                                                             
      ENDDO                                                               
c convective cloud                                                        
      CL(0)=0                                                             
C      IF (IC(1,1).GT.1) THEN                                             
       IF (CF(1,1).gt.0.0) THEN                                           
         CL(0)=1                                                          
         ACL(0)=CF(1,1)                                                   
C         LEV=(IC(1,1)+IC(1,2))/2                                         
         LEV=IC(1,1)                                                      
C cloud LWP is 2% super saturation in g/m^2                               
         CF(1,2)=0.01*SVP(T(LEV))/PR(LEV)*                                
     $          622*PDP(LEV-1)/GA                                         
         ALWP=CF(1,2)                                                     
         BASCL(0)=PR(IC(1,1))                                             
         TOPCL(0)=PR(IC(1,2))                                             
c quarter cloud except in bottom level                                    
         DO ICL=IC(1,1)+1,IC(1,2)                                         
c            PCLD(ICL)=1.5* LWP/10.0* ACL(0)                              
            PCLFR(ICL)=0.25*CF(1,1)                                       
            PQLWP(ICL)=ALWP                                               
         ENDDO                                                            
            PCLFR(IC(1,1))=CF(1,1)                                        
            PQLWP(IC(1,1))=ALWP                                           
       ENDIF                                                              
                                                                          
 456   CONTINUE                                                           
C                                                                         
C call nikos LW scheme                                                    
      CALL IRRAD(GAS,CKD,NLEV,PRES,PFLUX,TEMP,MMR,CL,ACL,LWCCL,           
     $          BASCL,TOPCL,FUP,FDWN,FNET,PDP,LLBLM,GA,TLEV)                   
                                                                          
!     ER hack for flux diffusion at large tau (replace FNETs)
      DO LHT=1,MXLEV
!     ABSLW1 in cm^2/g, GA in m/s^2, PFLUX in mbar (=1e3 g/cm/s^2)
         TAUCONST=(ABSLW1/GA/100.)  ! units of cm s^2/g
C         TAU=TAUCONST*PFLUX(LHT)*1.e3
C     REFPRES in Pa (=10 g/cm/s^2)
c         write(*,*) mod(2.0,1.0)
c            IF (mod(OPACIR_POWERLAW,1.0).eq.0.0) THEN
c              OPACIR_POWERLAW=int(OPACIR_POWERLAW) 
c        MODIF: Replacing the following line:
c         TAU=TAUCONST*PFLUX(LHT)*1.E3
c     &        *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)**(OPACIR_POWERLAW)
c        with a case specific exponent solution--MTR
           IF (OPACIR_POWERLAW.EQ.0) THEN
             TAU=TAUCONST*PFLUX(LHT)*1.E3
           ELSEIF (OPACIR_POWERLAW.EQ.1) THEN
         TAU=TAUCONST*PFLUX(LHT)*1.E3
     &        *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)
           ELSEIF (OPACIR_POWERLAW.EQ.2) THEN
         TAU=TAUCONST*PFLUX(LHT)*1.E3
     &   *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)
     &   *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)
           ELSEIF (OPACIR_POWERLAW.EQ.3) THEN
         TAU=TAUCONST*PFLUX(LHT)*1.E3
     &   *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)
     &   *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)
     &   *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)
           ELSE 
         TAU=TAUCONST*PFLUX(LHT)*1.E3
     &        *(1.e2*PFLUX(LHT)/OPACIR_REFPRES)**(OPACIR_POWERLAW)
           ENDIF
c          END MODIF MTR

         IF (TAU.GT.1.) THEN 
            TLIMIT=2.*LOG(1.e-2)/1.66
     &           /(10.**(-1.*OOM_IN/NL)-10.**(OOM_IN/NL))
         IF (TAU.GT.TLIMIT) THEN 
            IF (LHT.EQ.1) THEN
               GRAD=(TEMP(LHT+1)-TEMP(LHT))/(PRES(LHT+1)-PRES(LHT))*100.  !for mbar
               
c            MODIF: Same as above. Replacig ** with specific cases--MTR
                  IF (OPACIR_POWERLAW.EQ.0) THEN
                     EXPCORR=TAUCONST*1.E6*0.5
     &                *(2*PRES(LHT-1)-PRES(LHT)-PRES(LHT+1))
                  ELSEIF (OPACIR_POWERLAW.EQ.1) THEN
                     EXPCORR=TAUCONST*1.E6*0.5
     &              *(2*PRES(LHT-1)-PRES(LHT)-PRES(LHT+1))
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
                  ELSEIF (OPACIR_POWERLAW.EQ.2) THEN
                     EXPCORR=TAUCONST*1.E6*0.5
     &              *(2*PRES(LHT-1)-PRES(LHT)-PRES(LHT+1))
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
                  ELSEIF (OPACIR_POWERLAW.EQ.3) THEN
                     EXPCORR=TAUCONST*1.E6*0.5
     &              *(2*PRES(LHT-1)-PRES(LHT)-PRES(LHT+1))
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
                  ELSE
                     EXPCORR=TAUCONST*1.E6*0.5
     &              *(2*PRES(LHT-1)-PRES(LHT)-PRES(LHT+1))
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)**OPACIR_POWERLAW
                  ENDIF
c                 END MODIF--MTR
            ELSEIF (LHT.EQ.MXLEV) THEN
               GRAD=(TLEV(LHT)-TEMP(LHT))/(PFLUX(LHT)-PRES(LHT)/100.)
c                MODIF as above--MTR
c               EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
c     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)**OPACIR_POWERLAW
                IF (OPACIR_POWERLAW.EQ.0) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
                ELSEIF (OPACIR_POWERLAW.EQ.1) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
                ELSEIF (OPACIR_POWERLAW.EQ.2) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
                ELSEIF (OPACIR_POWERLAW.EQ.3) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
                ELSE
               EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)**OPACIR_POWERLAW
                ENDIF
            ELSE
               GRAD=(TEMP(LHT+1)-TEMP(LHT))/(PRES(LHT+1)-PRES(LHT))*100.
C               EXPCORR=TAUCONST*1.E6*0.5*(PRES(LHT)-PRES(LHT+1))
C     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)**OPACIR_POWERLAW
                  IF (OPACIR_POWERLAW.EQ.0) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*(PRES(LHT)-PRES(LHT+1))
                  ELSEIF (OPACIR_POWERLAW.EQ.1) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES) 
                  ELSEIF (OPACIR_POWERLAW.EQ.2) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
                  ELSEIF (OPACIR_POWERLAW.EQ.3) THEN
                   EXPCORR=TAUCONST*1.E6*0.5*PRES(LHT)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)                     
                  ELSE
               EXPCORR=TAUCONST*1.E6*0.5*(PRES(LHT)-PRES(LHT+1))
     &              *(1.E2*PFLUX(LHT)/OPACIR_REFPRES)**OPACIR_POWERLAW               
                  ENDIF

            ENDIF  ! units: K/mbar = 1e-3 K cm s^2 /g
!           boltzmann is 5.6704e-8 W/m^2 K^-4; want flux in W/m^2
C            FNETDIFF=3.0242e-10*(GRAD/TAUCONST)*((TLEV(LHT))**3) 
C           ONE MORE MODIF, AS ABOVE--MTR
C            FNETDIFF=3.0242E-10*GRAD*(TLEV(LHT))**3/TAUCONST
C     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)**OPACIR_POWERLAW
                IF (OPACIR_POWERLAW.EQ.0) THEN
                 FNETDIFF=3.0242E-10*GRAD*(TLEV(LHT))**3/TAUCONST
                ELSEIF (OPACIR_POWERLAW.EQ.1) THEN
                 FNETDIFF=3.0242E-10*GRAD*(TLEV(LHT))**3/TAUCONST
     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)
                ELSEIF (OPACIR_POWERLAW.EQ.2) THEN
                 FNETDIFF=3.0242E-10*GRAD*(TLEV(LHT))**3/TAUCONST
     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)
     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)
                ELSEIF (OPACIR_POWERLAW.EQ.3) THEN
                 FNETDIFF=3.0242E-10*GRAD*(TLEV(LHT))**3/TAUCONST
     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)
     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)
     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)
                ELSE
                FNETDIFF=3.0242E-10*GRAD*(TLEV(LHT))**3/TAUCONST
     &           *(OPACIR_REFPRES/PFLUX(LHT)/1.E2)**OPACIR_POWERLAW
                ENDIF
!           linear implementation, from tau=1 to 100
C            IF (TAU.LT.5) THEN
C               TFAC=(1./4.)*TAU-(1./4.)
C            ELSE
C               TFAC=1.
C            ENDIF
C            TFAC=(1.-(1.005)**(1.-TAU))
C            TFAC=1.-EXP(-1.66*(TAU-1.)/100.)
            TFAC=1.-EXP(-1.66*EXPCORR/43.7)
            FNET(LHT)=(1.0-TFAC)*FNET(LHT) + TFAC*FNETDIFF
C            FNET(LHT)=FNETDIFF
         ENDIF
         ENDIF
      END DO
!     IF (LHT.EQ.0) THEN FNET(LHT)=FBASEFLUX+FLS(1)/(1.0-SWALB)+FDWN(0)
      FNET(0)=FBASEFLUX+FLS(1)/(1.0-SWALB)
C      FNET(0)=FBASEFLUX+FLS(1)/(1.0-SWALB)+FDWN(0)

C Setup SW code                                                           
      IF (LBIN) THEN
        call BinaryFlux(SOLC,KOUNT,ITSPD)
        SOLC=SOLC*(1.0-TOAALB)
!        WRITE(88,*) SOLC
      ELSE
        SOLC=SOLC_IN*(1.0-TOAALB)
      ENDIF
          
!      LDIUR=.TRUE.   !Diurnally averaged if false                        
C      LDIUR=.FALSE.   !Diurnally averaged if false
C                                                                         
C       YCLOCK=3.14159      !time of day in radians                        
!       YCLOCK= (DOY-REAL(INT(DOY))*2.*3.14159                            
C                                                                         
C      CALL SOLANG(LDIUR,DOY,YCLOCK,ALAT,ALON,AMU0,RDAYL,CDISSEM)          
C         
C     globally averaged solar constant, vertical rays
      AMU0=1.0            
      PSOL=SOLC/4.             
      IF(.NOT.L1DZENITH) THEN 
         DDAY=FORCE1DDAYS 
CC       DDAY is rampup time from uniform to full heating
C        (1D for DDAYs, then linear to full in another DDAYs)
         IF(DAY.GT.DDAY) THEN
            DFAC=MIN(1.0,(DAY - DDAY)/DDAY)
            IF(.NOT.LDIUR) THEN
CC Modif for hemispheric forcing of Hot Jupiters
               AMU0=(1.0-DFAC)*AMU0 
     &              +DFAC*MAX(0.0,SIN(ALAT/360.*PI2)*SIN(SSLAT/360.*PI2)
     &                           +COS(ALAT/360.*PI2)*COS(SSLAT/360.*PI2)
     &                           *COS((ALON-SSLON)/360.*PI2))
               PSOL=(1.0-DFAC)*PSOL + DFAC*SOLC 
            ELSE
CC ER modif for diurnal forcing (a la Liu & Schneider 2010)
C               PSOL=(1.0-DFAC)*PSOL + DFAC*SOLC/PI*COS(ALAT/360.*PI2)
C ER modif for non-zero obliquity
               PSOL=(1.0-DFAC)*PSOL
     &              +DFAC*SOLC/PI*
     &              (SIN(ALAT/360.*PI2)*SIN(SSLAT/360.*PI2)*DLENGTH
     &              +COS(ALAT/360.*PI2)*COS(SSLAT/360.*PI2)*SIN(DLENGTH)
     &              )
            ENDIF
         ENDIF
      ENDIF

! KM For  planet with variable irradiation flux
!         IF(DAY.GT.(2*DDAY)) THEN
!            ECC=0.3
!         PSOL=SOLC*(1.+ECC*SIN((DAY-2*DDAY)/360.*PI2))  
!         ENDIF

C      ELSE
c call to L+H                                                              
C      ZSCT=SOLC*CDISSEM   ! Solar const * earth-sun distance              
C      PSOL=ZSCT*RDAYL     ! * fractional day length                       
C      ENDIF

c CALL SW scheme                                                          
      ZCARDI=MMR(2,1)                                                     
      CALL RADSW(T,H2O,O3,PCLFR,PQLWP,AMU0,ZCARDI,PSOL,SWDP,PFLUX,SWALB,FLS,GA,GASCON)

!KM: New Boundary condition for local energy conservation (verify, if SWALB>0)
!      FBASE=5.7  ! Jupiter-like: 5.7 W m^-2 from interior
!      FUP(0)=FDWN(0)+FLS(1)/(1.0-SWALB)+FBASE
!      FNET(0)=FLS(1)/(1.0-SWALB)+FBASE
                                                                          
c WORKS OUT FLUXES                                                        
C SW                                                                      
      RFLUXES(1,1,1)=PSOL*AMU0   ! SW down top                            
      RFLUXES(1,1,2)=FLS(1)/(1.0-SWALB)   ! SW down bottom                
      RFLUXES(1,2,1)=RFLUXES(1,1,1)-FLS(NRLEV+1)  ! SW up top             
      RFLUXES(1,2,2)=RFLUXES(1,1,2)*SWALB   ! SW up bottom                
C LW                                                                      
      RFLUXES(2,1,1)=FDWN(NLEV)   ! LW down top                           
      RFLUXES(2,1,2)=FDWN(0)       ! LW down bottom                       
      RFLUXES(2,2,1)=FUP(NLEV)       ! LW up top                          
      RFLUXES(2,2,2)=FUP(0)   ! LW up bottom

      write(*,*) 'ERROR! OLD CODE IN CNIKOS IS CALLED WITH WRONG NUMBER OF BINS'
      stop
C OUTPUT                                                                  
                                                                          
      GRCP=(GA/CPD)*24*3600 ! g/Cp *24*3600                               
                                                                          
C Calculate heating rates                                                 
C ER Modif: HTSW(LHT) should be FLS(LHT)-FLS(LHT-1), changed to that
                                                                          
      DO LHT=2,NRLEV         ! from bottom-most full sigma level to top                
C         HTSW(LHT)=GRCP*(FLS(LHT+1)-FLS(LHT))/SWDP(LHT)                   
         HTSW(LHT)=GRCP*(FLS(LHT)-FLS(LHT-1))/SWDP(LHT)                   
         HTLW(LHT)=-GRCP*(FNET(LHT-1)-FNET(LHT-2))/SWDP(LHT)              
      END DO                                                              
      IF (.NOT.(NEWTB.EQ.NEWTE)) THEN
         DO LHT=NEWTB,NEWTE
            LEV=NRLEV+1-LHT
            HTSW(LEV)=0.
            HTLW(LEV)=0.
         ENDDO
      ENDIF

!!      LFLUXDIAG=.TRUE. !! Switch on or off diagnostics in fort.63 below
C ER Modif: only write fort.63 every kountp timesteps
C     (kountp-1 b/c in cmorc nikos called when mod(kount,ntstep).eq.1)
      IF ((LFLUXDIAG).AND.(KOUTP.EQ.KOUNTP-1)) THEN
         WRITE(63,*) 'LATITUDE, LONGITUDE:',ALAT,ALON   
         WRITE(63,2010)'UPWARD FLUX (Wm-2)','DOWNWARD FLUX (Wm-2)',         
     $        'NET FLUX (Wm-2)', ' SW NET FLUX (Wm-2)'             
 2010    FORMAT(11X,A18,5X,A20,5X,A15,5x,A24)                               

C     ER Modif: output pressures in bar instead of mbar
         DO ILAY=NLEV,0,-1                                                  
            WRITE(63,2013)PFLUX(ILAY)/1.e3,FUP(ILAY),                            
     $           FDWN(ILAY),FNET(ILAY),FLS(ILAY+1)                          
     $           
 2013       FORMAT(2X,F12.3,5X,E12.5,10X,E12.5,10X,E12.5,10X,E12.5)         
         END DO                                                             
         
         WRITE(63,2023)'PRESSURE (bar)','HEATING RATES: SW (K/DAY)'          
     $        ,'LW (K/DAY)'                                                 
 2023    FORMAT(1X,A18,1X,A30,x,A10)                                        
         DO LHT=NRLEV,1,-1                                                  
            WRITE(63,2020) PR(LHT)/100.0/1.e3,HTSW(LHT),HTLW(LHT)                
 2020       FORMAT(5X,F10.3,X,E12.5,X,E12.5)                                
         END DO               
         WRITE(63,*)
      ENDIF                                                                    
                                                                          
      RETURN                                                              
      END                                                                 

C ER modif: commenting this all out, since it's unused              
c      SUBROUTINE SOLANG (LDIUR,DOY,YCLOCK,ALAT,alon,                      
c     :     AMU0,RDAYL,CDISSEM)                                            
cC**********************************************************               
cC             SUBROUTINE SOLANG                                           
cC**********************************************************               
cC Pm.M. deF   27-1-98                                                     
cC inputs                                                                  
cC LDIUR  ! Logical for diurnal average                                    
cC DOY    ! Julian day of year                                             
cC YCLOCK !                                                                
cC ALAT   ! Latitude in degrees                                            
cC ALON   ! Longitude in degrees                                           
cC outputs                                                                 
cC AMU0   ! Cosine of solar zenith angle                                   
cC RDAYL  ! Fractional day length                                          
cC CDISSEM  ! Suns relative distance as a fraction                         
c                                                                          
cC**** *SOLANG* - FOR SOLAR ZENITH ANGLE AND RELATIVE DAYLENGTH.           
cC                                                                         
cC     PURPOSE.                                                            
cC     --------                                                            
cC                                                                         
cC          THIS ROUTINE GIVES DIFFERENT RESULTS DEPENDING ON A LOGICAL    
cC     SWITCH. IF LDIUR IS TRUE ONE OBTAINS ACTUAL SOLAR ZENITH ANGLES     
cC     AND VALUES OF ONE OR ZERO DEPENDING ON THE SIGN OF THE FORMER. IF   
cC     LDIUR IS FALSE ONE GETS THE SAME ANSWERS AT ALL POINTS, I.E. MEAN   
cC     VALUE OF THE DAYTIME SOLAR ZENITH ANGLE AND RELATIVE LENGTH OF      
cC     THE DAY.                                                            
cC                                                                         
cC   LDIUR .... true: sun at time, false: diurnally averaged               
cC                                                                         
cC     ----------                                                          
cC                                                                         
c      REAL DOY,ALAT,RDAYL,AMU0,YCLOCK,ALON,XLON                           
cC                                                                         
c      REAL XLAT,YTIME,YEARL,API,ZC1YT,ZS1YT,ZC2YT,ZS2YT,CDISSEM           
c      INTEGER JDAY                                                        
c                                                                          
c      REAL ZMU0(128),ZRDAYL(128)                                          
c      LOGICAL LDIUR,LO                                                    
c      REAL ZCDIS(5),ZCEQT(5),ZCDEC(5)                                     
c      DATA ZCDIS/+1.000110,+0.034221,+0.001280,+0.000719,+0.000077/       
c      DATA CRAE/+0.1277E-02/                                              
c      DATA ZCDEC/+0.006918,-0.399912,+0.070257,-0.006758,+0.000907/       
c      DATA ZCEQT/+0.000075,+0.001868,-0.032077,-0.014615,-0.040849/       
c                                                                          
c                                                                          
c      YEARL=360.0   ! 360 DAY YEAR                                        
c      API=2.0*ASIN(1.0)                                                   
c      XLAT=API*alat/180.0                                                 
c      XLON=API*ALON/180.0                                                 
c      JDAY=INT(DOY)                                                       
c      YTIME = (REAL(JDAY-1)+YCLOCK/2./API)/YEARL*2.*API                   
cC                                                                         
cC diurnal cycle part                                                      
cC                                                                         
cC*    COMPUTATIONAL CONSTANTS.                                            
cC     ------------- ----------                                            
cC                                                                         
c      ZC1YT=COS(YTIME)                                                    
c      ZS1YT=SIN(YTIME)                                                    
c      ZC2YT=ZC1YT**2-ZS1YT**2                                             
c      ZS2YT=2.*ZS1YT*ZC1YT                                                
c      CDISSEM=ZCDIS(1)+ZCDIS(2)*ZC1YT+ZCDIS(3)*ZS1YT+ZCDIS(4)*ZC2YT       
c     *       +ZCDIS(5)*ZS2YT                                              
cC                                                                         
c      ZCRAE=CRAE*(CRAE+2.)                                                
cC     ------------------------------------------------------------------  
cC*         2.     SOLAR ANGLE AND OZONE/AEROSOL PARAMETERS COMPUTATIONS.  
cC                 ----- ----- --- ------------- ---------- -------------  
cC                                                                         
c 200  CONTINUE                                                            
cC                                                                         
cC*         2.1     INTRODUCE THE LATITUDE DEPENDENCY.                     
c      ZSIN = SIN(XLAT)                                                    
c      ZSQCST = SQRT(1.0-ZSIN**2)                                          
c      ZDECLI=ZCDEC(1)+ZCDEC(2)*ZC1YT+ZCDEC(3)*ZS1YT+ZCDEC(4)*ZC2YT        
c     *       +ZCDEC(5)*ZS2YT                                              
c      ZEQTIM=ZCEQT(1)+ZCEQT(2)*ZC1YT+ZCEQT(3)*ZS1YT+ZCEQT(4)*ZC2YT        
c     *       +ZCEQT(5)*ZS2YT                                              
c      ZZEN1=SIN(ZDECLI)                                                   
c      ZZEN2=COS(ZDECLI)*COS(YCLOCK+ZEQTIM)                                
c      ZZEN3=COS(ZDECLI)*SIN(YCLOCK+ZEQTIM)                                
cC                                                                         
c      ZTIM1 = ZZEN1 * ZSIN                                                
c      ZTIM2 =-ZZEN2 * ZSQCST                                              
c      ZTIM3 = ZZEN3 * ZSQCST                                              
cC                                                                         
cC     ---------------------------------------------------------------     
cC*         2.     COMPUTATIONS IF DIURNAL CYCLE "ON".                     
cC                 ------------ -- ------- ----- -----                     
cC                                                                         
cc      COSLON=1.0                                                         
cc      SINLON=0.0                                                         
c      COSLON=COS(XLON)                                                    
c      SINLON=SIN(XLON)                                                    
c      IF(LDIUR) THEN                                                      
c        AMU0=ZTIM1+ZTIM2*COSLON+ZTIM3*SINLON                              
c        IF(AMU0.GE.0.)THEN                                                
c          AMU0=AMU0                                                       
c          RDAYL=1.0                                                       
c        ELSE                                                              
c          AMU0=0.0                                                        
c          RDAYL=0.0                                                       
c        ENDIF                                                             
cC                 ------------------------------------                    
cC*         3.     COMPUTATIONS IF DIURNAL CYCLE "OFF".                    
cC                 ------------ -- ------- ----- ------                    
c      ELSE                                                                
c        DO 301 JLON=1,128                                                 
c          ZL=2.*API*(JLON-1.)/128.                                        
c          ZMU0(JLON)=ZTIM1+ZTIM2*COS(ZL)+ZTIM3*SIN(ZL)                    
c          IF(ZMU0(JLON).GE.0.)THEN                                        
c            ZMU0(JLON)=ZMU0(JLON)                                         
c            ZRDAYL(JLON)=1.0                                              
c          ELSE                                                            
c            ZMU0(JLON)=0.0                                                
c            ZRDAYL(JLON)=0.0                                              
c          ENDIF                                                           
c  301   CONTINUE                                                          
c        ZS1=SIGMA2(128,ZMU0(1),1)                                         
c        ZS2=SIGMA2(128,ZRDAYL(1),1)                                       
c        IF(ZS2.NE.0.) THEN                                                
c          ZS1=ZS1/ZS2                                                     
c          ZS2=ZS2/128.                                                    
c        END IF                                                            
c        AMU0=ZS1                                                          
c        RDAYL=ZS2                                                         
c      END IF                                                              
cC                                                                         
c      AMU0=CRAE/(SQRT(AMU0**2+ZCRAE)-AMU0)                                
c                                                                          
c      RETURN                                                              
c      END                                                                 
C*************************************                                    
C   FUNCTION SIGMA2                                                       
C*************************************                                    
C                                                                         
       REAL FUNCTION SIGMA2(N, SX, INCX)                                  
C                                                                         
       INTEGER N, INCX, I                                                 
       REAL SX(N*INCX)                                                    
C                                                                         
       SIGMA2 = 0.                                                        
       DO 1 I=1,N                                                         
          SIGMA2 = SIGMA2 + SX(I*INCX)                                    
 1     CONTINUE                                                           
C                                                                         
       END                                                                
c*************************************                                    
      REAL FUNCTION SVP(T)                                                
                                                                          
C Calculates the saturation vapour pressure given T, using Murray's       
C formulae. If t < 253 K assume saturation w.r.t. ice                     
C SVP returned in Pascals                                                 
                                                                          
      IMPLICIT NONE                                                       
                                                                          
      REAL T                                                              
                                                                          
      IF(T.GE.253.0)THEN                                                  
        SVP=611*EXP(17.27*(T-273.16)/(T-35.86))                           
      ELSE                                                                
        SVP=611*EXP(21.87*(T-273.16)/(T-7.66))                            
      ENDIF                                                               
                                                                          
      END                                                                 
                                                                          
