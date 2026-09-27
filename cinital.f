C************************************************************             
C                   SUBROUTINE INITAL                                     
C************************************************************             
      SUBROUTINE INITAL                                                   
C                                                                         
C     INITAL calls other initialisation routines.                         
C                                                                         
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'
C      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
C     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
C                                                                         
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
C     Needed so the gas/cloud opacity setup (called below, after
C     INISIMPRAD has read with_TiO_and_VO/opacity_method from fort.7)
C     sees the same values as ciniset.f.
C
      COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM(3), AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS(3), with_TiO_and_VO, opacity_method

      REAL SURFEMIS,ABSSW,ABSLW,ALBSW
      LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     & ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS
      REAL with_TiO_and_VO
      CHARACTER(len=6) :: opacity_method
      COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &   MAXTAU,MAXTAULOC,TCON,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS,GRAYCLDV,
     &   C_TO_O
      CHARACTER(30) :: AEROSOLMODEL
      CHARACTER(30) :: AEROSOLCOMP
      REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),MAXTAU,TCON(NL+1)
      REAL MOLEF(13)
      REAL MTLX, METALLICITY, C_TO_O
      INTEGER AERLAYERS
      LOGICAL DELTASCALE, HAZES, PICKET_FENCE_CLOUDS, GRAYCLDV
      COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN, TAULIMIT
      REAL :: OOM_IN, RFCOEFF_IN, BOTRELAXTIME, FBASEFLUX
      LOGICAL :: LPLOTMAP
      INTEGER :: NLPLOTMAP_IN, NTSTEP_IN, NSKIP_IN

      CALL INIGAU

!!KM Modif
      CALL INIVARPARAM
      CALL INISI
      CALL INIPHYS
      IF (LRESTIJ) THEN
        CALL INIRESIJ                                                     
      ELSE                                                                
        CALL INIRES                                                       
      ENDIF
      CALL INISTR

!! KM Modif
!      CALL INISURF
      CALL INISIMPRAD
C
C     Moved here (from ciniset.f) because with_TiO_and_VO/opacity_method
C     are only set by the INSIMPRAD namelist read just above; calling
C     this from ciniset.f ran before that read, so it always loaded the
C     '_witiovo' table regardless of fort.7 (self-corrected on the
C     second INISET/INITAL pass of a fresh start, but restarts only get
C     one pass and never picked up the correct table).
      WRITE(*,*) 'C to O ratio: ', C_TO_O
      CALL get_cloud_scattering_properties_wrapper
      CALL get_gas_opacity_corrk_wrapper(METALLICITY, C_TO_O, FBASEFLUX,
     &                                   GASCON, with_TiO_and_VO)

      CALL INIQS
      END                                                                 
