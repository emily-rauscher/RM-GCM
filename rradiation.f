C**********************************************************               
C             SUBROUTINE RADIATION                                             
C**********************************************************               
      SUBROUTINE RADIATION(TROPHT,IH)                                             
C                                                                         
C     RADIATION SCHEME DERIVED FROM PREVIOUS CMORC.F AND THE
C     TOON CODES (TOON ET AL 1989). THE SCHEME IS CURRENTLY DOUBLE GRAY
C     AND APPLIES THE TWO-STREAM APPROXIMATION WITH QUADRATURE IN THE
C     VISIBLE AND HEMISPHERIC MEAN IN THE INFRARED. 
                          
C     It passes the pressure of the full sigma levels and the surface                                   
C     to the Radiation scheme temperatures from TG and TSTAR                                  
C                                                                         
C     Determines model resolution                                         
      use omp_lib
      include 'params.i'


C     Sets basic constants, especially those needed for array dimensions

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

C     Basic planetary parameters for run plus information about
C     vertical grid structure

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

C     Legendre polynomials and information about gaussian latitudes
C
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC
C
C
C     Array ordering in GRIDP must correspond to that in SPECTR.
C     Real arrays: multi-level arrays are 2-dimensional.
C     the variables have been renamed to coincide with
C     variable names in bgcm5 DGRMLT
C
C swapped round VLNG and UNLG to check
      COMMON/GRIDP/ CHIG(IGC,NL),SFG(IGC,NL),UG(IGC,NL),VG(IGC,NL)
     :              ,TTVD(IGC,NL),QTVD(IGC,NL),TG(IGC,NL)
     :              ,TRAG(IGC,NL,NTRAC)
     :              ,PLG(IGC),TYBL(IGC),TXBL(IGC)
     :              ,SPG(IGC),VPG(IGC),TTRD(IGC,NL)
     :              ,TNLG(IGC,NL),TRANLG(IGC,NL,NTRAC),UNLG(IGC,NL)
     :              ,VNLG(IGC,NL),TTLR(IGC,NL),UTRAG(IGC,NL,NTRAC)
     :              ,TTCR(IGC,NL),VTRAG(IGC,NL,NTRAC)
     :              ,UTVD(IGC,NL),VTVD(IGC,NL)
     :         ,ASSBL(IGC),ASHBL(IGC),ASLBL(IGC),ARRCR(IGC),ARRLR(IGC)
     :         ,arflux(igc,6),asfld(igc,6),acld(igc,4)
     :         ,SSBL(IGC),SHBL(IGC),SLBL(IGC),RRCR(IGC),RRLR(IGC)
     :         ,rflux(igc,6),sfld(igc,6),cld(igc,4)
C
C
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN

       LOGICAL LPLOTMAP

C     Constant arrays and variables associated with time and vertical
C     differencing. Also counters.

      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)
     +              ,BEGDOY,DOY
C
            COMMON/PHYS/  CCR,RCON,DTBUOY,TSLA,TSLB,TSLC,TSLD,CUT1,CUT2
     :              ,TSTAR(IGC,JG),QSTAR(IGC,JG),FRAD(JG,NHEM)
     :              ,TSTARO(IGC,JG),TDEEPO(IGC,JG),smstar(igc,jg)
     :              ,tdeep(igc,jg),hsnow(igc,jg),sqstar(igc,jg)
     :              ,SALB(IGC,JG),SBAL(IGC,JG),BLCD(IGC)
     :              ,SVEGE(IGC,JG),CD,DRAG,BLVAD,BLA,BLRH,BLVB(IGC)
     :              ,AKVV,AKTV,AKQV,ESCONA,ESCONB,EPSIQ,CTQ,CCC
     : ,ctqi,sdsn,shcs,shcsp,shcsn,skse,sksn,slhf,sd1,sd2,sdw
     :        ,ssmc,sdsnd,sasnow,saice,shsstar,shsmax
     :     ,LOC,LNOICE,LOLDBL,LCOND,LNNSK
     :              ,NLCR,CURHM,AKTC,AKQC,CUBMT,CBADJT,CBADJP
     :              ,SKAP(NL),SK(NLM),FWS(NL),CLR(NL),FB(NLM)
     :              ,TTDC(NL),QTDC(NL),TTMC(NL),QTMC(NL),TC(NL),QC(NL)
     :              ,CTCR(NL,NHEM),CTLR(NL,NHEM)
     :              ,LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ
     :              ,LSL,NAVRD,NAVWT,DELT2C,SHCO,SHCI,ITSLL,ITSLO,NCUTOP
                LOGICAL LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ,LSL,LOC
     :       ,LNOICE,LOLDBL,LCOND,LNNSK

       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM, AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS,TEMPERATURE_INTERNAL,TEMPERATURE_IRRAD

       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS

       CHARACTER(30) :: AEROSOLMODEL
       CHARACTER(30) :: AEROSOLCOMP
       REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),TCON(NL+1)
       LOGICAL DELTASCALE
       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &   MAXTAU,MAXTAULOC,TCON,AERSOLCOMP,MTLX,MOLEF,AERLAYERS

      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC
     +              ,RNTAPO
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG
     +              ,LSPO(NL),LGPO(NL)
     $              ,LSHIST,LMINIH
      LOGICAL LSHIST,LMINIH
      LOGICAL LSPO,LGPO



C
C     Setup moisture variables by equivilencing them to
C     Tracer No. 1
C
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))
C
      COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)
     :     ,SNET(IGC,JG),RRFLUX(IGC,JG,6)
     :     ,TTSW(IGC,NL),TTLW(IGC,NL)
      COMMON/GSG/GSG(IGC,JG)
      REAL htnet
      COMMON /RADHT/ HTNET(NHEM,JG,MG,NL)
      REAL TAVE(IGP)

       REAL PR(NL+1),T(NL+1),PRFLUX(nl+1),htlw(nl+1),htsw(nl+1)
       real PRB2T(NL+1),adum


      integer ifirst                ! If =1, first time reading o3
                                    ! and h2o (2 months' worth).

      real amfrac                   ! fraction through month
      integer ichange               ! =1 when in process of month change
*---------------------------
      integer ifirstcol             ! =1 first time through column
                                    ! calculation (open new file).
      real p0
      real ps                       ! sfc pressure (used in
                                    ! interpolation from climatology
                                    ! to model).
      integer im                    ! Pointer for array plg (for
                                    ! getting sfc pressure).

C     Array to hold fluxes at top and bottom of atmosphere
C     1st index - flux 1=SW, 2=LW
C     2nd index - Direction 1=DN, 2=UP
C     3rd index - Where 1=TOP, 2=SURFACE
      real fluxes(2,2,2)

c     The following for parallel testing --MTR
      integer TID, NTHREADS
      double precision test_wctime

      save                          ! Want to keep things like dcompl.

      DATA IFIRST/1/
      data ifirstcol/1/

      RHSCL=288.0*GASCON/GA
      CHRF=86400.*WW*CT   ! factor to non-dimensionalise heating rates

c  Skipping part set up here.
c  Radiation scheme only called every nskip longitudes
c  nskip must divide exactly into mg for longitude.
c  ntstep is the number of timesteps to skip.

      ntstep=NTSTEP_IN
      nskip=NSKIP_IN

      MALSKY_INDEX_NUMBER = MALSKY_INDEX_NUMBER + 1

c --------------------------------- Now start the radiation bit.
c
C loop over hemispheres
      IOFM=0

      DO 800 ihem=1,nhem
        IF (mod(kount,ntstep).eq.0) then
          ilast=0
          DO i=1,mg
            im=i+iofm
            idocalc=0
            IF ((i.eq.1).or.(i-ilast.ge.nskip)) then
              idocalc=1
            ELSE
              IF (LNNSK) THEN
                imp=im+1
                IF (imp.gt.(mg+iofm)) imp=1+iofm
                  imm=im-1
                  IF (((gsg(im,jh).gt.0.).and.(gsg(imp,jh).eq.0.)).or.
     $               ((gsg(im,jh).eq.0.).and.(gsg(imp,jh).gt.0.)).or.
     $               ((gsg(im,jh).gt.0.).and.(gsg(imm,jh).eq.0.)).or.
     $               ((gsg(im,jh).eq.0.).and.(gsg(imm,jh).gt.0.))) THEN
                    idocalc=1
                  ENDIF
                ENDIF
              ENDIF
              IF (idocalc.eq.1) then

c --------------------------------------- Now set rest of column.
                DO LD=1,NL    ! Start of loop over column.
                  L=NL-LD+2  ! Reverse index (Morc goes bottom up).
                  PR(LD)=SIGMA(LD)*PLG(im)*P0 ! Pressure
                  PRB2T(L)=PR(LD)
                  T(LD)=TG(im,ld)*CT ! Temperature
                  AEROPROF(LD)=0.0
                ENDDO

                AEROPROF(NL+1)=0.0
                PRB2T(1)=PLG(im)*P0
                PR(NL+1)=PLG(im)*P0
                T(NL+1)=((FBASEFLUX+rrflux(IM,JH,1))/5.6704e-8)**0.25

! This could also just be the ground temperature... decision to be made
c ----------------------------------------------------- And alat1

                alat1=alat(JH)*REAL(-(ihem*2.)+3)
                IF ((LFLUXDIAG).AND.(KOUNTP-KOUTP.LT.NTSTEP_IN)) THEN
                  IF(JH.EQ.1.AND.IHEM.EQ.1.AND.I.EQ.1) THEN
                    REWIND(63) !! Rewind file for fluxes in nikosrad
                    REWIND(62) ! rwnd file for ancillary RT results

                    IF (PORB.NE.0) THEN
                      SSLON=(1./PORB-1.)*KOUNT*360./ITSPD
                      SSLON=MOD(SSLON,360.)
                    ELSE
                        SSLON=0.
                    ENDIF

                    SSLAT=ASIN(SIN(OBLIQ*PI/180.)*SIN(PI2*KOUNT/ITSPD/PORB))*180./PI
                    WRITE(63,2021) DAY,SSLON,SSLAT
                    WRITE(62,2021) DAY,SSLON,SSLAT
 2021               FORMAT('DAY:',F7.2,', SUBSTELLAR LON,LAT:',2F7.2)
                    WRITE(63,*)
                    WRITE(62,*)''

                  ENDIF
                ENDIF

!! BOTTOM SHORT WAVE ALBEDO SET IN FORT.7 INISIMPRAD,
! Note to future modelers-this is not location or wavelenght dependent
             SWALB=ALBSW

! PR AND T ARE THE TEMPERATURE AND PRESSURE AT THE SIGMA LEVELS
! AND BOTTOM BOUNDARY, AS USED BY THE DYNAMICAL CODE.
! TO COMPUTE HEATING RATES AT THESE CENTERS, WE NEED TO DEFINE
! LAYER EDGES AT WHICH FLUXES ARE COMPUTED, PRFLUX.
             DO LD    = 1,NL-1
                 PRFLUX(LD+1)=(pr(LD)+pr(LD+1))/2.
             ENDDO
             PRFLUX(NL+1)=PR(NL+1)
             PRFLUX(1)=pr(1)*0.5
!             write(*,*)'fluxes',fluxes
C cloud cf and ic passed. fluxes returned.
C which is net flux at TOA in profile
C Call radiation scheme
                  alon=REAL(i-1)/REAL(mg)*360.0

!     PR in pascals for layer boundaries (NL+1), T in Kelvin for layer
!     centers + one layer for the bottom boundary. The top is n=1, the
!     bottom is n=NL+1


! MTR: THE FOLLOWING LOOP IS ONLY FOR THE OLDER VERSION OF THE CODE
! WHERE DISTRIBUTIONS ARE HARDWIRED (ROMAN & RAUSCHER 2017)

            IF((AEROSOLS).AND.(AEROSOLMODEL.NE.'Global')) THEN
!            IF(AEROSOLS.AND.) THEN
!           AEROSOLS TIME!
!           Extract a single column from the array AER4LAT(NLEV,LON,HEM)
!            AEROPROF=AERO4LAT(:,mg,ihem)
              DO  LD=1,NL +1
              AEROPROF(LD)=TAUAEROSOL(LD,i,ihem,ih)
             ENDDO

            ENDIF
 !           call calc_radheat(pr,t,prflux,alat1,alon,htlw,htsw,DOY,cf,
 !    $                 ic,fluxes,swalb,kount,itspd)
            call calc_radheat(pr,t,prflux,alat1,alon,htlw,htsw,DOY,cf,ic,fluxes,swalb,kount,itspd)
           pr=prb2t

                  PNET(IM,JH)=fluxes(1,1,1)-fluxes(1,2,1)+fluxes(2,1,1)-fluxes(2,2,1)
                  SNET(IM,JH)=fluxes(1,1,2)-fluxes(1,2,2)+fluxes(2,1,2)-fluxes(2,2,2)


                  rrflux(im,jh,1)=fluxes(1,1,2)
                  rrflux(im,jh,2)=fluxes(1,2,2)
                  rrflux(im,jh,3)=fluxes(2,1,2)
                  rrflux(im,jh,4)=fluxes(2,2,2)
                  rrflux(im,jh,5)=fluxes(1,1,1)-fluxes(1,2,1)
                  rrflux(im,jh,6)=fluxes(2,2,1)

                  DO l=nl,1,-1
c  bottom heating rate is zero in morecret
                     LD=NL+1-L
                     IM=I+IOFM
                     HTNETO=HTNET(IHem,JH,I,LD)

                     htnet(ihem,jh,i,ld)=(htlw(l+1)+htsw(l+1))


                     TTRD(IM,LD)=(HTNETO
     $                    +HTNET(IHEM,JH,I,LD))/(CHRF*2.)


                     IF ((i-ilast.gt.1).and.(nskip.gt.0)) then
                        write(*,*),i,last
                        DO j=ilast+1,i-1
                           a=REAL(j-ilast)/REAL(i-ilast)
                           b=1.-a

                 write(*,*),'CANNOT SKIP LONGITUDES IN PARALLEL!! ABORT'
                 write(*,*),'Please set nskip=0 in fort.7'
                 STOP
                           HTNETO=HTNET(IHEM,JH,J,LD)
                           htnet(ihem,jh,j,ld)=a*htnet(ihem,jh,i,ld)+
     $                          b*htnet(ihem,jh,ilast,ld)
                           im=j+iofm
                           TTRD(IM,LD)=(HTNETO
     $                          +HTNET(IHEM,JH,J,LD))/(CHRF*2.)
!                    write(*,*)'TTRD',LD,TTRD(IM,LD)
                            IF (l.eq.nl) then
                              pnet(im,jh)=a*pnet(i+iofm,jh)+
     $                             b*pnet(ilast+iofm,jh)
                              snet(im,jh)=a*snet(i+iofm,jh)+
     $                             b*snet(ilast+iofm,jh)

                              DO k=1,6
                                 rrflux(im,jh,k)=a*rrflux(i+iofm,jh,k)
     $                                +b*rrflux(ilast+iofm,jh,k)
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
!                  stop
                  ilast=i
C end of conditional execution of morcrette code
               ENDIF

            ENDDO

            IF (nskip.ne.0) then
               write(*,*),'CANNOT SKIP LONGITUDES IN PARALLEL!! ABORT'
               write(*,*),'Please set nskip=0 in fort.7'
               STOP
            ELSE
                ilast=mg
            ENDIF
            IF (ilast.ne.mg) then

                DO j=ilast+1,mg
                  a=REAL(j-ilast)/REAL(mg+1-ilast)
                  b=1.-a
                  im=j+iofm
                  DO l=nl,1,-1
                     ld=nl+1-l
                     HTNETO=HTNET(IHEM,JH,J,LD)

                     htnet(ihem,jh,j,ld)=a*htnet(ihem,jh,1,ld)+
     $                    b*htnet(ihem,jh,ilast,ld)
                     TTRD(IM,LD)=(HTNET(IHEM,JH,J,LD)
     $                    +HTNETO)/(CHRF*2.)
                   write(*,*) 'Now to the print statement'
                     IF (l.eq.nl) then
                        pnet(im,jh)=a*pnet(1+iofm,jh)+
     $                       b*pnet(ilast+iofm,jh)
                        snet(im,jh)=a*snet(1+iofm,jh)+
     $                       b*snet(ilast+iofm,jh)
                        DO k=1,6
                           rrflux(im,jh,k)=a*rrflux(1+iofm,jh,k)
     $                          +b*rrflux(ilast+iofm,jh,k)
                        ENDDO
                     ENDIF
                  ENDDO
!           stop
               ENDDO
            ENDIF

         ELSE                   ! ntstep requirement
            DO i=1,mg
               DO LD=1,NL
                  im=i+IOFM
                  TTRD(im,LD)=(htnet(ihem,jh,i,ld))/CHRF
               ENDDO
            ENDDO
         ENDIF
         IOFM=MGPP
 800  CONTINUE                  ! end of loop over hemispheres
      IF (LSHORT.AND.(KOUNT.eq.1)) then
         DO l=1,nl
            DO i=1,igc
               ttrd(i,l)=ttrd(i,l)*2.
            ENDDO
         ENDDO
      ENDIF


      !IF (mod(MALSKY_INDEX_NUMBER,50000) .eq. 0) THEN
      !  DO LD=1,NL    ! Start of loop over column.
      !      write(*,*) T(LD), ','
      !  ENDDO
      !  write(*,*)
      !END IF

      RETURN
      END
