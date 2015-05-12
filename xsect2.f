C**********************************************************               
C             SUBROUTINE XSECT2                                            
C**********************************************************               
      SUBROUTINE XSECT2                                             
C                                                                         
C
C        
C     This subroutine for personal outputs
C     K. Menou                                                                
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
     +              ,MF,MFP,JZF,NF,NFP                                    
     +              ,AKAP,GA,GASCON,RADEA,WW,P0,PFAC,EZ,AIOCT             
     +              ,RD,RV,CPD,CLATNT                                     
     +              ,LRSTRT,LSHORT,LTVEC,LSTRETCH                         
     +              ,LFLUX                                                
     +              ,LBALAN,LRESTIJ                                       
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
     +              ,LNOISE                                               
      COMPLEX EZ,AIOCT                                                    
      LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ                 
     +       ,LFLUX,LNOISE                                                
     +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                           
C                                                                         
C                                                                         
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BM1(IDE),AK(NNP),AQ(NL2),G(NL2),TAU(NL2)              
     +              ,KOUNT,KITS,KSTART,KTOTAL,KRUN,BEGDAY,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,CTRA(NTRAC),KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C                                                                         
C                                                                         
C     Array ordering in GRIDP must correspond to that in SPECTR.          
C     Real arrays: multi-level arrays are 1-dimensional.                  
C                                                                         
      COMMON/GRIDP/ CHIG(IGD),SFG(IGD),UG(IGD),VG(IGD)                    
     *              ,ZG(IGD),DG(IGD),TG(IGD)                              
     +              ,TRAG(IGD,NTRAC)                                      
     *              ,PLG(IGC),PJG(IGC),PMG(IGC)                           
     *              ,SPG(IGC),VPG(IGC),EG(IGD)                            
     +              ,TNLG(IGD),TRANLG(IGD,NTRAC),FUG(IGD),FVG(IGD)        
     +              ,UTG(IGD),UTRAG(IGD,NTRAC)                            
     +              ,VTG(IGD),VTRAG(IGD,NTRAC),FVGT(IGD),FUGT(IGD)        
     $              ,GRPAD(NGRPAD)                                        
C                                                                         
C                                                                         
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
C                       

C                                                                         
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 
      
       LOGICAL LPLOTMAP
                                                  
C     Switches counters and constants controlling type and frequency of   
C     model output                                                        
C                                                                         
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
      COMMON/GRIDSS/ZG1(IGD,JG),DG1(IGD,JG),UG1(IGD,JG),VG1(IGD,JG),      
     :              TG1(IGD,JG),SPG1(IGC,JG)                              
C                                                                         
C                                                                         
C     Restoration fields and timescale                                    
C                                                                         
      COMMON/RESTOR/ZRES(IGN),DRES(IGN),TRES(IGN),SPRES(IGM),DAMP         
C                                                                         
      INTEGER LU(IGG),LT(IGG),LVT(IGG)                                    
c                                                                         
C KM Modif needed to pass rrflux array (OLR in fort.29)--ER Modif to fort.64
      COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     &               ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     &               ,TTSW(IGC,NL),TTLW(IGC,NL)                           

C
      integer nmaxdims,nall                                               
      parameter(nmaxdims=4,nall=10)                                       
      integer ndim,nvar,natts(nall),nattsvar(nall),                       
     :        vdims(nall),vadims(nmaxdims,nall),                          
     :        ndims(nall)                                                 
      character*200 dimname(nall),varname(nall),                          
     :              attdimname(2,nmaxdims,nall),                          
     :              attvarname(2,nmaxdims,nall)                           
      integer nc1,iddim1(nall),idvar1(nall)                               
      common/netcdfvars/nc1,iddim1,idvar1                                 
c                                                                         
      REAL LUG(MG,JGG,NL),LTG(MG,JGG,NL)                                  
      LOGICAL LFIRST                                                      
      INTEGER INETCOUNT                                                   
      SAVE LFIRST,INETCOUNT                                               
      DATA LFIRST/ .TRUE. /                                               
      DATA INETCOUNT/0/                                                   


C     Array to store fields until loops are done
      REAL*8  TOUTST(NL,JG*2,MG,3)
      REAL*8  TOUTLAT,TOUTLON
C------ Emily modif 1/3
      REAL*8  TOUTU(NL,3)
      REAL*8  TOUTT(NL,6)
C     TOUTU(L,1) is max u, TOUTU(L,2) is min, TOUTU(L,3) is zonal average of u for L
C     TOUTT(L,*) are T at: sub,anti-stellar, E,W limb (at eqtr), N,S pole for L
      REAL*8  TOUTE(NL,2)
      REAL*8  TOUTSP(JG*2,MG)
C----- KM modif for OLR flux i fort.29      
      REAL*8 TOUTFLUX(MG,JG*2)
C------
      REAL SSLON,SSLAT
C--- erin for binary flux output                                                                                                  
      REAL BINFLUX

C PGPLOT map variables
C      REAL*4 TRMAT(6)
C      REAL*4 BRIGHT
C      REAL*4 CONTRA
C      REAL*4 PLOTMAP(MG,JG*2),UPLOT(MG,JG*2),VPLOT(MG,JG*2)
C      REAL*4 VMIN, VMAX
C      INTEGER NLPLOTMAP
C      INTEGER C1,C2,NC

C      VMIN=1e30
C      VMAX=-1e30

C------ Emily modif 2/3 
C     Initialize TOUTU, TOUTT
      DO 50 L=1,NL
         DO 51 K=1,3
            TOUTU(L,K)=0.
 51      CONTINUE
         DO 52 K=1,6
            TOUTT(L,K)=0.
 52      CONTINUE
         DO 53 K=1,2
            TOUTE(L,K)=0.
 53      CONTINUE
 50   CONTINUE


C------

C                                                                         
C     Output is wanted. Read grid point fields from stream 24             
C                                                                         
      REWIND 26 
      REWIND 50
      REWIND 64

           
      WRITE (26,106) JG*2,MG,NL
      WRITE (50,101) JG*2,MG
      WRITE (64,106) JG*2,MG,NL
 106  FORMAT(3I5)     
 101  FORMAT(2I5)

                                              
      RMG=1./REAL(MG)                                                     
!!      WRITE (2,104) DAY,DOY                                               
!!  104 FORMAT(/' CROSS SECTIONS FOR DAY, DOY',2F6.1)                       
      DO 10 J=1,JG                                                        
C                                                                         
C        Loop over latitudes to read gridpoint data, dimensionalise       
C        and calculate zonal means                                        
C                                                                         
         SEC=ALAT(J)/57.29578                                             
         DMASSJ=RADEA*RADEA*COS(SEC)*(PI2/MG)*(PI/2./JG)/GA
C         SEC=10.*CV/COS(SEC)                                              
         SEC=CV/COS(SEC)                                              
C         READ(24) ZG,DG,UG,VG,TG,SPG                                     
         DO I=1,IGD                                                       
            ZG(I)=ZG1(I,J)                                                
            DG(I)=DG1(I,J)                                                
            UG(I)=UG1(I,J)                                                
            VG(I)=VG1(I,J)                                                
            TG(I)=TG1(I,J)                                                
         END DO                                                           
         DO I=1,IGC                                                       
            SPG(I)=SPG1(I,J)                                              
         END DO                                                           

C     Loop over longitudes for surface pressure
C
         DO 13 IHEM=1,NHEM
            IP=(IHEM-1)*MGPP
            DO 31 I=1,MG
               IF (IHEM.EQ.1) THEN
                  TOUTSP(J,I)=SPG(IP+I)
                  IF (NHEM.EQ.1) THEN
                     TOUTSP(JG*2-(J-1),I)=SPG(IP+I)
                  ENDIF
               ELSE
                  TOUTSP(JG*NHEM-(J-1),I)=SPG(IP+I)
               ENDIF
 31            CONTINUE
 13         CONTINUE

C                                                                         
C        Loop for [U] and [T] sections                                    
C                                                                         
         DO L=1,NL                                                        
            DMASS=DMASSJ*DSIGMA(L)
         DO IHEM=1,NHEM                                                   
               IF (IHEM.EQ.1) THEN                                        
                  JJ=J                                                    
               ELSE                                                       
                  JJ=JGGP-J                                               
               ENDIF                                                      
            IOF=MGPP*(IHEM-1)+(L-1)*IGC                                   
            DO I=1,MG                                                     
!                  TOUTST(L,JG*NHEM-(JJ-1),I,1)=0.1*SEC*UG(I+IOF)
!                  TOUTST(L,JG*NHEM-(JJ-1),I,2)=0.1*SEC*VG(I+IOF)
!                  TOUTST(L,JG*NHEM-(JJ-1),I,3)=CT*(TG(I+IOF)+T0(L)) 

                  UB=UG(IOF+I)*SEC
                  VB=VG(IOF+I)*SEC
                  TB=(TG(IOF+I)+T0(L))*CT

                  EK=DMASS*(TOUTSP(J,I)+1.)*0.5*(UB*UB+VB*VB)

                  TOUTST(L,JJ,I,1)=UB
                  TOUTST(L,JJ,I,2)=VB
                  TOUTST(L,JJ,I,3)=TB

                  TOUTE(L,1)=TOUTE(L,1)+EK
                  TOUTE(L,2)=TOUTE(L,2)+EK
     &                 +DMASS*(TOUTSP(J,I)+1.)*(GASCON/AKAP)*TB



!!               LUG(I,JJ,L)=0.1*SEC*UG(I+IOF)                              
!!            END DO                                                        
!!            DO I=1,MG                                                     
!!               LTG(I,JJ,L)=CT*(TG(I+IOF)+T0(L))                           
            END DO                                                        
        end do                                                            
        end do                                                            

   10 CONTINUE                                                            
C                                                                         
CC      call writevar2(nc1,idvar1(1),lug,                                 
CC     :                     1,mg,1,jgg,1,nl,inetcount,inetcount)         
CC      call writevar2(nc1,idvar1(2),ltg,                                 
CC     :                     1,mg,1,jgg,1,nl,inetcount,inetcount)         


      DO 20 L=1,NL   
       DO 21 I=1,MG 
          TOUTLON= 360.0*(I-1)/MG
        DO 22 J=1,JG*2
           TOUTLAT=ALAT(J)
                  WRITE (26,107) TOUTLON,TOUTLAT,L,
     &  TOUTST(L,J,I,1),TOUTST(L,J,I,2),TOUTST(L,J,I,3)
              IF (L.EQ.NL) WRITE (50,102) TOUTLON,TOUTLAT,TOUTSP(J,I)
   22            CONTINUE         
   21    CONTINUE                 
   20 CONTINUE                    

 107              FORMAT(2E13.5,I4,2E13.5)
 102              FORMAT(3E13.5)

C-------- Emily Modif 3/3
      DO 40 L=1,NL
         DO 41 I=1,MG
            TOUTU(L,1)=MAX(TOUTU(L,1),TOUTST(L,JG,I,1),
     &                     TOUTST(L,JG+1,I,1))
            TOUTU(L,2)=MIN(TOUTU(L,2),TOUTST(L,JG,I,1),
     &                     TOUTST(L,JG+1,I,1))
 41      CONTINUE
C         TOUTU(L,1)=MAX(TOUTST(L,JG,1:MG,1),TOUTST(L,JG+1,1:MG,1))
C         TOUTU(L,2)=MIN(TOUTST(L,JG,1:MG,1),TOUTST(L,JG+1,1:MG,1))
         TOUTU(L,3)=SUM(TOUTST(L,JG:JG+1,1:MG,1))/MG/2.
         TOUTT(L,1)=SUM(TOUTST(L,JG:JG+1,1,3))/2.
         TOUTT(L,2)=SUM(TOUTST(L,JG:JG+1,MG/2+1,3))/2.
         TOUTT(L,3)=SUM(TOUTST(L,JG:JG+1,MG/4+1,3))/2.
         TOUTT(L,4)=SUM(TOUTST(L,JG:JG+1,MG*3/4+1,3))/2.
         TOUTT(L,5)=SUM(TOUTST(L,1,1:MG,3))/MG
         TOUTT(L,6)=SUM(TOUTST(L,JG*2,1:MG,3))/MG
         WRITE (27,109) TOUTU(L,1),TOUTU(L,2),TOUTU(L,3)
 109     FORMAT(3E13.5)
         WRITE (28,110) TOUTT(L,1),TOUTT(L,2),TOUTT(L,3),TOUTT(L,4),
     &        TOUTT(L,5),TOUTT(L,6)
 110     FORMAT(6E13.5)
         WRITE (52,111) TOUTE(L,1),TOUTE(L,2)
 111     FORMAT(2E13.5)
 40   CONTINUE

C--------


C---------- KM Modif for OLR output in fort.29 -- ER modif: fort.64
      IOFM=0

       DO IHEM=1,NHEM
       DO I=1,MG 
C          TOUTLON= 360.0*(I-1)/MG
          IM=I+IOFM
        DO J=1,JG
C           JJ=J + JG*(IHEM-1)
C           TOUTLAT=ALAT(JJ)

C                  WRITE (29,111) TOUTLON,TOUTLAT, RRFLUX(IM,J,6) 
C 111              FORMAT(3E13.5)

           IF (IHEM.EQ.1) THEN
           TOUTFLUX(I,J)=RRFLUX(IM,J,6)
           ELSE
           TOUTFLUX(I,2*JG-J+1)=RRFLUX(IM,J,6)
           ENDIF
           
        ENDDO         
       ENDDO
       IOFM=MGPP
       ENDDO

       DO I=1,MG 
          TOUTLON= 360.0*(I-1)/MG
        DO J=1,JG*2
           TOUTLAT=ALAT(J)
                  WRITE (64,109) TOUTLON,TOUTLAT, TOUTFLUX(I,J)
        ENDDO         
       ENDDO


C Define coordinate range of graph and draw axes     
C      CALL PGSCI(1)
C      CALL PGENV(1., float(n), 200., 340.,  0,  0)
C      CALL PGENV(1., float(n), -8., 8.,  0,  0)
C Label the axes (note use of \u and \d for raising exponent).
C      CALL PGLAB('Latitude grid', 'Temperature', 'Latitudinal Profiles')

!       CALL PGWNAD(0.0, 1.0+nx, 0.0, 1.0+ny)palett
!       CALL pgimag(f)
     
C       TRMAT(1) = 0.0
C       TRMAT(2) = 1.0/(MG)
C       TRMAT(3) = 0.0
C       TRMAT(4) = 0.0
C       TRMAT(5) = 0.0
C       TRMAT(6) = 1.0/(JG*2)

!        NLPLOTMAP=NL
!        NLPLOTMAP=1
!        NLPLOTMAP=(NL+1)/2  ! Plot mid-level
C       NLPLOTMAP=NLPLOTMAP_IN

C      CALL PGQCIR(C1, C2)
C      NC = MAX(0, C2-C1+1)
C      WRITE (*,*) 'Number of color indices used for image: ', NC
C      IF (NC .LT.8) THEN 
C         WRITE (*,*) 'Not enough colors available on this device'
C         STOP
C      ELSE
C         WRITE (*,*)
C      END IF

!       write(*,*) LPLOTMAP

C       IF(LPLOTMAP) THEN

C      CALL PGPAGE
!      CALL SETVP
!      CALL PGWNAD(0.0, (1.0+MG), 0.0, (1.0+JG*2))
C Set up the color map.
C      BRIGHT = 0.5      
C      CONTRA  = 1.0
C      CALL PALETT(2, CONTRA, BRIGHT)

C       DO L=1,NL   
C       DO I=1,MG 
C          TOUTLON= 360.0*(I-1)/MG
C        DO J=1,JG*2
C           TOUTLAT=ALAT(J)

C           IF (L.EQ.NLPLOTMAP) THEN
C              PLOTMAP(I,J)=REAL(TOUTST(L,J,I,3))
C              UPLOT(I,J)=REAL(TOUTST(L,J,I,1))
C              VPLOT(I,J)=REAL(TOUTST(L,J,I,2))
!              PLOTMAP(J,I)=REAL(cos(TOUTLON)*cos(TOUTLAT))
!              write(*,*) 'PLOTMAP', J, I, PLOTMAP(J,I)
C              VMIN=MIN(PLOTMAP(I,J),VMIN)
C              VMAX=MAX(PLOTMAP(I,J),VMAX)
C           ENDIF
C         ENDDO
C        ENDDO
C       ENDDO

C      write(*,*) VMIN, VMAX
    

C Draw the map with PGIMAG.
C      CALL PGIMAG(PLOTMAP,MG,JG*2,1,MG,1,JG*2,VMIN,VMAX,TRMAT)
!!      CALL PGSCH(2.0)

!      CALL PGSAH(1,45.0,0.3)
!      CALL PGVECT(UPLOT, VPLOT, MG, JG*2, 1, MG, 1, JG*2, 1.0, 0, 
!     & TRMAT,0.0)

!!       CALL PGSCH(1.0)
      

C  Annotate the plot
C      CALL PGMTXT('t',1.0,0.0,0.0,
C     & 'TEMPERATURE MAP AND VELOCITY VECTORS')

C      ENDIF

      IF (PORB.NE.0) THEN 
         SSLON=(1./PORB-1.)*KOUNT*360./ITSPD  ! in degrees
         SSLON=MOD(SSLON,360.)         
      ELSE
         SSLON=0.
      ENDIF
      SSLAT=OBLIQ*SIN(PI2*KOUNT/ITSPD/PORB)

      WRITE (26,105) DAY,SSLON,SSLAT
      WRITE (50,105) DAY,SSLON,SSLAT
      WRITE (64,105) DAY,SSLON,SSLAT
 105  FORMAT(/' OUTPUTS FOR DAY ',F10.4,', SUBSTELLAR LON, LAT:',2F8.3)

      IF (mod(KOUNT,ITSPD)=0) THEN
 201  FORMAT(E13.5)
        call BinaryFLux(KOUNT,BINFLUX)
        write(89,201) BINFLUX
      ENDIF 
      
      RETURN                                                              
      END                                                                 
                                                                          

C&&&&&&&&&&&&&&& ADDED PGPLOT ROUTINE &&&&&&&&&&&&&&

c      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used byC PGIMAG.
C-----------------------------------------------------------------------
c      INTEGER TYPE
c      REAL*4 CONTRA, BRIGHT
C
c      REAL*4 GL(2), GR(2), GG(2), GB(2)
c      REAL*4 RL(9), RR(9), RG(9), RB(9)
c      REAL*4 HL(5), HR(5), HG(5), HB(5) 
c      REAL*4 WL(10), WR(10), WG(10), WB(10)
c      REAL*4 AL(20), AR(20), AG(20), AB(20)
C
c      DATA GL /0.0, 1.0/
c      DATA GR /0.0, 1.0/
c      DATA GG /0.0, 1.0/
c      DATA GB /0.0, 1.0/

c      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
c      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
c      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/ 
c      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/

c      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
c      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
c      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
c      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/

c      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
c      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
c      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
c      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
c      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
c     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
c      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
c     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
c      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
c     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
c      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
c     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

c      IF (TYPE.EQ.1) THEN
C        -- gray scale
c         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
c      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
c         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
c      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
c         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
c      ELSE IF (TYPE.EQ.4) THEN
C        -- weird IRAF
c         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
c      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
c         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
c      END IF
c      END

c      SUBROUTINE SETVP
C-----------------------------------------------------------------------
C Set the viewport, allowing margins around the edge for annotation.C (This is similar in effect to PGVSTD, but has different margins.)C The routine determines the view-surface size and allocates marginsC as fractions of the minimum of width and height.
C-----------------------------------------------------------------------
c      REAL*4 D, VPX1, VPX2, VPY1, VPY2
c      CALL PGSVP(0.0, 1.0, 0.0, 1.0)
c      CALL PGQVP(1, VPX1, VPX2, VPY1, VPY2)
c      D = MIN(VPX2-VPX1, VPY2-VPY1)/40.0
c      VPX1 = VPX1 + 5.0*D
c      VPX2 = VPX2 - 2.0*D
c      VPY1 = VPY1 + 8.0*D
c      VPY2 = VPY2 - 2.0*D
c      CALL PGVSIZ(VPX1, VPX2, VPY1, VPY2)
c      END
