      SUBROUTINE INI_NETCDF                                               
C                                                                         
C     Determines model resolution                                         
C                                                                         
      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
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
C                                                                         
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
c                                                                         
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
c                                                                         
      real xcoord(mg),ycoord(jg*nhem),zcoord(nl),tcoord(10000)            
      character*200 fname1                                                
c                                                                         
      do i=1,mg                                                           
         xcoord(i)=(i-1)*360.0/real(mg)                                   
      end do                                                              
c                                                                         
      do j=1,jg*nhem                                                      
         ycoord(j)=alat(j)                                                
      end do                                                              
c                                                                         
      do l=1,nl                                                           
         zcoord(l)=1000.0*sigma(l)                                        
      end do                                                              
c                                                                         
      itime=(ktotal-kstart)/kountp+1                                      
      print*,' Number of output times = ',itime                           
      do i=1,itime                                                        
         tcoord(i)=i                                                      
      end do                                                              
c                                                                         
      call setup_nc(mg,jg*nhem,nl,itime,                                  
     :              nmaxdims,nall,                                        
     :              ndim,nvar,natts,nattsvar,vdims,                       
     :              vadims,ndims,                                         
     :              dimname,varname,attdimname,                           
     :              attvarname)                                           
c                                                                         
c changing to my directory structure - marc 12/3/03                       
c      fname1='/home/swsvalde/igcm.nc'                                    
      fname1='/home/marc/igcm.nc'                                         
      ifname1=lnsig(fname1)                                               
c                                                                         
      call ininc(fname1(1:ifname1),                                       
     :           nmaxdims,ndim,nvar,                                      
     :           natts,nattsvar,                                          
     :           vdims,vadims,ndims,                                      
     :           dimname,varname,                                         
     :           attdimname,attvarname,                                   
     :           nc1,iddim1,idvar1)                                       
c                                                                         
      print*,' ready to write rest '                                      
c                                                                         
      call writedim(nc1,iddim1(1),xcoord)                                 
      call writedim(nc1,iddim1(2),ycoord)                                 
      call writedim(nc1,iddim1(3),zcoord)                                 
      call writedim(nc1,iddim1(4),tcoord)                                 
c                                                                         
      print*,' Finished writing dimensions '                              
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      SUBROUTINE END_NETCDF                                               
c                                                                         
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
c                                                                         
      print*,' calling netcdf '                                           
      call closenc(nc1)                                                   
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine setup_nc(nlon,nlat,nl,ntime,                             
     :                 nmaxdims,nall,                                     
     :                 ndim,nvar,natts,nattsvar,vdims,                    
     :                 vadims,ndims,                                      
     :                 dimname,varname,attdimname,                        
     :                 attvarname)                                        
      implicit none                                                       
      integer nlon,nlat,nl,ntime,nmaxdims,nall                            
      integer ndim,nvar,natts(nall),nattsvar(nall),                       
     :        vdims(nall),vadims(nmaxdims,nall),                          
     :        ndims(nall)                                                 
      character*200 dimname(nall),varname(nall),                          
     :              attdimname(2,nmaxdims,nall),                          
     :              attvarname(2,nmaxdims,nall)                           
c                                                                         
      integer loc_dim                                                     
c                                                                         
c     This sets up a file similar to .pc files                            
c                                                                         
      ndim=0                                                              
c                                                                         
      ndim=ndim+1                                                         
      dimname(ndim)='longitude'                                           
      ndims(ndim)=nlon                                                    
      natts(ndim)=2                                                       
      attdimname(1,1,ndim)='long_name'                                    
      attdimname(2,1,ndim)='longitude'                                    
      attdimname(1,2,ndim)='units'                                        
      attdimname(2,2,ndim)='degrees_east'                                 
c                                                                         
      ndim=ndim+1                                                         
      dimname(ndim)='latitude'                                            
      ndims(ndim)=nlat                                                    
      natts(ndim)=2                                                       
      attdimname(1,1,ndim)='long_name'                                    
      attdimname(2,1,ndim)='latitude'                                     
      attdimname(1,2,ndim)='units'                                        
      attdimname(2,2,ndim)='degrees_north'                                
c                                                                         
      ndim=ndim+1                                                         
      dimname(ndim)='p'                                                   
      ndims(ndim)=nl                                                      
      natts(ndim)=2                                                       
      attdimname(1,1,ndim)='units'                                        
      attdimname(2,1,ndim)='mb'                                           
      attdimname(1,2,ndim)='positive'                                     
      attdimname(2,2,ndim)='down'                                         
c                                                                         
      ndim=ndim+1                                                         
      dimname(ndim)='time'                                                
      ndims(ndim)=ntime                                                   
      natts(ndim)=2                                                       
      attdimname(1,1,ndim)='long_name'                                    
      attdimname(2,1,ndim)='time'                                         
      attdimname(1,2,ndim)='units'                                        
      attdimname(2,2,ndim)='days or months'                               
c                                                                         
      nvar=1                                                              
      varname(nvar)='zonal_wind'                                          
      vdims(nvar)=4                                                       
      vadims(1,nvar)=loc_dim('longitude',dimname,nall)                    
      vadims(2,nvar)=loc_dim('latitude',dimname,nall)                     
      vadims(3,nvar)=loc_dim('p',dimname,nall)                            
      vadims(4,nvar)=loc_dim('time',dimname,nall)                         
      nattsvar(nvar)=2                                                    
      attvarname(1,1,nvar)='long_name'                                    
      attvarname(2,1,nvar)=                                               
     :   'Zonal Wind'                                                     
      attvarname(1,2,nvar)='units'                                        
      attvarname(2,2,nvar)='ms-1'                                         
c                                                                         
      nvar=nvar+1                                                         
      varname(nvar)='temp'                                                
      vdims(nvar)=4                                                       
      vadims(1,nvar)=loc_dim('longitude',dimname,nall)                    
      vadims(2,nvar)=loc_dim('latitude',dimname,nall)                     
      vadims(3,nvar)=loc_dim('p',dimname,nall)                            
      vadims(4,nvar)=loc_dim('time',dimname,nall)                         
      nattsvar(nvar)=2                                                    
      attvarname(1,1,nvar)='long_name'                                    
      attvarname(2,1,nvar)=                                               
     :   'Temperature'                                                    
      attvarname(1,2,nvar)='units'                                        
      attvarname(2,2,nvar)='K'                                            
c                                                                         
      nvar=nvar+1                                                         
      varname(nvar)='convrain'                                            
      vdims(nvar)=3                                                       
      vadims(1,nvar)=loc_dim('longitude',dimname,nall)                    
      vadims(2,nvar)=loc_dim('latitude',dimname,nall)                     
      vadims(3,nvar)=loc_dim('time',dimname,nall)                         
      nattsvar(nvar)=2                                                    
      attvarname(1,1,nvar)='long_name'                                    
      attvarname(2,1,nvar)=                                               
     :   'Convective Rainfall'                                            
      attvarname(1,2,nvar)='units'                                        
      attvarname(2,2,nvar)='unknown'                                      
c                                                                         
      nvar=nvar+1                                                         
      varname(nvar)='lscalerain'                                          
      vdims(nvar)=3                                                       
      vadims(1,nvar)=loc_dim('longitude',dimname,nall)                    
      vadims(2,nvar)=loc_dim('latitude',dimname,nall)                     
      vadims(3,nvar)=loc_dim('time',dimname,nall)                         
      nattsvar(nvar)=2                                                    
      attvarname(1,1,nvar)='long_name'                                    
      attvarname(2,1,nvar)=                                               
     :   'Large Scale Rainfall'                                           
      attvarname(1,2,nvar)='units'                                        
      attvarname(2,2,nvar)='unknown'                                      
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine netout2                                                  
C                                                                         
C     Determines model resolution                                         
C                                                                         
      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
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
      PARAMETER(RAD=180./PI)                                              
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
C                                                                         
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
c                                                                         
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
      REAL LUG(MG+1,JGG),LTG(MG+1,JGG)                                    
      CHARACTER*14 TEXT1                                                  
c                                                                         
      INTEGER INETCOUNT                                                   
      SAVE INETCOUNT                                                      
      DATA INETCOUNT/0/                                                   
c                                                                         
      INETCOUNT=INETCOUNT+1                                               
c                                                                         
      rewind 24                                                           
      do j=1,jg                                                           
         sec=alat(j)/rad                                                  
         sec=cv/cos(sec)                                                  
                                                                          
         read(24)assbl,ashbl,aslbl,arrcr,arrlr,                           
     :           arflux,asfld,acld,                                       
     :           ssbl,shbl,slbl,rrcr,rrlr,                                
     :           rflux,sfld,cld                                           
         do ihem=1,nhem                                                   
            iof=(ihem-1)*mgpp                                             
            if (ihem.eq.1) then                                           
               jj=j                                                       
            else                                                          
               jj=jggp-j                                                  
            end if                                                        
            do i=1,mg                                                     
               lug(i,jj)=rrcr(i+iof)                                      
            end do                                                        
            lug(mg+1,jj)=lug(1,jj)                                        
            do i=1,mg                                                     
               ltg(i,jj)=rrlr(i+iof)                                      
            end do                                                        
            ltg(mg+1,jj)=ltg(1,jj)                                        
         end do                                                           
      end do                                                              
      rewind 24                                                           
c                                                                         
      call writevar2(nc1,idvar1(3),lug,                                   
     :                     1,mg,1,jgg,inetcount,inetcount,0,0)            
      call writevar2(nc1,idvar1(4),ltg,                                   
     :                     1,mg,1,jgg,inetcount,inetcount,0,0)            
c                                                                         
      return                                                              
      end                                                                 
