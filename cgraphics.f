      SUBROUTINE INI_GRAPHICS                                             
c                                                                         
      rar=0.45                                                            
c                                                                         
      call disset(23)                                                     
      call winini1                                                        
      call xdevic                                                         
      call xdspac(1.0)                                                    
      call xpspac(0.0,1.9/1.5,0.0,1.0)                                    
      call xspace(2,3,0.0,xlim,ylim)                                      
      call xpscal(1.0,1.0)                                                
      call zjinit                                                         
      call zstret(2.*rar,rar)                                             
      call zstmol(rar)                                                    
      call zstamp(2.*rar)                                                 
      call zstlnt(30.,30.,80.)                                            
      call xpages(-1,0)                                                   
c                                                                         
      call ini_coast                                                      
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine ini_coast                                                
      parameter(idt1=1290,idt2=176000)                                    
      integer inseg                                                       
      real cstpnt                                                         
      common /coasts/ inseg(idt1),cstpnt(idt2,2)                          
c                                                                         
      print*,' Reading coastline data '                                   
      open(unit=10,file='/home/swsvalde/qplot/zjmap.dat')                 
      call zcstrd(10,inseg,idt1,cstpnt,idt2)                              
      close(unit=10)                                                      
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine plotfields(ipage)                                        
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
      COMMON/GRIDSS/ZG1(IGD,JG),DG1(IGD,JG),UG1(IGD,JG),VG1(IGD,JG),      
     :              TG1(IGD,JG),SPG1(IGC,JG)                              
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
      INTEGER LU(JGG,NL),LT(JGG,NL),LVT(JGG,NL)                           
      REAL LUG(MG+1,JGG),LTG(MG+1,JGG)                                    
      REAL LZUG(JGG,NL),LZTG(JGG,NL)                                      
      CHARACTER*14 TEXT1                                                  
c                                                                         
      do j=1,jg                                                           
         sec=alat(j)/rad                                                  
         sec=cv/cos(sec)                                                  
         do ihem=1,nhem                                                   
            iof=(ihem-1)*mgpp                                             
            if (ihem.eq.1) then                                           
               jj=j                                                       
            else                                                          
               jj=jggp-j                                                  
            end if                                                        
            iof1=iof+(nl-1)*igc                                           
            do i=1,mg                                                     
               lug(i,jj)=sec*ug1(i+iof1,j)                                
            end do                                                        
            lug(mg+1,jj)=lug(1,jj)                                        
            do l=1,nl                                                     
               iof1=iof+(l-1)*igc                                         
               asum=0.0                                                   
               do i=1,mg                                                  
                  asum=asum+ug1(i+iof1,j)                                 
               end do                                                     
               lzug(jj,nl+1-l)=sec*asum/real(mg)                          
            end do                                                        
            iof1=iof+(nl-1)*igc                                           
            do i=1,mg                                                     
               ltg(i,jj)=ct*(tg1(i+iof1,j)+t0(nl))-273.16                 
            end do                                                        
            ltg(mg+1,jj)=ltg(1,jj)                                        
            do l=1,nl                                                     
               iof1=iof+(l-1)*igc                                         
               asum=0.0                                                   
               do i=1,mg                                                  
                  asum=asum+tg1(i+iof1,j)                                 
               end do                                                     
               lztg(jj,nl+1-l)=ct*(asum/real(mg)+t0(l))-273.16            
            end do                                                        
         end do                                                           
      end do                                                              
      rewind 24                                                           
c                                                                         
      text1='Day= 999999.99'                                              
      write(text1(6:14),'(f9.2)')day                                      
c                                                                         
      call xnwpic                                                         
      call frame_graphics(1)                                              
      call contour_graphics(ltg,mg+1,jgg,-0.1,1,1,                        
     :     'Temperature at bottom level '//text1)                         
      call xnwpic                                                         
      call frame_graphics(1)                                              
      call contour_graphics(lug,mg+1,jgg,-0.1,0,1,                        
     :     'Zonal Wind at bottom level '//text1)                          
c                                                                         
      call xnwpic                                                         
      call frame_graphics(11)                                             
      call contour_graphics(lztg,jgg,nl,-0.1,1,11,                        
     :     'Zonal Cross_Section of Temperature '//text1)                  
      call xnwpic                                                         
      call frame_graphics(11)                                             
      call contour_graphics(lzug,jgg,nl,-0.1,0,11,                        
     :     'Zonal Cross_Section of Zonal Wind '//text1)                   
c                                                                         
      if (ipage.eq.1) call xpages(-1,0)                                   
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
c                                                                         
      subroutine plotfields2(ipage)                                       
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
      REAL LUG(MG+1,JGG),LTG(MG+1,JGG)                                    
      CHARACTER*14 TEXT1                                                  
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
      text1='Day= 999999.99'                                              
      write(text1(6:14),'(f9.2)')day                                      
c                                                                         
      call xnwpic                                                         
      call frame_graphics(1)                                              
      call contour_graphics(ltg,mg+1,jgg,-0.1,1,1,                        
     :     'Convective Rainfall at bottom level '//text1)                 
      call xnwpic                                                         
      call frame_graphics(1)                                              
      call contour_graphics(lug,mg+1,jgg,-0.1,0,1,                        
     :     'Large scale Rainfall '//text1)                                
c                                                                         
      if (ipage.eq.1) call xpages(-1,0)                                   
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine frame_graphics(iptyp)                                    
      parameter(idt1=1290,idt2=176000)                                    
      integer inseg                                                       
      real cstpnt                                                         
      common /coasts/ inseg(idt1),cstpnt(idt2,2)                          
c                                                                         
      call drawbox(iptyp,0.45,0)                                          
c                                                                         
      if (iptyp.ne.11) then                                               
         call xlcol(3,4,4,-8)                                             
         call zcstpl(inseg,idt1,cstpnt,idt2)                              
      end if                                                              
c                                                                         
      call xclear                                                         
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine drawbox(iptyp,rar,izoom)                                 
c                                                                         
      if (iptyp.eq.0) then                                                
         call zststl(0)                                                   
         call xmap(-0.95,0.95,-1.05,0.45)                                 
         call zpenup(-2.*rar,-1.0*rar)                                    
         call zpendn( 2.*rar,-1.0*rar)                                    
         call zpendn( 2.*rar,     rar)                                    
         call zpendn(-2.*rar,     rar)                                    
         call zpendn(-2.*rar,-1.0*rar)                                    
      else if (iptyp.eq.11) then                                          
         call zststl(0)                                                   
         call xmap(-1.5,1.5,-0.4,1.2)                                     
         call xlcol(8,4,4,0)                                              
         do l=1,6                                                         
            call zpenup(-1.0,0.2*(l-1))                                   
            call zpendn( 1.0,0.2*(l-1))                                   
         end do                                                           
         do j=1,7                                                         
            call zpenup(-1.0+(j-1)/3.0,0.0)                               
            call zpendn(-1.0+(j-1)/3.0,1.0)                               
         end do                                                           
         call xlcol(1,4,4,0)                                              
         call zpenup(-1.0,0.0)                                            
         call zpendn(-1.0,1.0)                                            
         call zpendn( 1.0,1.0)                                            
         call zpendn( 1.0,0.0)                                            
         call zpendn(-1.0,0.0)                                            
      else                                                                
         call xmap(-0.95,0.95,-0.85,0.65)                                 
         call zststl(iabs(iptyp))                                         
         if (iptyp.eq.4) then                                             
            call zstcen(10.0,90.0)                                        
            call zsthta(30.0)                                             
         else if (iptyp.eq.-4) then                                       
            call zstcen(10.0,-90.0)                                       
            if (izoom.eq.0) then                                          
               call zsthta(-30.0)                                         
            else                                                          
               call zsthta(-15.0)                                         
            end if                                                        
         else                                                             
            call zstcen(0.,0.)                                            
         end if                                                           
         call xlcol(8,4,4,0)                                              
         call zxyref                                                      
         call xlcol(1,4,4,0)                                              
         call zoutln                                                      
      end if                                                              
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine contour_graphics(zg,ix,iy,cinc,itype,iptyp,              
     :           text1)                                                   
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
C======================================================================   
C                                                                     =   
C     Do contour plot - line drawing method                           =   
C                                                                     =   
C======================================================================   
      real xcoord(mg+1),ycoord(jgg)                                       
      real cnt(100)                                                       
      integer iwk((mg+1)*jgg)                                             
      character*(*) text1                                                 
      character*43 text2                                                  
c                                                                         
      if (iptyp.eq.1) then                                                
         do i=1,ix                                                        
            xcoord(i)=(i-1)*pi2/real(mg)                                  
         end do                                                           
         do j=1,iy                                                        
            ycoord(j)=pi/180.0*alat(j)                                    
         end do                                                           
      else if (iptyp.eq.11) then                                          
         do i=1,ix                                                        
            xcoord(i)=alat(i)/90.0                                        
         end do                                                           
         do l=1,iy                                                        
            ycoord(l)=1.0-sigma(iy+1-l)                                   
         end do                                                           
      end if                                                              
c                                                                         
      call maxmin(zg,ix,iy,fmin,fmax)                                     
c                                                                         
      ncnt=100                                                            
      call setcv(FMIN,FMAX,CINC,CNT,NCNT,ITYPE,FINC)                      
      text2='min =  999.999, max=  999.999, ci=  999.999'                 
      write(text2(7:14),'(f8.3)')fmin                                     
      write(text2(22:29),'(f8.3)')fmax                                    
      write(text2(36:43),'(f8.3)')finc                                    
c                                                                         
      do ival=1,ncnt                                                      
        if (cnt(ival).gt.0.1*FINC) then                                   
            call xlcol(2,4,4,0)                                           
         else if (cnt(ival).lt.-0.1*FINC) then                            
            call xlcol(4,4,4,0)                                           
         else                                                             
            call xlcol(6,4,4,0)                                           
         end if                                                           
         CALL ZGCNTR (ZG,xcoord,ycoord,IWK,IX,IX,IY,CNT(IVAL))            
      end do                                                              
      call xchmag(0.017)                                                  
      if (iptyp.eq.1) then                                                
         xoff=-0.9                                                        
         yoff=-0.45                                                       
      else                                                                
         xoff=-1.0                                                        
         yoff=0.0                                                         
      end if                                                              
      ilen=lnsig1(text1)                                                  
      call xcharl(xoff,yoff-0.06,text1(1:ilen))                           
      ilen=lnsig1(text2)                                                  
      call xcharl(xoff,yoff-0.12,text2(1:ilen))                           
      call xclear                                                         
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      subroutine maxmin(zg,ix,iy,fmin1,fmax1)                             
      real zg(ix*iy)                                                      
c                                                                         
      fmin1=zg(1)                                                         
      fmax1=zg(1)                                                         
      do i=2,ix*iy                                                        
         if (zg(i).lt.fmin1) fmin1=zg(i)                                  
         if (zg(i).gt.fmax1) fmax1=zg(i)                                  
      end do                                                              
c                                                                         
      return                                                              
      end                                                                 
c                                                                         
      SUBROUTINE SETCV (FMIN,FMAX,CINC,CV,NCNT,ITYPE,FINC)                
C======================================================================   
C                                                                     =   
C  Set up contour levels for a data array with known min and max      =   
C  values. Three options, depending on sign and magnitude of input    =   
C  contour interval, CINC :                                           =   
C     a) CINC.GT.0.       : Use input contour interval CINC.          =   
C     b) CINC.LT.-1.      : Contour interval is (-CINC) in normalised =   
C                           units, defined such that max. abs. value  =   
C                           is 1000 units.                            =   
C     c) -1.LT.CINC.LT.0. : Contour interval is (-CINC) in normalised =   
C                           units, defined such that next power of 10 =   
C                           greater than max. abs. value is 1 unit.   =   
C  In (c) the contour interval is halved if the max. abs. value is    =   
C  less than half the next power of ten, or divided by 5 if less than =   
C  1/5th of it.                                                       =   
C                                                                     =   
C   If itype.eq.0 then (b) and (c) is calculated using maximum value  =   
C   If itype.eq.1 then (b) and (c) uses range                         =   
C======================================================================   
C                                                                         
      REAL CV(NCNT)                                                       
C                                                                         
      FAMAX=AMAX1(ABS(FMAX),ABS(FMIN))                                    
      FRANGE=FMAX-FMIN                                                    
      IF (FRANGE.LE.FAMAX*1.E-10) THEN                                    
         FINC=0.0                                                         
         NCNT=0                                                           
         RETURN                                                           
      END IF                                                              
c                                                                         
      IF (CINC.GT.0.) THEN                                                
         FINC=CINC                                                        
      ELSE IF (CINC.LE.-1.) THEN                                          
         FINC=-FAMAX*CINC/1000.                                           
      ELSE IF (CINC.LT.0.) THEN                                           
         IF (ITYPE.EQ.0) THEN                                             
            IEXP=NINT(ALOG10(FAMAX)*0.99999)                              
            FA10=10.**IEXP                                                
            FRACT=FAMAX/FA10                                              
         ELSE                                                             
            IEXP=NINT(ALOG10(FRANGE)*0.99999)                             
            FA10=10.**IEXP                                                
            FRACT=FRANGE/FA10                                             
         END IF                                                           
         IF (FRACT.GT.0.5) THEN                                           
            FINC=-FA10*CINC                                               
         ELSE IF (FRACT.GT.0.2) THEN                                      
            FINC=-FA10*CINC*0.5                                           
         ELSE                                                             
            FINC=-FA10*CINC*0.2                                           
         END IF                                                           
      END IF                                                              
c      print*,cinc,iexp,fmin,fmax,frange,fract,fa10,finc                  
C                                                                         
      ILO=INT(FMIN/FINC)                                                  
      IF (FMIN.GE.0.) ILO=ILO+1                                           
      CLV=FLOAT(ILO-1)*FINC                                               
      I=0                                                                 
 10   CLV=CLV+FINC                                                        
      I=I+1                                                               
      IF ((I.GT.NCNT).OR.(CLV.GE.FMAX)) GO TO 20                          
      CV(I)=CLV                                                           
      GO TO 10                                                            
 20   NCNT=I-1                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
c                                                                         
      subroutine end_graphics                                             
c                                                                         
      call xgrend                                                         
c                                                                         
      return                                                              
      end                                                                 
