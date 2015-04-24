C**********************************************************               
C             SUBROUTINE O3INTERP                                         
C**********************************************************               
      subroutine O3INTERP(o3clim,o3mod,ps)                                
c                                                                         
c This subroutine takes in a profile from climatology (o3clim) and        
c interpolates to the model vertical grid (nl.gt.15), computing the       
c model profile by constraining its column to agree with climatology.     
c The model profile is passed back in o3mod.                              
c N.B. The o3clim passed to the routine must be on the pressure levels    
c defined in the routine (ocliml).                                        
c                                                                         
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
C     Array ordering in SPECTR must correspond to that in GRIDP.          
C                                                                         
      COMMON/SPECTR/Z(IGB),D(IGB),T(IGB),TRA(IGB,NTRAC),SP(IGA),GS(IGA)   
     :              ,SPA(IGA),VP(IGA),DTE(IGB),TT(IGB)                    
     :              ,TRAT(IGB,NTRAC),DT(IGB),ZT(IGB)                      
     :              ,ZMI(IGB),DMI(IGB),TMI(IGB)                           
     :              ,TRAMI(IGB,NTRAC),SPMI(IGA)                           
      COMPLEX Z,D,T,TRA,SP,GS,SPA,VP,DTE,TT,TRAT,DT,ZT                    
     :       ,ZMI,DMI,TMI,TRAMI,SPMI                                      
C                                                                         
c                                                                         
      real o3clim(15),ocliml(15),halflev(14),deltap(15),ohdcol(15)        
      real o3mod(nl),dpmod(nl),modcol(nl)                                 
      real ps, pressu                                                     
      logical ifirst                                                      
      save                                                                
      data ocliml/1.0,3.0,10.0,30.0,50.0,70.0,100.0,150.0,200.0,          
     &              250.0,300.0,400.0,500.0,700.0,850.0/                  
      data ifirst/.true./                                                 
c                                                                         
      rhoo3=2.14 ! density of O3 at STP                                   
      GR=GA*rhoo3 !g*density of O3 at STP                                 
      if (ifirst) then                                                    
        ifirst=.false.                                                    
        do l=1,14                                                         
          halflev(l)=(ocliml(l+1)+ocliml(l))/2.0    ! in mb               
        enddo                                                             
        deltap(1)=halflev(1)                        ! Top                 
        do l=2,14                                                         
          deltap(l)=halflev(l)-halflev(l-1)         ! In between          
        enddo                                                             
        deltap(15)=1.0e3-halflev(14)                ! Bottom              
        do l=1,15                                                         
          deltap(l)=deltap(l)*100.0          ! Convert all to Pa          
        enddo                                                             
      endif                                                               
                                                                          
                               ! Compute overhead column of climatology   
      col=0.0                  ! for this particular profile.             
      do l=1,15                                                           
        col=col+o3clim(l)*deltap(l)                                       
        ohdcol(l)=col*1.0e5/(GR)    ! In DU                               
      enddo                                                               
                                                                          
      do l=1,nl                    ! model delta pressure                 
        dpmod(l)=dsigma(l)*ps      ! in Pa                                
      enddo                                                               
                                                                          
                                                                          
      sumrdp=0.0     ! MMR*dp summed down column                          
                                                                          
      do 269 l=1,nl                                                       
        IF (l.LT.NL)THEN                                                  
           pressu=sigmah(l)*ps/100.0   ! in mb                            
        ELSE                                                              
           pressu=ps/100.0                                                
        ENDIF                                                             
        do 270 lo=1,14                                                    
          if (pressu.lt.halflev(1)) then                                  
                                                                          
            modcol(l)=(pressu/halflev(1))*ohdcol(1)  ! Model column (DU)  
                                               ! forced from climatology  
            o3mod(l)=(modcol(l)*GR/1.0e5-sumrdp)/dpmod(l)                 
            goto 271                                                      
                                                                          
          elseif (pressu.eq.halflev(lo)) then                             
                                                                          
            modcol(l)=ohdcol(lo)                                          
            o3mod(l)=(modcol(l)*GR/1.0e5-sumrdp)/dpmod(l)                 
            goto 271                                                      
                                                                          
          elseif (pressu.gt.halflev(lo).and.pressu.lt.halflev(lo+1))      
     &         then                                                       
                                                                          
            fraction=(pressu-halflev(lo))/(halflev(lo+1)-halflev(lo))     
            modcol(l)=ohdcol(lo)+fraction*(ohdcol(lo+1)-ohdcol(lo))       
            o3mod(l)=(modcol(l)*GR/1.0e5-sumrdp)/dpmod(l)                 
            goto 271                                                      
                                                                          
          elseif (pressu.ge.halflev(14)) then                             
                                                                          
            fraction=(pressu-halflev(14))/(1.0e3-halflev(14))             
            modcol(l)=ohdcol(14)+fraction*(ohdcol(15)-ohdcol(14))         
            o3mod(l)=(modcol(l)*GR/1.0e5-sumrdp)/dpmod(l)                 
            goto 271                                                      
                                                                          
          endif                                                           
                                                                          
 270  continue                                                            
 271  sumrdp=sumrdp+o3mod(l)*dpmod(l)                                     
 269  continue                                                            
      RETURN                                                              
      END                                                                 
