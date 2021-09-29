      subroutine radout
!
!     **************************************************************
!     *  Purpose             :  Prints out radiative transfer      *
!     *                         results.                           *
!     **************************************************************

      include 'rcommons.h'
      real t(NZ),p(NZ)      
      p=p_aerad
      t=t_aerad
      write(LUNORAD,*) 'radout: u0 = ',u0
      write(LUNORAD,*) 'radout: alb(toa)=',alb_toai
      write(LUNORAD,*) 'radout: alb(tom)=',alb_tomi
      write(LUNORAD,*) 'solar energy absorbed=',fnetbs(nlayer)-fnetbs(1)
      write(LUNORAD,*) 'ir energy absorbed=',fnetbi(nlayer)-fnetbi(1)
!
!       PRINT OUT THE INPUT VARIABLES
!
        WRITE(LUNOPRT,560)
! 560    FORMAT(" RADOUT:",/,      &
!               " j   p(j)    press(j)    t(j) ",      &
!               "     tt(j)  rdh2o(j)  ctot       firu     fird  ",
!               &
!               "    fsLu     fsLd   ")
 560   Format( " p         t        firu     fird      fsLu     fsLd")
!
        DO 565 J = 1, NVERT
           ctot = 0.
           do ig = 1, NGROUP
             do i = 1, NRAD
               ctot = ctot + caer(i,j,ig)
               write(*,*) ctot, caer(i,j,ig)
             enddo
           enddo
           stop


!           WRITE(LUNOPRT,562) J,P(J),PRESS(J),T(J),TT(J),      &
!                        RDH2O(J),Ctot,fupbi(j),fdownbi(j),      &
!                        fupbs(j),fdownbs(j)
           WRITE(LUNOPRT,569) p(j), t(j), fupbi(j),fdownbi(j),      
     &                   fupbs(j),fdownbs(j)
 562       FORMAT(I4,11(1PE9.2))
 569       FORMAT(6(1PE9.2))
 565    CONTINUE
        WRITE(LUNOPRT,563) NLAYER,PRESS(NLAYER),TT(NLAYER)
 563    FORMAT(I3,9X,1PE9.2,9X,1PE9.2)
!
!       PRINT OUT THE CALCULATED VARIABLES
!
        WRITE(LUNOPRT,526)
 526    FORMAT("***delta scaled***:",/      
     &      "  solnet     xirdown     irup       u0",      
     &      "       opd(1)     opd(8)     opd(16)    opd(26)    opd(38)")
        WRITE(LUNOPRT,530) SOLNET,XIRDOWN,XIRUP,      
     &          U0,OPD(1,NLAYER),OPD(8,NLAYER),OPD(16,NLAYER),      
     &          OPD(26,NLAYER),OPD(38,NLAYER)
 530    FORMAT(9(1PE10.3,1X))
!
        WRITE(LUNOPRT,536)
 536    FORMAT("***delta scaled***:"/,      
     &     "              j    heats",      
     &     "       heati        heat         wol         gol        opd",
     &     "        taul")
!
        DO 550 J         =  1,NVERT
           WRITE(LUNOPRT,540) J,HEATS(J),HEATI(J),      
     &                   HEAT(J),WOL(8,J),GOL(8,J),OPD(8,J),TAUL(8,J)
 540       FORMAT(10X,I4,7(1PE11.3,1X))
 550    CONTINUE

      i_more_outpt = 0
      if (i_more_outpt.ne.0) then
        j = NVERT/2. + 1.5
        write(LUNORAD,780) j-1
 780    format(//,'RADOUT: unscaled optical properties for layer',i3,/,
     &   '                ----------- cloud --------------- ',      
     &   ' ----------cloud and gas---------',/      
     &   '   i     wave        tau         w0         g0 ',      
     &   '       tau         w0         g0       ')
        do 790 i = 1, nwave
          index = nprobi(i,1)
          w0cloud = taus(i,j)/taua(i,j)
          write(LUNORAD,781) i,wave(i),taua(i,j),w0cloud,goL(index,j),
     &         utauL(index,j),uw0(index,j),ug0(index,j)
 781      format(1x,i4,7(1pe11.2))
          do 790 ii = 2, nprobi(i,2)
            index = index + 1
            write(LUNORAD,784) utauL(index,j),uw0(index,j),ug0(index,j)
 784        format(49x,3(1pe11.2))
 790    continue
        write(LUNORAD,785) nwave+1,wave(nwave+1)
 785    format(1x,i4,1pe11.2)
      endif

      write(LUNORAD,566)
 566  format(//,' RADOUT: Top of atmosphere radiative fluxes:',/      
     &    '       -------------------solar----------------- ',      
     &    ' --infra-red (top of model)--'/,      
     &    '    i    wave         up       down     albedo   ',     
     &    '       wave     up')

      write(LUNORAD,567) (i,wave(i),fsLu(i),fsLd(i),alb_toa(i),      
     &   .5*(wave(i+nsol)+wave(i+1+nsol)),firu(i),i=1,nir)
 567  format((1x,i4,3(1pe11.2),0pf9.3,2x,2(1pe11.2)))
      write(LUNORAD,568) (i,wave(i),fsLu(i),fsLd(i),alb_toa(i),      
     &             i=nir+1,nsol)
 568  format ((1x,i4,3(1pe11.2),0pf9.3))

      write(LUNORAD,782) tslu,tsld,alb_toai,tiru
 782  format(' totals:        ',2(1pe11.2),0pf9.3,11x,1pe11.2)

      write(LUNORAD,556)
 556  format(//,' RADOUT: total unscaled optical depth',/      
     &    '    i    wave        opd')
      write(LUNORAD,537) (i,wave(i),uopd(i,nlayer),i=1,nwave)
!      print*, 'optical depth = ', uopd(9,nlayer)
 537  format((1x,i4,1p2e10.2,0p))

      return
      end

