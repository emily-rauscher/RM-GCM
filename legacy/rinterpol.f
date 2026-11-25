      subroutine interpol(v,x,u,nvx,nu,result)

!...Linear interpolation

      integer :: nvx, nu
      logical :: if_inc
      real, dimension(nvx) :: v, x
      real, dimension(nu) :: u, result

      if_inc = .false.
      if( x(2) > x(1) ) if_inc = .true.

      if( if_inc ) then

        do i = 1, nu

          jj_x = 0

          do j = 1, nvx-1
            if( u(i) >= x(j) .and. u(i) < x(j+1) ) jj_x = j
          enddo

          if( jj_x > 0 ) then
            result(i) = v(jj_x) + (v(jj_x+1)-v(jj_x)) * (u(i)-x(jj_x)) /
     & (x(jj_x+1)-x(jj_x))
          else if( u(i) < x(1) ) then
            result(i) = v(1) + (v(2)-v(1)) * (u(i)-x(1)) / (x(2)-x(1))
          else
            result(i) = v(nvx) + (v(nvx)-v(nvx-1)) * (u(i)-x(nvx)) /
     & (x(nvx)-x(nvx-1))
          endif

        enddo

      else

        do i = 1, nu

          jj_x = 0

          do j = 1, nvx-1
            if( u(i) <= x(j) .and. u(i) > x(j+1) ) jj_x = j
          enddo

          if( jj_x > 0 ) then
            result(i) = v(jj_x) + (v(jj_x+1)-v(jj_x)) * (u(i)-x(jj_x)) /
     & (x(jj_x+1)-x(jj_x))
          else if( u(i) > x(1) ) then
            result(i) = v(1) + (v(2)-v(1)) * (u(i)-x(1)) / (x(2)-x(1))
          else
            result(i) = v(nvx) + (v(nvx)-v(nvx-1)) * (u(i)-x(nvx)) /
     & (x(nvx)-x(nvx-1))
          endif

        enddo

      endif


      end subroutine interpol

