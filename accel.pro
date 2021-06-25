pro accel,levplot,psfile=psfile

if not keyword_set(psfile) then begin
  device,true_color=24,decomposed=0,retain=2
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20;,bits_per_pixel=24
  !p.font=0
  !x.thick=6
  !y.thick=6
  !p.thick=6
  !p.charsize=2
endelse

  ga=15.
  gascon=3523.
  radea=9.e7
  oom=6
  p0=1.e3

  print,'Check variables and need both fort.26 and fort.50'

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

   lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
   u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])

   pres=p0*get_sigma(oom,nlev)

   vertvel,ga,gascon,radea,oom,p0,/psurf,/noplot,om=omega

   udev=fltarr(nlat,nlon,nlev)
   vdev=fltarr(nlat,nlon,nlev)
   odev=fltarr(nlat,nlon,nlev)
   zlat=(!pi/180.)*reform(lat[*,0,0])
   hacc=fltarr(nlat,nlev)
   vacc=fltarr(nlat,nlev)
   vcor=fltarr(nlat,nlev)
   for ilev=0,nlev-1 do begin
      toss=fltarr(nlat)
      for ilat=0,nlat-1 do begin
         udev[ilat,*,ilev]=u[ilat,*,ilev]-mean(u[ilat,*,ilev])
         vdev[ilat,*,ilev]=v[ilat,*,ilev]-mean(v[ilat,*,ilev])
         odev[ilat,*,ilev]=omega[ilat,*,ilev]-mean(omega[ilat,*,ilev])
         toss[ilat]=(cos(zlat[ilat]))^2 * mean(udev[ilat,*,ilev]*vdev[ilat,*,ilev])
         vcor[ilat,ilev]=mean(udev[ilat,*,ilev]*odev[ilat,*,ilev])
      endfor
      hacc[*,ilev]=-1./radea/(cos(zlat))^2 * deriv(zlat,toss)
   endfor
   for ilat=0,nlat-1 do vacc[ilat,*]=-1. * deriv(pres,reform(vcor[ilat,*]))

   if levplot ge 0 then begin
      hp=reform(hacc[1:nlat-2,levplot])
      vp=reform(vacc[1:nlat-2,levplot])
      zp=zlat[1:nlat-2]*180./!pi
      amin=min([hp,vp],max=amax)
      plot,hp,zp,yr=[-90,90],xr=[amin,amax],xtit='Acceleration [m/s^2]',ytit='Latitude [degrees]',linestyle=2,/ystyle
      oplot,vp,zp,linestyle=1
      oplot,hp+vp,zp
      oplot,[0,0],[-90,90]
      al_legend,['Horizontal','Vertical','Total'],linestyle=[2,1,0],charsize=1.5
   endif else begin
      hp=reform(hacc[nlat/2,*])
      vp=reform(vacc[nlat/2,*])
      amin=min([hp,vp],max=amax)
      plot,hp,pres,/ylog,yr=[p0,min(pres)],xtit='Acceleration [m/s^2]',ytit='Pressure [bar]',$
           linestyle=2,/ystyle
      oplot,vp,pres,linestyle=1
      oplot,hp+vp,pres
      oplot,[0,0],[p0,min(pres)]
      al_legend,['Horizontal','Vertical','Total'],linestyle=[2,1,0],charsize=1.5,/bottom
   endelse

   if keyword_set(psfile) then begin
      !p.font=-1
      !x.thick=1
      !y.thick=1
      !p.thick=1
      !p.charsize=1
      psclose
   endif

end
