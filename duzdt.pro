pro duzdt,oom,p0,ww,radea=radea,ga=ga,components=components,bfield=bfield,tdrag_min=tdrag_min,$
          psfile=psfile,nmaglevs=nmaglevs

; p0 in bar
; ww in rad/s
; radea in m
; ga in m/s

  if not keyword_set(nmaglevs) then nmaglevs=7
  print,'Using fort.26 and fort.50'

  gascon=3523.                              ; J/kg/K
  if not keyword_set(radea) then radea = 1.e8 ; m
  if not keyword_set(ga) then ga=8.           ; m/s
  print,'           CHECK: gascon, radea, ga:',gascon,radea,ga

; Read File with 3D fields
  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

;  ;get surface presssures
  openr,18,'fort.50'
  readf,18,nlat,nlon
  ab=fltarr(3,nlat,nlon)
  readf,18,ab
  close,18

  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])
  sp=reform(ab[2,*,*])
  sp=(sp+1.)*p0

  sigma=get_sigma(oom,nlev)
  pres=p0*sigma
;  dP=fltarr(nlev)
;  for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1]) ; to match swdp
;  dP[0]=0.5*(pres[1])
;  dP[nlev-1]=0.5*(p0-pres[nlev-2])

  vertvel,ga,gascon,radea,oom,p0,/psurf,/noplot,om=om   ; om output in bar/s

  uzt=fltarr(4,nlat,nlev)       ; four components, 0:merid. mom. by z-mean circ., 1:vert. mom. adv. by z-mean circ.,
                                ;                  2:merid. eddy-mom. convergence, 3:vert. eddy-mom. converg.
  uz=total(u,2)/nlon  ;(nlat,nlev)
  vz=total(v,2)/nlon  ;(nlat,nlev)
  zom=total(om,2)/nlon
  ue=fltarr(nlat,nlon,nlev)
  ve=fltarr(nlat,nlon,nlev)
  ome=fltarr(nlat,nlon,nlev)
  for ilev=0,nlev-1 do begin
     for ilat=0,nlat-1 do begin
        ue[ilat,*,ilev]=u[ilat,*,ilev]-uz[ilat,ilev]
        ve[ilat,*,ilev]=v[ilat,*,ilev]-vz[ilat,ilev]
        ome[ilat,*,ilev]=om[ilat,*,ilev]-zom[ilat,ilev]  ; bar/s
     endfor
  endfor
  ze=total(ue*ve,2)/nlon
  zuom=total(ue*ome,2)/nlon

  toss1=fltarr(nlat,nlev)
  toss2=fltarr(nlat,nlev)
  for ilev=0,nlev-1 do begin
     toss1[*,ilev]=deriv(lat[*,0,0]*!pi/180.,uz[*,ilev]*cos(lat[*,0,0])) ; m/s
     toss2[*,ilev]=deriv(lat[*,0,0]*!pi/180.,ze[*,ilev]*(cos(lat[*,0,0]))^2.) ; m^2/s^2
  endfor
  for ilat=0,nlat-1 do begin
     ; meridional momentum advection by zonal-mean circulation
     uzt[0,ilat,*] = vz[ilat,*] * (2.*ww*sin(lat[ilat,0,0]*!pi/180.) $
                                   - toss1[ilat,*]/radea/cos(lat[ilat,0,0]*!pi/180.)) ; m/s^2
     ; vertical momentum advection by zonal-mean circulation
     uzt[1,ilat,*] = -1.*zom[ilat,*]*deriv(pres,uz[ilat,*])  ; m/s^2
     ; meridional eddy-momentum convergence
     uzt[2,ilat,*] = -1.*toss2[ilat,*]/radea/(cos(lat[ilat,0,0]*!pi/180.))^2 ; m/s^2
     ; vertical eddy-momentum convergence
     uzt[3,ilat,*] = -1.*deriv(pres,zuom[ilat,*])  ; m/s^2
  endfor

  if keyword_set(bfield) then begin
     if not keyword_set(tdrag_min) then tdrag_min=0.
     uztm=fltarr(nlat,nlev)
     for ilat=0,nlat-1 do begin
        ang_f=abs(cos(!pi/2.-lat[ilat,0,0]*!pi/180.))
        for ilev=0,nlev-1 do begin
           rhocgs=reform(sp[ilat,*]*sigma[ilev]/gascon/temp[ilat,*,ilev])*1.e-3 ;cgs
           tdrag=tmag(rhocgs,reform(temp[ilat,*,ilev]),bfield)/ang_f ; in sec
           for ilon=0,nlon-1 do tdrag[ilon]=max([tdrag[ilon],tdrag_min])
           uztm[ilat,ilev]=-1.*total(u[ilat,*,ilev]/tdrag)/nlon
        endfor
     endfor
  endif

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   if keyword_set(components) then begin
      xsz=50
      ysz=30
   endif else begin
      xsz=25
      ysz=25
   endelse
   psopen,filename,/enc,/color,xsize=xsz,ysize=ysz
   !p.font=0
   !x.thick=4
   !y.thick=4
   !p.thick=4
   !p.charsize=2
endif

if not keyword_set(components) then begin

   flat=reform(lat[*,0,0])
   tuzt=total(uzt,1)

   tuzt[0,*]=0.
   tuzt[nlat-1,*]=0.

   umax=max(tuzt,min=umin)
   nlevels=25
   step=(umax-umin)/nlevels
   ulevels=indgen(nlevels)*step+umin

   loadct,10,/silent 
   contour, tuzt,flat,pres, /cell_fill ,LEVELS=ulevels,/ystyle,ymargin=[4,4],$
               yr=[max(pres),min(pres)],/ylog,xtitle='Latitude [degrees]',$
               ytitle='Pressure [bar]',max_value=umax,min_value=umin,/xstyle
   colorbar,position=[0.23,0.965,0.9,0.98],range=[umin,umax],format='(f8.2)',charsize=2
   if keyword_set(bfield) then begin
      uztm[0,*]=0.
      uztm[nlat-1,*]=0.
      umax=max(uztm,min=umin)
      absmax=max([abs(umax),abs(umin)])
      maxrange=2.*absmax        ;from -absmax to +absmax
      nlevels=nmaglevs
;      step=(umax-umin)/nlevels
      step=maxrange/(nlevels-1)
      mylevels=indgen(nlevels)*step-absmax
      contour, uztm,flat,pres, /follow,/over,levels=mylevels,c_charsize=1.5
      contour, uztm,flat,pres,/over,levels=[0.],c_thick=8
      print,'contour levels:',mylevels
      print,'Min., max. mag d(uz)/dt:',minmax(uztm)
   endif else contour,tuzt,flat,pres,/overplot,levels=[0.],color=255,/closed
   loadct,0,/silent

endif else begin

   loadct,10,/silent 

   flat=reform(lat[*,0,0])
   uzt[*,0,*]=0.
   uzt[*,nlat-1,*]=0.

   umax=max(uzt[0,*,*],min=umin)
   nlevels=25
   step=(umax-umin)/nlevels
   ulevels=indgen(nlevels)*step+umin
   contour, reform(uzt[0,*,*]),flat,pres, /cell_fill ,LEVELS=ulevels,/ystyle,position=[0.1,0.6,0.45,0.9],$
               yr=[max(pres),min(pres)],/ylog,xtitle='Latitude [degrees]',$
               ytitle='Pressure [bar]',max_value=umax,min_value=umin,/xstyle
;   contour,reform(uzt[0,*,*]),flat,pres,/overplot,levels=[0.],color=255,/closed
   colorbar,position=[0.125,0.92,0.425,0.95],range=[umin,umax],format='(f8.2)',/top,charsize=1.25

   umax=max(uzt[1,*,*],min=umin)
   nlevels=25
   step=(umax-umin)/nlevels
   ulevels=indgen(nlevels)*step+umin
   contour, reform(uzt[1,*,*]),flat,pres, /cell_fill ,LEVELS=ulevels,/ystyle,position=[0.55,0.6,0.9,0.9],$
               yr=[max(pres),min(pres)],/ylog,xtitle='Latitude [degrees]',$
               ytitle='Pressure [bar]',max_value=umax,min_value=umin,/xstyle,/noerase
;   contour,reform(uzt[1,*,*]),flat,pres,/overplot,levels=[0.],color=255,/closed
   colorbar,position=[0.575,0.92,0.875,0.95],range=[umin,umax],format='(f8.2)',/top,charsize=1.25

   umax=max(uzt[2,*,*],min=umin)
   nlevels=25
   step=(umax-umin)/nlevels
   ulevels=indgen(nlevels)*step+umin
   contour, reform(uzt[2,*,*]),flat,pres, /cell_fill ,LEVELS=ulevels,/ystyle,position=[0.1,0.1,0.45,0.4],$
               yr=[max(pres),min(pres)],/ylog,xtitle='Latitude [degrees]',$
               ytitle='Pressure [bar]',max_value=umax,min_value=umin,/xstyle,/noerase
;   contour,reform(uzt[2,*,*]),flat,pres,/overplot,levels=[0.],color=255,/closed
   colorbar,position=[0.125,0.42,0.425,0.45],range=[umin,umax],format='(f8.2)',/top,charsize=1.25

   umax=max(uzt[3,*,*],min=umin)
   nlevels=25
   step=(umax-umin)/nlevels
   ulevels=indgen(nlevels)*step+umin
   contour, reform(uzt[3,*,*]),flat,pres, /cell_fill ,LEVELS=ulevels,/ystyle,position=[0.55,0.1,0.9,0.4],$
               yr=[max(pres),min(pres)],/ylog,xtitle='Latitude [degrees]',$
               ytitle='Pressure [bar]',max_value=umax,min_value=umin,/xstyle,/noerase
;   contour,reform(uzt[3,*,*]),flat,pres,/overplot,levels=[0.],color=255,/closed
   colorbar,position=[0.575,0.42,0.875,0.45],range=[umin,umax],format='(f8.2)',/top,charsize=1.25

   loadct,0,/silent

endelse

if keyword_set(psfile) then begin
   !p.font=-1
   !x.thick=1
   !y.thick=1
   !p.thick=1
   !p.charsize=1
   psclose
   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
   spawn,'convert plot.png eps3:plot.eps'
endif

end
