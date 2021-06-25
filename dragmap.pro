pro dragmap,level,dtime,psfile=psfile

; plots map of KE dissipated at level=level
; reads fort.26
;
; level=0 is top, level=nlev-1 is bottom
; dtime in seconds
;
; if n_elements(level,dtime) gt 1 then plot,KE diss,pressure
;   and print total over atmosphere

p0=1.e2  ;bar
oom=5
radea=1.e8
ga=8.  ;m/s^2

if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0,retain=2
  !p.font=-1
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20;,bits_per_pixel=24
  !p.font=0
endelse

  !x.style=1
  !p.charsize=1.5

latoff=0
lonoff=0

; Read File with 3D fields

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

if n_elements(level) eq 1 then begin
; Pick vertical level (lower layer is nlev-1)
   nl=level                     ; from function command line
    
   lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl])
   u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl])
   temp=reform(xy[5,*,*,nl])

;mass elements (followed dsigma calc in cinisi.f)
   dm=fltarr(nlat,nlon)
;use relative areas (so still will get W/m^2)
   if level eq nlev-1 then dA=p0*(1.-10.^(-1.*oom/nlev))*(2.*!pi/nlon)*(!pi/nlat)/ga/4./!pi else begin
      if level ne 0 then begin
         dA=p0*(10.^(-1.*oom*(nlev-level-1)/nlev)-10.^(-1.*oom*(nlev-level)/nlev)) $
            *(2.*!pi/nlon)*(!pi/nlat)/ga/4./!pi
      endif else dA=p0*(10.^(-1.*oom*(nlev-1)/nlev))*(2.*!pi/nlon)*(!pi/nlat)/ga/4./!pi
   endelse
   for ilat=0,nlat-1 do begin
      dm[ilat,*]=replicate(cos(lat[ilat,0,0]*!pi/180.)*dA,nlon)
   endfor
   dm*=1.e5                     ;for units conversion to kg

   ke=0.5*dm*(u^2+v^2)          ; J/m^2

   kerate=ke/dtime              ; W/m^2
   print,minmax(kerate)

   flon=lon*360./lon(0,nlon-1)
   flat=lat*90./lat(0,nlat-1)
  
   MAP_SET, /Miller_CYLINDRICAL,0.,0.,/ISOTROPIC ;, TITLE=' Temperature and Velocity Map'
   loadct, 4,/silent
   contour, kerate,flon,flat, /overplot,/cell_fill,/closed ,NLEVELS=45 ;, min_value=1000, max_value=2000
   colorbar,position=[0.1,0.07,0.90,0.1],range=[min(kerate),max(kerate)],format='(e8.1)',charsize=2
   map_grid, /label ,color=255

   Print, 'Total KE loss rate:',total(kerate),' W/m^2'

endif else begin

   if n_elements(level) ne n_elements(dtime) then begin
      print,'n_elements(level) must equal n_elements(dtime)'
      abort
   endif
   lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*])
   u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*])
   temp=reform(xy[5,*,*,*])
   kerate=fltarr(nlev)
   pres=p0*get_sigma(oom,nlev)
   for i=0,n_elements(level)-1 do begin
      dm=fltarr(nlat,nlon)
      if level[i] eq nlev-1 then dA=p0*(1.-10.^(-1.*oom/nlev))*(2.*!pi/nlon)*(!pi/nlat)/ga else begin
         if level[i] ne 0 then begin
            dA=p0*(10.^(-1.*oom*(nlev-level[i]-1)/nlev)-10.^(-1.*oom*(nlev-level[i])/nlev)) $
               *(2.*!pi/nlon)*(!pi/nlat)/ga
         endif else dA=p0*(10.^(-1.*oom*(nlev-1)/nlev))*(2.*!pi/nlon)*(!pi/nlat)/ga
      endelse
      for ilat=0,nlat-1 do dm[ilat,*]=replicate(cos(lat[ilat,0,0]*!pi/180.)*dA,nlon)*radea^2
      dm*=1.e5                  ;for units conversion to kg
      ke=total(0.5*dm*(u[*,*,level[i]]^2+v[*,*,level[i]]^2))  ;J
      kerate[i]=ke/dtime[i]           ; W
   endfor
   plot,kerate,pres,yr=[p0,min(pres)],/ylog,xtit='KE dissipation [W]',ytit='Pressure [bar]'
   Print, 'Total KE loss rate:',total(kerate),' W'
endelse

loadct,0,/silent
  !x.style=0
  !p.charsize=1.

if keyword_set(psfile) then psclose

end
