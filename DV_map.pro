pro DV_map,levplot,latoff,lonoff,qtp=qtp,psfile=psfile,moreinfo=moreinfo,vfrac=vfrac,prange=prange

; plots divergence field (or relative vorticity, for qtp ne 0)
; reads from fort.51
; if moreinfo set, will also read fort.26 and overplot contours for T
;    and vectors for wind

if not keyword_set(psfile) then begin
  device,true_color=24,decomposed=0,retain=2
  !p.font=-1
  !p.charsize=1.5
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20;,bits_per_pixel=24
  !p.font=0
  !p.charsize=2
endelse

  !x.style=1


openr,17,'fort.51'
readf,17,nlat,nlon,nlev
if levplot ge nlev then print,'ERROR: levplot > nlev'
xy=fltarr(5,nlat,nlon,nlev)
readf,17,xy
close,17

lon=reform(xy[0,*,*,levplot]) & lat=reform(xy[1,*,*,levplot])
div=reform(xy[3,*,*,levplot]) & vort=reform(xy[4,*,*,levplot])

if keyword_set(moreinfo) then begin
   openr,18,'fort.26'
   readf,18,nlat,nlon,nlev
   ab=fltarr(6,nlat,nlon,nlev)
   readf,18,ab
   close,18
   u=reform(ab[3,*,*,levplot]) & v=reform(ab[4,*,*,levplot]) & temp=reform(ab[5,*,*,levplot])
endif

if not keyword_set(qtp) then begin
   qtp=div 
   ptitle='Divergence field'
endif else begin
   qtp=vort
   ptitle='Relative vorticity'
endelse

flon=lon*360./lon(0,nlon-1)
flat=lat*90./lat(0,nlat-1)
if not keyword_set(latoff) then latoff=0.
if not keyword_set(lonoff) then lonoff=0.

print,minmax(div),minmax(vort)

loadct,4
if not keyword_set(psfile) then MAP_SET,/Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC, TITLE=ptitle $
  else MAP_SET, /Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC

if not keyword_set(prange) then prange=[min(qtp),max(qtp)]

if prange[0] ne prange[1] then $
   colorbar,position=[0.1,0.08,0.90,0.12],range=prange,format='(E8.1)',charsize=1.5

contour,qtp,flon,flat,/overplot,/cell_fill,nlevels=55,min_value=prange[0],max_value=prange[1]

contour,qtp,flon,flat,/overplot,levels=[0.],color=0

if keyword_set(moreinfo) then begin
   ;temperature contours
   contour,temp,flon,flat,/overplot,nlevels=5,color=255,thick=6
   ;wind vectors
   if not keyword_set(vfrac) then vfrac=0.9
   partvelvec,u,v,flon,flat,/over,fraction=vfrac,color=0
endif

map_grid,/label,color=255,charsize=1.5

print,'L',levplot+1

if keyword_set(psfile) then psclose

end
