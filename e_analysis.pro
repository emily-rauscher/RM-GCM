pro e_analysis,nl
  
; calculate conserved energy: D(0.5*v^2+c_p T +gz)/Dt = 0, taking into account
; heating, to isolate local distribution of KE dissipation

  R=3523.
  akap=0.286
  cp=R/akap
  g=8.
  tgr=2900.
  p0=1.e2
  oom=5
  radea=1.e8

  print,'-----WARNING: Double-check parameters and make sure using correct fort.26, fort.50, and fort.63'

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*])
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])

  altitude,R,g,TGR,p0=p0,oom=oom,z=z
  vertvel,g,R,radea,oom,p0,/psurf,ww=ww,/noplot
  eset=fltarr(nlat,nlon,nlev)

; (include vertical velocity?)
  eset=0.5*(u^2+v^2+ww^2) + cp*temp + g*z  ; J/kg

; calculate advective flux of eset:
  feset=fltarr(nlat,nlon,nlev)  ; J/kg/s
  for ilev=0,nlev-1 do begin
     for ilat=0,nlat-1 do begin
        coslat=cos(lat[ilat,0,0]*!pi/180.)
        feset[ilat,*,ilev]+=reform(u[ilat,*,ilev])* $  ; u*
                            deriv(!pi/180.*reform(lon[ilat,*,ilev])*radea*coslat,reform(eset[ilat,*,ilev])) ; d(eset)/dx
     endfor
     for ilon=0,nlon-1 do begin
        feset[*,ilon,ilev]+=reform(v[*,ilon,ilev])* $  ;v*
                            deriv(!pi/180.*reform(lat[*,ilon,ilev])*radea,reform(eset[*,ilon,ilev])) ; d(eset)/dy
     endfor
  endfor
  for ilat=0,nlat-1 do begin
     for ilon=0,nlon-1 do begin
        feset[ilat,ilon,*]+=reform(ww[ilat,ilon,*])* $ ;w*
                            deriv(reform(z[ilat,ilon,*]),reform(eset[ilat,ilon,*])) ; d(eset)/dz
     endfor
  endfor

  qrad=fltarr(nlat,nlon,nlev)
  entry = create_struct( 'lat', 0.0, $
                         'lon', 0.0, $
                         'p0',0.0, $
                         'fpres', fltarr(nlev), $
                         'sflux', fltarr(nlev), $ ; SW flux
                         'flux', fltarr(nlev), $  ; actual IR flux
                         'rflux', fltarr(nlev), $ ; cirrad flux
                         'pres', fltarr(nlev), $
                         'sw', fltarr(nlev), $ ;sw heating
                         'lw', fltarr(nlev))   ;lw heating
  nent=(nlat-1)*nlon
  rates=replicate(entry,nent)
  read1=strarr(2)
  data1=fltarr(5,nlev+1)
  read4=''
  data2=fltarr(3,nlev+1)
  read5=''
  openr,18,'fort.63'
  for i=0,nent-1 do begin
     readf,18,read1
     rates[i].lat=float(strmid(read1[0],23,6))
     rates[i].lon=float(strmid(read1[0],47,6))
     readf,18,data1
     rates[i].fpres=reform(data1[0,0:nlev-1])
     rates[i].rflux=reform(data1[1,0:nlev-1]-data1[2,0:nlev-1])
     rates[i].flux=reform(data1[3,0:nlev-1])
     rates[i].sflux=reform(data1[4,0:nlev-1])
     readf,18,read4
     readf,18,data2
     readf,18,read5
     rates[i].pres=reform(data2[0,0:nlev-1]) ; in bar (or mbar, if keyword set)
     rates[i].p0=data2[0,nlev]
     rates[i].sw=reform(data2[1,0:nlev-1])/24./3600.*cp    ; in J/s/kg
     rates[i].lw=reform(data2[2,0:nlev-1])/24./3600.*cp    ; in J/s/kg
;    fix to deal with rates too small the exponent gets incorrectly dropped (making them big again):
     test=where(data2[1,*] lt min(abs(data2[2,0:nlev-1])))
     if test[0] ne -1 then rates[i].sw[test[0]:nlev-1]=replicate(0.,nlev-test[0])
     ;find nearest lat, lon
     dlat=lat[*,0,0]-rates[i].lat
     toss=min(abs(dlat),lati)
     dlon=lon[0,*,0]-rates[i].lon
     toss=min(abs(dlon),loni)
;     print,lat[lati,0,0],rates[i].lat,'   ',lon[0,loni,0],rates[i].lon
     qrad[lati,loni,*]=rates[i].sw+rates[i].lw  ; J/s/kg
  endfor
  close,18        

; fix for missing latitude
zq=where(qrad[0,*,*] ne 0)
if zq[0] eq -1 then begin
   eset=eset[1:nlat-1,*,*]
   feset=feset[1:nlat-1,*,*]
   qrad=qrad[1:nlat-1,*,*]
   nnlat=nlat-1
   blat=1
endif else begin
   nnlat=nlat
   blat=0
endelse

diss=feset-qrad
ndiss=diss/eset 
print,'Total relative dissipation: ',total(ndiss)


device,true_color=24,decomposed=0,retain=2

if nl ge 0 then begin
   qmin=min(ndiss[*,*,nl],max=qmax)
   qtp=fltarr(nnlat,nlon+1)
   qtp[*,0:nlon-1]=ndiss[*,*,nl]
   qtp[*,nlon]=ndiss[*,0,nl]
   
   flon=fltarr(nnlat,nlon+1)
   flon[*,0:nlon-1]=lon[blat:nlat-1,*,nl]
   flon[*,nlon]=reform(lon[blat:nlat-1,0,nl]+360.)
   flat=fltarr(nnlat,nlon+1)
   flat[*,0:nlon-1]=lat[blat:nlat-1,*,nl]
   flat[*,nlon]=reform(lat[blat:nlat-1,0,nl])

   latoff=0
   lonoff=0

   MAP_SET, /Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC ;, TITLE=' Temperature and Velocity Map'
   loadct, 4,/silent
   cbottom=55.
   nlevels=45.
   step=(qmax-qmin)/nlevels
   mylevels=indgen(nlevels)*step+qmin
   ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

   contour, qtp,flon,flat, /overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels
   contour,qtp,flon,flat,/overplot,levels=[0.]

   if qmin ne qmax then colorbar,position=[0.1,0.07,0.90,0.1],range=[qmin,qmax],format='(e8.1)',$
                                 charsize=2,bottom=cbottom,ncolors=255-cbottom
   loadct,0,/silent
endif else begin
   loadct, 4,/silent
   zave=fltarr(nnlat,nlev)
   pres=p0*get_sigma(oom,nlev)
   flev=fltarr(nnlat,nlev)
   flat=fltarr(nnlat,nlev)
   for ilat=0,nnlat-1 do begin
      flev[ilat,*]=pres
      flat[ilat,*]=reform(lat[ilat+blat,0,*])
      for ilev=0,nlev-1 do begin
         zave[ilat,ilev]=mean(ndiss[ilat,*,ilev])/nlon
      endfor
   endfor
   umax=max(zave,min=umin)
   nlevels=25
   step=(umax-umin)/nlevels
   ulevels=indgen(nlevels)*step+umin
   contour, zave,flat,flev, /cell_fill ,LEVELS=ulevels,/ystyle,ymargin=[4,4],$
            yr=[max(flev),min(flev)],/ylog,xtitle='Latitude [degrees]',$
            ytitle='Pressure [bar]',max_value=umax,min_value=umin,ytickname=ylabels,/xstyle
   colorbar,position=[0.23,0.965,0.9,0.98],range=[umin,umax],format='(e8.1)',charsize=2
   loadct, 0,/silent
   contour,zave,flat,flev,/overplot,levels=[0.]
endelse

;stop

end
