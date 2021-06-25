PRO albedoproj,inparm,latoff,lonoff,nl,grid=GRID,psfile=psfile,spherical_alb=spherical_alb,mu0corr=mu0corr,toaswfluxes=toafluxes,lons=lons,lats=lats,prange=prange,ctcol=ctcol

;Program reads in fort.26 and fort.62 to map out the top of the atmosphere
;albedoes SW, upwelling LW flux, and pressure where tau=input tau.

lonset=0.
;latoff=0.0
;lonoff=0.0
filler=''
location=strarr(4)
nlat=1.
nlon=1.
nlev=1.
lat=1.
lon=1.
sslat=1.
sslon=1.
data=1.0
data2=[1.0,1.0]
rotat=0.0 

;restore the clear case
;restore,f='/u/Exocomputer1/mtroman/RM-GCM/Kepler7b/T31/case40clear/clear30irtau.sav'
;qtpclear=qtp



if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0,retain=2
;  set_plot, 'x'
  window, 0, xsize=x_size_window, ysize=y_size_window, retain=2
  !p.font=-1
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=23.5,ysize=20,bcolor=1;,bits_per_pixel=24
  !p.font=0
   print,'psfile set'
endelse

     OpenR, 17, 'fort.26'
     ReadF, 17, nlat,nlon,nlev
     Close, 17
     
     nlev      = nlev+1.
     toaa      = fltarr(nlon,nlat)
     toma      = fltarr(nlon,nlat)
     IRflux    = fltarr(nlon,nlat) 
     mu0       = fltarr(nlon,nlat)
     outgoing  = fltarr(nlon,nlat)
     outgoingIR= fltarr(nlon,nlat)
     incoming  = fltarr(nlon,nlat)
     toafluxes = fltarr(nlon,nlat,2)
     lats      = fltarr(nlat)
     lons      = fltarr(nlon)
     fullarr   = fltarr(12,nlev)
     pbars     = fltarr(nlon,nlat,nlev)
     tauintsw  = fltarr(nlon,nlat,nlev)
     tauintlw  = fltarr(nlon,nlat,nlev)
     tmisw     = fltarr(nlon,nlat,nlev)
     tmilw     = fltarr(nlon,nlat,nlev)
     pi0s      = fltarr(nlon,nlat,nlev)
     directsw  = fltarr(nlon,nlat,nlev)
     i=1
     j=1


     OpenR, 17, 'fort.62'
     ReadF, 17, filler,sslat,sslon,format=['(A34,F5.2,F5.2)']
;Loop over lats, loop over lons
    FOR i  =  0,nlat-1 DO BEGIN
      FOR j  =  0,nlon-1 DO BEGIN
        readF, 17, filler
        readF, 17, filler,lat,filler,lon,format=['(A21,F9.5,A12,F9.5)']
        IF (lat lt 0.) THEN ii = nlat-(i+1.)/2. ELSE ii = i/2.
;        IF (lon eq 0.) then stop
        lats(ii)=lat
        lons(j)=lon
        readF, 17, filler,data, format=['(A39,E23.22)']
        mu0(j,ii) = data
        readF, 17, filler,data, format=['(A30,E23.22)']
        toaa(j,ii) = data
        readF, 17, filler,data, format=['(A25,E23.22)']
        toma(j,ii) = data
        readF, 17, filler
        readF, 17, filler,data, format=['(A45,E30.10)']
        IRflux(j,ii)  = data
        readF, 17, filler
        readF, 17, data2
        toafluxes(j,ii,*)=data2
        readF, 17, filler
        readF, 17, fullarr
        pbars(j,ii,*)       =    fullarr(0,*) 
        tauintsw(j,ii,*)    =    fullarr(3,*)
        tauintlw(j,ii,*)    =    fullarr(4,*)
        pi0s(j,ii,*)        =    fullarr(5,*) 
        tmisw(j,ii,*)       =    fullarr(10,*)
        tmilw(j,ii,*)       =    fullarr(11,*)
        directsw(j,ii,*)    =    fullarr(9,*)
      ENDFOR
    ENDFOR
     CLOSE,17	
     ;Now shift things, nominally so that the sub-stellar point is centered
     sf=lonset/360.*nlon 
     lons  =  shift(lons,sf)
     mu0   =  shift(mu0,sf,0)
     toaa  =  shift(toaa,sf,0)
     toma  =  shift(toma,sf,0)
     IRflux=  shift(IRflux,sf,0)
     pbars =  shift(pbars,sf,0,0)
     tauintsw  = shift(tauintsw,sf,0,0)
     tauintlw  = shift(tauintlw,sf,0,0)
     tmisw     = shift(tmisw,sf,0,0)
     tmilw     = shift(tmilw,sf,0,0)
     pi0s      = shift(pi0s,sf,0,0) 
     directsw  = shift(directsw,sf,0,0)
     tauintsw  = tauintsw(*,*,nl)
     tauintlw  = tauintlw(*,*,nl)
     tmisw     = tmisw(*,*,nl)
     tmilw     = tmilw(*,*,nl)
     pi0s      = pi0s(*,*,nl)
     directsw  = directsw(*,*,nl)


     mu0big  = rebin(mu0,nlon*10.,nlat*10.,/sample)
     toaabig  = rebin(toaa,nlon*10.,nlat*10.,/sample)
     tomabig  = rebin(toma,nlon*10.,nlat*10.,/sample)
     IRfluxbig  = rebin(IRflux,nlon*10.,nlat*10.,/sample)
     pi0big     = rebin(pi0s,nlon*10.,nlat*10.,/sample)
     directswbig   = rebin(directsw,nlon*10.,nlat*10.,/sample)
     tauintswbig = rebin(tauintsw,nlon*10.,nlat*10.,/sample)
     tauintlwbig = rebin(tauintlw,nlon*10.,nlat*10.,/sample)
;     tmiswbig = rebin(tmisw(*,*,nl),nlon*10.,nlat*10.)
 ;    tmilwbig = rebin(tmilw(*,*,nl),nlon*10.,nlat*10.)     
   param = mu0
   if inparm eq 1 then param  =  mu0
   if inparm eq 2 then param  =  toaa
   if inparm eq 3 then param  =  toma
   if inparm eq 4 then param  =  IRflux
   if inparm eq 5 then param  =  pi0s
   if inparm eq 6 then param  =  tmisw
   if inparm eq 7 then param  =  tmilw
   if inparm eq 8 then param  = directsw
   if inparm eq 9 then param  = tauintsw
   if inparm eq 10 then param = tauintlw


   if inparm eq 1 then print, 'mu0'
   if inparm eq 2 then print, 'toaa'
   if inparm eq 3 then print, 'toma'
   if inparm eq 4 then print, 'IRflux'
   if inparm eq 5 then print, 'pi0s'
   if inparm eq 6 then print, 'tmisw'
   if inparm eq 7 then print, 'tmilw'
   if inparm eq 8 then print, 'directsw'
   if inparm eq 9 then print, 'tauintsw'
   if inparm eq 10 then print, 'tauintlw'

   if inparm eq 1 then title='mu0'
   if inparm eq 2 then title= 'Albedo'
   if inparm eq 3 then title= 'Albedo'
   if inparm eq 4 then title= 'IRflux'
   if inparm eq 5 then title= 'pi0s'
   if inparm eq 6 then title= 'tmisw'
   if inparm eq 7 then title= 'tmilw'
   if inparm eq 8 then title= 'directsw'
   if inparm eq 9 then title= 'tauintsw'
   if inparm eq 10 then title= 'tauintlw'

   qtp=fltarr(nlon+1,nlat)
   qtp[0:nlon-1,*]=param
   qtp[nlon,*]=reform(param[0,*])

   if keyword_set(mu0corr) then begin
   title='Albedo '+cgsymbol('times')+' '+cgsymbol('mu')+'!I0'
   qtp[0:nlon-1,*]=param*(mu0>1e-5)
   qtp[nlon,*]=reform(param[0,*])*(mu0[0,*]>1e-5)
   endif
   
   flon=fltarr(nlon+1)
   flon[0:nlon-1]=lons ;rebin(lons,nlon,nlat)
   flon[nlon,*]=lons(0) ;reform(lons[0]+360.)
   flat=lats ;fltarr(nlat)
;   flat[0:nlon-1]=rebin(nlat)
;   flat[nlon,*]=reform(lats[0])
   qmin=min(qtp,max=qmax)
   if keyword_set(prange) then begin
    qmin=0.0;prange(0)
    qmax=1.0;prange(1)
   endif
if keyword_set(psfile) then begin
   Polyfill, [1,1,0,0,1], [1,.2,.2,1,1], /NORMAL, COLOR=1
   MAP_SET,latoff,lonoff,rotat,/ortho,/ISOTROPIC,/noerase,position=[0.08,0.4,0.48,0.9]
endif else MAP_SET,latoff,lonoff,rotat,/ortho,/ISOTROPIC 

nlevels=250.
cbottom=5
step=(qmax-qmin)/nlevels
mylevels=indgen(nlevels)*step+qmin
ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

;cgloadct,39
cgloadct,3
;gamma_ct,.2
cgcontour,qtp,flon,flat,/overplot,/cell_fill,/closed,levels=mylevels
;stop
;qtpcloud=(qtp-qtpclear)*100.

;qmin=min(qtpcloud,max=qmax)
;nlevels=250.
;cbottom=55
;step=(qmax-qmin)/nlevels
;mylevels=indgen(nlevels)*step+qmin
;ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
;cgcontour,qtpcloud,flon,flat,/overplot,/cell_fill,/closed,levels=mylevels
;if keyword_set(grid) then map_grid,/label,color=0
colorbar,range=[qmin,qmax],format='(F12.3)',charsize=1,color=255,$
         position=[0.10, 0.3, 0.46, 0.34],bottom=cbottom,$
               ncolors=255-cbottom,title=title

stop
if keyword_set(psfile) then begin
;   Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=1
   MAP_SET,latoff,lonoff,rotat,/ortho,/ISOTROPIC,/noerase,position=[0.52,0.4,0.92,0.9]
endif else MAP_SET,latoff,lonoff,rotat,/ortho,/ISOTROPIC

nlevels=250.
cbottom=5
   qmin=min(qtp,max=qmax)
   if keyword_set(prange) then begin
    qmin=prange(0)
    qmax=prange(1)
   endif
step=(qmax-qmin)/nlevels
mylevels=indgen(nlevels)*step+qmin
ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

;cgloadct,39
cgloadct,3
cgcontour,qtp,flon,flat,/overplot,/cell_fill,/closed,levels=mylevels
;stop
;qtpcloud=(qtp-qtpclear)*100.

;qmin=min(qtpcloud,max=qmax)
;nlevels=250.
;cbottom=55
;step=(qmax-qmin)/nlevels
;mylevels=indgen(nlevels)*step+qmin
;ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
;cgcontour,qtpcloud,flon,flat,/overplot,/cell_fill,/closed,levels=mylevels
if keyword_set(grid) then map_grid,/label,color=255,londel=30;GRID_LONGITUDE=15,grid_latitude=15,londel=15
colorbar,range=[qmin,qmax],format='(F12.3)',charsize=1,color=255,$
         position=[0.54, 0.3, 0.90, 0.34],bottom=cbottom,$
               ncolors=255-cbottom,title=title


;xyouts,0.735,0.93,'Temperature [K]',charsize=2,/normal,color=255
if keyword_set(psfile) then begin
   psclose
   !p.font=-1
endif

latedgep=findgen((nlat/2.)+1)/(nlat/2.)*89.
 latedges=fltarr(nlat+1)
 latedges(0:nlat/2)=reverse(latedgep)
 latedges(nlat/2+1:*)=latedgep(1:*)*(-1.)

desired=findgen(nlon)+.5
lonedgesx=interpolate(lons,desired)
lonedgesx(nlon-1)=360.-lonedgesx(0)
lonedges=findgen(nlon+1)
lonedges(1:*)=lonedgesx
lonedges(0)=-1.*lonedgesx(0)
latser=latedges*!dtor
lonser=lonedges*!dtor
r=1.
 gridarea=fltarr(nlon,nlat)
  for i=0,nlon-1 do begin
    for j=0,nlat-1 do begin
       gridarea(i,j)=r^2.*(lonser(i+1)-lonser(i))*(sin(latser(j))-sin(latser(j+1)))
       outgoing(i,j)=toafluxes(i,j,0)*gridarea(i,j)
       outgoingIR(i,j)=IRflux(i,j)*gridarea(i,j)
       incoming(i,j)=toafluxes(i,j,1)*gridarea(i,j)   
    endfor
  endfor
       spherical_alb= total(outgoing)/total(incoming) 
       fracgrid=gridarea/(total(gridarea)/1.)
   spheralb=total(fracgrid*toaa)*2. 
;print,'spherical albedo old: ',spheralb

print,'Spherical albedo: ',spherical_alb 
print,'Outging (IR+vis):',total(outgoingIR)+total(outgoing)
print,'incoming (vis):',total(incoming)
print,'difference (Out - in):',total(outgoingIR)+total(outgoing)-total(incoming)
print,'percent difference (out-in)/in:',(total(incoming)-(total(outgoingIR)+total(outgoing)))/(total(outgoingIR)+total(outgoing))*100.

print,'outgoing/ingoing',(total(outgoingIR)+total(outgoing))/(total(incoming))*100.
print,'param average:', qtp*fracgrid
stop


end

