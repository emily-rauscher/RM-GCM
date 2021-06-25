pro apegen,levplot,qplot=qplot,psfile=psfile

;qplot=0 (not set): plot q_nf*(1-T_m/T)
;qplot=1: plot (1-T_m/T)
;qplot=2: plot q_nf
;qplot=3: plot ke
;qplot=4: plot log(abs(tdrag)) for Q_fr=abs(0.1*Q_nf)
;qplot=5: plot q_fr*(1-T_m/T) for uniform drag with taudrag (set below)
;
;NEED TO HAVE FORT.26, FORT.54, FORT.50, AND ARRAYS SET BELOW

if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0,retain=2
;  set_plot, 'x'
  window, 0, xsize=x_size_window, ysize=y_size_window, retain=2
  !p.font=-1
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20;,bits_per_pixel=24
  !p.font=0
endelse

  !x.style=1
  !p.charsize=1.5

;for R58w:
resttt=[511.,524.,539.,555.,570.,585.,598.,618.,653.,690.,$
730.,769.,817.,890.,963.,1049.,1125.,1190.,1258.,1299.,$
1345.,1387.,1437.,1477.,1497.,1780.,1782.,1787.,1791.,1793.,$
1796.,1801.,1808.]
redtep=[replicate(1000.,12),994.,956.,918.,$
880.,842.,804.,766.,728.,690.,652.,614.,576.,538.,replicate(0.,8)]
RESTIM = [0.014,0.017,0.019,0.022,0.025,0.031,0.040,0.049,0.060,0.067,$
0.074,0.080,0.11,0.15,0.21,0.29,0.40,0.54,0.67,1.7,$
4.3,9.7,22.,48.,120.,replicate(1.e30,8)]
GA = 9.42
GASCON = 4593.
RADEA = 9.44e7
AKAP = 0.321
WW = 2.06e-5
OOM = 5.34363
p0=220.612
restim*=(2.*!pi/ww)
taudrag=(2.*!pi/ww)*100.        ;uniform drag times at P_orb/100.

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17
  lon=reform(xy[0,*,*,levplot]) & lat=reform(xy[1,*,*,levplot]) 
  u=reform(xy[3,*,*,levplot]) & v=reform(xy[4,*,*,levplot]) 
  temp=reform(xy[5,*,*,levplot])
  t3d=reform(xy[5,*,*,*])

  OpenR, 18, 'fort.54'
  ReadF, 18, nlat,nlon,nlev	  
	  ab=fltarr(4,nlat,nlon,nlev)   
  ReadF, 18, ab
  Close, 18
  ke=reform(ab[3,*,*,levplot])
  ke*=p0*1.e5

  OpenR, 19, 'fort.50'
  ReadF, 19, nlat,nlon	  
	  ij=fltarr(3,nlat,nlon)   
  ReadF, 19, ij
  Close, 19
  sp=reform(ij[2,*,*])
  sp=(sp+1.)*p0*1.e5

  sigma=get_sigma(oom,nlev)

  flon=lon*360./lon(0,nlon-1)
  flat=lat*90./lat(0,nlat-1)

;3D mass elements (boundary calcs rough)
dp=fltarr(nlat,nlon,nlev)
dm=fltarr(nlat,nlon,nlev)
for i=1,nlev-2 do dp[*,*,i]=sp*(sigma[i+1]-sigma[i-1])/2.
dp[*,*,0]=sp*sigma[1]/2. ;boundary is sig=0
dp[*,*,nlev-1]=sp*(1.-sigma[nlev-2])/2. ;boundary is sig=1
for i=0,nlev-1 do dm[*,*,i]=dp[*,*,i]*(1./ga)*radea^2 $
   *cos(lat*!pi/180.)*(2.*!pi/nlon)*(!pi/nlat)

;3D q_nf
teq=fltarr(nlat,nlon,nlev)
qnf=fltarr(nlat,nlon,nlev)
nonheat=where(redtep eq 0)  ;for inert lower levels
if nonheat[0] eq -1 then heatind=nlev-1 else heatind=nonheat[0]-1
for i=0,heatind do begin
   for j=0,nlon-1 do begin
      if (lon[0,j] ge 90.) and (lon[0,j] le 270.) then begin
         teq[*,j,i]=resttt[i]
         qnf[*,j,i]=dm[*,j,i]*(gascon/akap)*$
                    (teq[*,j,i]-t3d[*,j,i])/restim[i]
      endif else begin
         for k=0,nlat-1 do begin
            teq[k,j,i]=( (resttt[i])^4 + $
                         ((resttt[i]+redtep[i])^4 $
                          -(resttt[i])^4) $
                         *cos(lat[k,j]*!pi/180.)*cos(lon[k,j]*!pi/180.) $
                       )^0.25
            qnf[k,j,i]=dm[k,j,i]*(gascon/akap)*$
                       (teq[k,j,i]-t3d[k,j,i])/restim[i]
         endfor
      endelse
   endfor
endfor
if nonheat[0] ne -1 then begin
   for i=nonheat[0],nlev-1 do begin
      teq[*,*,i]=resttt[i]
      qnf[*,*,i]=0.
   endfor
endif

;1.-(T_m/T) 3D efficiency factor
tm=1./(total(dm[*,*,0:heatind]/t3d[*,*,0:heatind])/total(dm[*,*,0:heatind]))
print,'T_m:',tm
print,'atm int of efficiency:',total((1.-tm/t3d[*,*,0:heatind])*dm[*,*,0:heatind])/total(dm[*,*,0:heatind]) ;should =0
tvar=1.-tm/t3d

print,'atm net nf ape gen [W]:',total(tvar*qnf)
print,'atm net nf heating [W]:',total(qnf)

  loadct,4

;psopen,'plot',/enc,/color;,/landscape
;!P.MULTI=[0,2,2]
;quant=reform((tvar)[*,0,*])
;contour,quant,reform(lat[*,0]),sigma,/cell_fill,/closed,nlevels=45,$
;        yr=[1.,min(sigma)],/ylog,ytit='lon=0'
;colorbar,position=[0.05,0.93,0.45,0.95],range=[min(quant),max(quant)],$
;         format='(e8.1)',divisions=4;,charsize=1.5
;contour,quant,reform(lat[*,0]),sigma,/overplot,levels=[0.],color=255,thick=4
;quant=reform((tvar)[*,24,*])
;contour,quant,reform(lat[*,0]),sigma,/cell_fill,/closed,nlevels=45,$
;        yr=[1.,min(sigma)],/ylog,ytit='lon=90'
;colorbar,position=[0.55,0.93,0.95,0.95],range=[min(quant),max(quant)],$
;         format='(e8.1)',divisions=4;,charsize=1.5
;contour,quant,reform(lat[*,0]),sigma,/overplot,levels=[0.],color=255,thick=4
;quant=reform((tvar)[*,48,*])
;contour,quant,reform(lat[*,0]),sigma,/cell_fill,/closed,nlevels=45,$
;        yr=[1.,min(sigma)],/ylog,ytit='lon=180'
;colorbar,position=[0.05,0.5,0.45,0.53],range=[min(quant),max(quant)],$
;         format='(e8.1)',divisions=4;,charsize=1.5
;contour,quant,reform(lat[*,0]),sigma,/overplot,levels=[0.],color=255,thick=4
;quant=reform((tvar)[*,72,*])
;contour,quant,reform(lat[*,0]),sigma,/cell_fill,/closed,nlevels=45,$
;        yr=[1.,min(sigma)],/ylog,ytit='lon=270'
;colorbar,position=[0.55,0.5,0.95,0.53],range=[min(quant),max(quant)],$
;         format='(e8.1)',divisions=4;,charsize=1.5
;contour,quant,reform(lat[*,0]),sigma,/overplot,levels=[0.],color=255,thick=4
;psclose
;!P.MULTI=0

;stop

  MAP_SET, /Miller_CYLINDRICAL,0,0,/ISOTROPIC ;,TITLE=' Temperature 

if keyword_set(qplot) then begin
   if qplot eq 1 then begin
      contour,tvar[*,*,levplot],flon,flat,/overplot,/cell_fill,/closed,nlevels=45
      colorbar,position=[0.1,0.07,0.90,0.1],range=[min(tvar[*,*,levplot]),max(tvar[*,*,levplot])],$
               format='(f6.3)',charsize=2
      contour,tvar[*,*,levplot],flon,flat,/overplot,levels=[0.],color=255,thick=5
      ;print,minmax(tvar)
   endif
   if qplot eq 2 then begin
      contour,qnf[*,*,levplot],flon,flat,/overplot,/cell_fill,/closed,nlevels=45
      colorbar,position=[0.1,0.07,0.90,0.1],range=[min(qnf[*,*,levplot]),max(qnf[*,*,levplot])],$
               format='(e8.1)',charsize=2
      contour,qnf[*,*,levplot],flon,flat,/overplot,levels=[0.],color=255,thick=5
      ;print,minmax(qnf)
   endif
   if qplot eq 3 then begin
      contour,ke,flon,flat,/overplot,/cell_fill,/closed,nlevels=45
      colorbar,position=[0.1,0.07,0.90,0.1],range=[min(ke),max(ke)],$
               format='(e8.1)',charsize=2
      contour,ke,flon,flat,/overplot,levels=[0.],color=255,thick=5
      print,minmax(ke)
   endif
   if qplot eq 4 then begin
      tdrag=10.*(ke/qnf[*,*,levplot]) ;drag times to get dissipation=0.1*agen
      ;(negative drag time means works against heating)
      tdraglog=alog10(abs(tdrag))
      contour,tdraglog,flon,flat,/overplot,/cell_fill,/closed,nlevels=45
      colorbar,position=[0.1,0.07,0.90,0.1],range=[min(tdraglog),max(tdraglog)],$
               format='(f5.2)',charsize=2
      contour,tdrag,flon,flat,/overplot,levels=[0.],color=255,thick=5
      print,minmax(tdrag)
   endif
   if qplot eq 5 then begin
      agenfr=(ke/taudrag)*tvar[*,*,levplot]
      contour,agenfr,flon,flat,/overplot,/cell_fill,/closed,nlevels=45
      colorbar,position=[0.1,0.07,0.90,0.1],range=[min(agenfr),max(agenfr)],$
               format='(e8.1)',charsize=2
      contour,agenfr,flon,flat,/overplot,levels=[0.],color=255,thick=5
      print,minmax(agenfr)
   endif
endif else begin
   agen=tvar[*,*,levplot]*qnf[*,*,levplot]
   contour,agen,flon,flat,/overplot,/cell_fill,/closed,nlevels=45
   colorbar,position=[0.1,0.07,0.90,0.1],range=[min(agen),max(agen)],$
            format='(e8.1)',charsize=2
   contour,agen,flon,flat,/overplot,levels=[0.],color=255,thick=5
   print,minmax(agen)
endelse

  map_grid,/label,color=255

  print,'pressure:',p0*sigma[levplot]

if not keyword_set(psfile) then begin
  set_plot, 'x'
endif else begin
   psclose
endelse

stop

  loadct,0

end
