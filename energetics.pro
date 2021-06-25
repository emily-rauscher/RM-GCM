pro energetics,p0,bfield,ga,radea,specific=specific,paperplot=paperplot,horizontal=horizontal,old=old
;to calculate the available and unavailable enthalpies, the
;kinetic energies, and the conversion rates locally throughout the
;model

;actual, use keyword for specific

;(so far) ignoring variation in p_surf

;  device,true_color=24,decomposed=0,retain=2

if keyword_set(old) then begin
   ;for run 58w (or 58y, or r64-66):
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
   WW = 2.06e-5                 ;rad/sec
   OOM = 5.34363
   p0=220.612
   restim*=(2.*!pi/ww)          ;convert to seconds
   tgr=1813.
   tfrc=fltarr(33)

   ;for runs with mag drag:
   tfrc=[6.6,19,32,44,56,69,81,94,110,120,130,140,150,160,170,180,190,200,210,$
         220,230,470,710,950,1200,1400,1700,1900,2100,2400,2600,2900,3100.] ;r64
   ;tfrc/=10. ;r65
   tfrc/=100.                   ;r66
   drag=1                       ;0 for Cua for RM10, 1 for CuDaD for PMR

   ;tfrc*=(2.*!pi/ww)  ;convert to seconds
endif

gascon=3523.
akap=2./7

OpenR, 17, 'fort.26'
ReadF, 17, nlat,nlon,nlev	  
   xy=fltarr(6,nlat,nlon,nlev)   
ReadF, 17, xy
Close, 17
lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) 
temp=reform(xy[5,*,*,*])

sigma=get_sigma(oom,nlev)
pres=p0*sigma
dP=fltarr(nlev)
for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1]) ; to match swdp
dP[0]=0.5*(pres[1])
dP[nlev-1]=0.5*(p0-pres[nlev-2])
dmass=dP[nl]*1.e5/ga*radea^2    ; kg

if not keyword_set(old) then begin

   openr, 18, 'fort.50'
   readf, 18, nlat,nlon
   ab=fltarr(3,nlat,nlon)
   readf, 18, ab
   close, 18
   sp=reform(ab[2,*,*])
   sp = (sp+1.)*p0*1.e5         ;SI

   rhocgs=fltarr(nlat,nlon,nlev)
   tdrag=fltarr(nlat,nlon,nlev)
   udrag=fltarr(nlat,nlon,nlev)
   heating=fltarr(nlat,nlon,nlev)

   for ilat=0,nlat-1 do begin
      ang_f=abs(cos(!pi/2.-xy[1,ilat,0,0]*!pi/180.))
      for ilev=0,nlev-1 do begin
         rhocgs[ilat,*,ilev]=reform(sp[ilat,*]*sigma[ilev]/gascon/xy[5,ilat,*,ilev])*1.e-3 ;cgs
         tdrag[ilat,*,ilev]=tmag(reform(rhocgs[ilat,*,ilev]),reform(xy[5,ilat,*,ilev]),$
                                 bfield,simp=usesimp)/ang_f ; in sec
      endfor
   endfor
endif

;mass elements
dm=fltarr(nlat,nlon,nlev)
for ilev=0,nlev-1 do begin
   for ilat=0,nlat-1 do begin
      dm[ilat,*,ilev]=replicate(cos(lat[ilat,0,0]*!pi/180.),nlon)*(2.*!pi/nlon)*(!pi/nlat)*dmass[ilev]
   endfor
endfor

;3D q_nf
if keyword_set(old) then begin
   teq=fltarr(nlat,nlon,nlev)
   qnf=fltarr(nlat,nlon,nlev)
   nonheat=where(redtep eq 0)   ;for inert lower levels
   if nonheat[0] eq -1 then heatind=nlev-1 else heatind=nonheat[0]-1
   for i=0,heatind do begin
      for j=0,nlon-1 do begin
         if (lon[0,j,i] ge 90.) and (lon[0,j,i] le 270.) then begin
            teq[*,j,i]=resttt[i]
            qnf[*,j,i]=dm[*,j,i]*(gascon/akap)*$
                       (teq[*,j,i]-temp[*,j,i])/restim[i]
         endif else begin
            for k=0,nlat-1 do begin
               teq[k,j,i]=( (resttt[i])^4 + $
                            ((resttt[i]+redtep[i])^4 $
                             -(resttt[i])^4) $
                            *cos(lat[k,j,i]*!pi/180.)*cos(lon[k,j,i]*!pi/180.) $
                          )^0.25
               qnf[k,j,i]=dm[k,j,i]*(gascon/akap)*$
                          (teq[k,j,i]-temp[k,j,i])/restim[i]
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
endif else begin
   heatind=nlev-1
   ; NEED TO PULL INFO OUT OF FORT.63, MAYBE FROM HEATING.PRO?
endelse
print,'Total radiative heating in J/s:',total(qnf)

;Carnot efficiency T_w and T_c:
qpos=where(qnf gt 0.,complement=qneg)
tw=total(temp[qpos]*qnf[qpos])/total(qnf[qpos])
tc=total(temp[qneg]*qnf[qneg])/total(qnf[qneg])
print,'Carnot T_w, T_c, eta:',tw,tc,(tw-tc)/tw

;1.-(T_m/T) 3D efficiency factor
;tm=1./(total(dm[*,*,0:heatind]/temp[*,*,0:heatind])/total(dm[*,*,0:heatind]))
tm=1./(total(dm/temp)/total(dm))
print,'T_m:',tm
print,'atm int of efficiency:',total((1.-tm/temp[*,*,0:heatind])*dm[*,*,0:heatind])/total(dm[*,*,0:heatind]) ;should =0
;tvar=1.-tm/temp

;energies: (derivatives will be for total, not just T component, of u
;and a, but will keep this anyway)
uT=(gascon/akap)*tm*alog(temp/tm)
aT=(gascon/akap)*[temp-tm]-uT
;u=(gascon/akap)*tm*alog(temp/tm)
ke=0.5*(u^2+v^2)


;integrate hydrostatic for geopotential:
;gz=fltarr(nlat,nlon,nlev)
;set altitude of first level (up from base=p0, where T=TGR)
;gz[*,*,nlev-1]=gascon*0.5*(temp[*,*,nlev-1]+TGR)*alog(1./sigma[nlev-1])
;integrate hydrostatic to solve for higher levels
;for i=nlev-2,0,-1 do begin
;   gz[*,*,i]=gz[*,*,i+1]+gascon*0.5*(temp[*,*,i]+temp[*,*,i+1])*alog(sigma[i+1]/sigma[i])
;endfor

;integrate continuity for omega (=0 at p=0)
; careful: increasing the lat index is going south, not north
om=fltarr(nlat,nlon,nlev)
for ilev=0,nlev-1 do begin
   if ilev gt 0 then dP=p0*(sigma[ilev]-sigma[ilev-1]) else dP=p0*sigma[0]
   for ilon=0,nlon-1 do begin
      ilo=ilon-1
      ilonn=ilon+1
      if ilon eq 0 then ilo=nlon-1
      if ilon eq nlon-1 then ilonn=0
      ;first do poles (should be ilat -> ilat,ilon+nlon/2), with flow now negative)
      ;careful with +nlon/2
      if ilon+nlon/2 gt nlon-1 then iilon=ilon-nlon/2 else iilon=ilon+nlon/2
      if ilev gt 0 then begin
         om[0,ilon,ilev]=(-0.5/radea)$
                          *[(u[0,ilonn,ilev]-u[0,ilo,ilev])/(cos(lat[0,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(-1.*v[0,iilon,ilev]-v[1,ilon,ilev])/(!pi/nlat)] *dP + om[0,ilon,ilev-1]
         om[nlat-1,ilon,ilev]=(-0.5/radea)$
                          *[(u[nlat-1,ilonn,ilev]-u[nlat-1,ilo,ilev])$
                            /(cos(lat[nlat-1,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[nlat-2,ilon,ilev]+v[nlat-1,iilon,ilev])/(!pi/nlat)] *dP + om[nlat-1,ilon,ilev-1]
      endif else begin
         om[0,ilon,ilev]=(-0.5/radea) $
                          *[(u[0,ilonn,ilev]-u[0,ilo,ilev])/(cos(lat[0,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(-1.*v[0,iilon,ilev]-v[1,ilon,ilev])/(!pi/nlat)] *dP 
         om[nlat-1,ilon,ilev]=(-0.5/radea) $
                          *[(u[nlat-1,ilonn,ilev]-u[nlat-1,ilo,ilev])$
                            /(cos(lat[nlat-1,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[nlat-2,ilon,ilev]+v[nlat-1,iilon,ilev])/(!pi/nlat)] *dP 
      endelse
      for ilat=1,nlat-2 do begin
         if ilev gt 0 then om[ilat,ilon,ilev]=(-0.5/radea)$
                          *[(u[ilat,ilonn,ilev]-u[ilat,ilo,ilev])/(cos(lat[ilat,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[ilat-1,ilon,ilev]-v[ilat+1,ilon,ilev])/(!pi/nlat)] *dP + om[ilat,ilon,ilev-1] $
         else om[ilat,ilon,ilev]=(-0.5/radea) $
                          *[(u[ilat,ilonn,ilev]-u[ilat,ilo,ilev])/(cos(lat[ilat,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[ilat-1,ilon,ilev]-v[ilat+1,ilon,ilev])/(!pi/nlat)] *dP 
      endfor
   endfor
endfor
om*=1.e5  ;to convert units to: (kg/m/s^2)/s


;creations (non frictional)
Cu=(tm/temp)*qnf
Ca=qnf-Cu

;conversions
Cake=fltarr(nlat,nlon,nlev) ;=omega (RT/P)
;CuTaT=om*tm*gascon*1.e-5  ;for units
for ilev=0,nlev-1 do begin
   Cake[*,*,ilev]=gascon*temp[*,*,ilev]*om[*,*,ilev]/(p0*sigma[ilev]*1.e5)
;   CuTaT[*,*,ilev]=CuTaT[*,*,ilev]/(p0*sigma[ilev])
;   for ilon=0,nlon-1 do begin
;      ilo=ilon-1
;      ilonn=ilon+1
;      if ilon eq 0 then ilo=nlon-1
;      if ilon eq nlon-1 then ilonn=0
;      if ilon+nlon/2 gt nlon-1 then iilon=ilon-nlon/2 else iilon=ilon+nlon/2
      ;poles first (+nlon/2, geopotential [now omega] so no sign change)
      ;CaTke[0,ilon,ilev]=(0.5/radea)*[u[0,ilon,ilev]*(gz[0,ilonn,ilev]-gz[0,ilo,ilev])$
      ;                                /(cos(lat[0,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
      ;                                +v[0,ilon,ilev]*(gz[0,iilon,ilev]-gz[1,ilon,ilev])/(!pi/nlat)]
      ;CaTke[nlat-1,ilon,ilev]=(0.5/radea)*[u[nlat-1,ilon,ilev]*(gz[nlat-1,ilonn,ilev]-gz[nlat-1,ilo,ilev])$
      ;                                     /(cos(lat[nlat-1,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
      ;                                     +v[nlat-1,ilon,ilev]*(gz[nlat-2,ilon,ilev]-gz[nlat-1,iilon,ilev])/(!pi/nlat)]
;      for ilat=1,nlat-2 do begin
         ;CaTke[ilat,ilon,ilev]=(0.5/radea)*[u[ilat,ilon,ilev]*(gz[ilat,ilonn,ilev]-gz[ilat,ilo,ilev])$
         ;                                   /(cos(lat[ilat,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
         ;                                   +v[ilat,ilon,ilev]*(gz[ilat-1,ilon,ilev]-gz[ilat+1,ilon,ilev])/(!pi/nlat)]
;      endfor
;   endfor
endfor

;dissipation
test=where(tfrc ne 0)
if test[0] ne -1 then begin
   ds=fltarr(nlat,nlon,nlev)
   for i=0,nlev-1 do ds[*,*,i]=ke[*,*,i]/tfrc[i]
endif else ds=fltarr(nlat,nlon,nlev)
;contribution to u and a
CuD=(tm/temp)*ds
CaD=(1.-tm/temp)*ds

if not keyword_set(specific) then begin
   uT*=dm
   aT*=dm
   ke*=dm
;   Cu*=dm  (already in units of J/s)
;   Ca*=dm  (already in units of J/s)
;   CuTaT*=dm
   Cake*=dm
   ds*=dm
   CuD*=dm
   CaD*=dm
endif else begin
   Cu/=dm
   Ca/=dm
endelse

;CuT+=CuTD
;CaT+=CaTD

print,'Total energies (uT, aT, v^2) in J:',total(ut),total(aT),total(ke)
print,'Total rates of u, a source terms in J/s:',total(Cu),total(Ca)
print,'Total conversion rate of aT->ke in J/s:',-1.*total(Cake)
print,'Total rate of dissipation (v^2)/tfrc in J/s:',total(ds)
print,'Total contribution dissipation would have to u and a in J/s:',total(CuD),total(CaD)
print,''

vert=fltarr(nlat,nlon,nlev)
for i=0,nlat-1 do begin
   for j=0,nlon-1 do begin
      vert[i,j,*]=p0*sigma
   endfor
endfor

night=where(lon[*,*,*] ge 90. and lon[*,*,*] lt 270.)
upper=where(vert[0,0,*] lt 1.)

print,'day vs. night: '
toss=total(Ca[night])
print,'  APE, ',total(Ca)-toss,toss
toss=total(Cu[night])
print,'  UPE gen:',total(Cu)-toss,toss
toss=total(ds[night])
print,'  diss:',total(ds)-toss,toss
toss=total(CaD[night])
print,'  APE, ',total(CaD)-toss,toss
toss=total(CuD[night])
print,'  UPE diss gen:',total(CuD)-toss,toss

print,'upper (above 1 bar) vs. lower: '
toss=total(Ca[*,*,upper])
print,'  APE, ',toss,total(Ca)-toss
toss=total(Cu[*,*,upper])
print,'  UPE gen:',toss,total(Cu)-toss
toss=total(ds[*,*,upper])
print,'  diss:',toss,total(ds)-toss
toss=total(CaD[*,*,upper])
print,'  APE, ',toss,total(CaD)-toss
toss=total(CuD[*,*,upper])
print,'  UPE diss gen:',toss,total(CuD)-toss

if keyword_set(horizontal) then begin

;device,true_color=24,decomposed=0,retain=2

psopen,'Cua',/enc,/color,bits_per_pixel=8
!p.charsize=1.5
!p.font=0
!p.thick=4

loadct,4
;   nl=13
nl=horizontal
   print,'P=',p0*sigma[nl],'bar'

  flon=reform(lon[*,*,0])*360./lon[0,nlon-1,0]
  flat=reform(lat[*,*,0])*90./lat[nlat-1,0,0]

if drag eq 0 then begin

umax=max(Cu[*,*,nl],min=umin)
if abs(umax) gt abs(umin) then umin=-1.*umax else umax=-1.*umin
nlevels=45
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
MAP_SET, /Miller_CYLINDRICAL,0.,180.,/ISOTROPIC,position=[0.1,0.5,0.9,0.90]
contour,reform(Cu[*,*,nl]),flon,flat,/cell_fill,levels=mylevels,min_value=umin,max_value=umax,/overplot,/closed
map_grid,/label,charsize=1,latlab=355.;,color=255
colorbar,range=[umin,umax],position=[0.15,0.52,0.85,0.54],format='(e8.1)',divisions=4,charsize=1.25
xyouts,10.,70.,'a) d(UPE)/dt',charsize=1.5,color=255
umax=max(Ca[*,*,nl],min=umin)
if abs(umax) gt abs(umin) then umin=-1.*umax else umax=-1.*umin
nlevels=45
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
MAP_SET, /Miller_CYLINDRICAL,0.,180.,/ISOTROPIC,position=[0.1,0.1,0.9,0.50],/noerase
contour,reform(Ca[*,*,nl]),flon,flat,/cell_fill,levels=mylevels,min_value=umin,max_value=umax,/overplot,/closed
map_grid,/label,charsize=1,latlab=355.;,color=255
colorbar,range=[umin,umax],position=[0.15,0.12,0.85,0.14],format='(e8.1)',divisions=4,charsize=1.25
xyouts,10.,70.,'b) d(APE)/dt',charsize=1.5,color=255

endif else begin

umax=max(CuD[*,*,nl],min=umin)
absmax=max([abs(umax),abs(umin)])
range=umax-umin
maxrange=2.*absmax  ;from -absmax to +absmax
rat=range/maxrange
cbottom=(absmax+umin)/maxrange*255.
nlevels=44.*rat+1.  ;so that delta color always = 255./44.
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
MAP_SET, /Miller_CYLINDRICAL,0.,180.,/ISOTROPIC,position=[0.1,0.5,0.9,0.90]
contour,reform(CuD[*,*,nl]),flon,flat,/cell_fill,min_value=umin,max_value=umax,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels,/overplot,/closed
map_grid,/label,charsize=1,latlab=355.;,color=255
colorbar,range=[umin,umax],position=[0.15,0.52,0.85,0.54],format='(e8.1)',bottom=cbottom,ncolors=256*rat,divisions=4,charsize=1.25
xyouts,10.,70.,'a) d(UPE)/dt',charsize=1.5,color=255
umax=max(CaD[*,*,nl],min=umin)
absmax=max([abs(umax),abs(umin)])
range=umax-umin
maxrange=2.*absmax  ;from -absmax to +absmax
rat=range/maxrange
cbottom=(absmax+umin)/maxrange*255.
nlevels=44.*rat+1.  ;so that delta color always = 255./44.
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
MAP_SET, /Miller_CYLINDRICAL,0.,180.,/ISOTROPIC,position=[0.1,0.1,0.9,0.50],/noerase
contour,reform(CaD[*,*,nl]),flon,flat,/cell_fill,min_value=umin,max_value=umax,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels,/overplot,/closed
map_grid,/label,charsize=1,latlab=355.;,color=255
colorbar,range=[umin,umax],position=[0.15,0.12,0.85,0.14],format='(e8.1)',bottom=cbottom,ncolors=256*rat,divisions=4,charsize=1.25
xyouts,10.,70.,'b) d(APE)/dt',charsize=1.5,color=255

endelse

;window,0
;umin=min(CuD[*,*,nl],max=umax)
;  MAP_SET, /Miller_CYLINDRICAL,0.,0.,/ISOTROPIC,title='UPE'
;  contour, reform(CuD[*,*,nl]),flon,flat, /overplot,/cell_fill,/closed ,NLEVELS=45 ;, min_value=1000, max_value=2000
;  colorbar,  position=[0.1,0.07,0.90,0.1], format='(e8.1)', charsize=2,range=[umin,umax]
;  map_grid, /label ,color=255
;window,1
;umin=min(CaD[*,*,nl],max=umax)
;  MAP_SET, /Miller_CYLINDRICAL,0.,0.,/ISOTROPIC,title='APE'
;  contour, reform(CaD[*,*,nl]),flon,flat, /overplot,/cell_fill,/closed ,NLEVELS=45 ;, min_value=1000, max_value=2000
;  colorbar,  position=[0.1,0.07,0.90,0.1], format='(e8.1)', charsize=2,range=[umin,umax]
;  map_grid, /label ,color=255
;loadct,0

loadct,0
psclose

endif

if keyword_set(paperplot) then begin

;  device,true_color=24,decomposed=0,retain=2
;  window, 1, xsize=600, ysize=800, retain=2

loadct,4


print,minmax(Cu),minmax(Ca)

psopen,'Cua',/enc,/color,bits_per_pixel=8
!p.charsize=1.5
!p.font=0
!p.thick=4
!y.tickname=['100','10','1','0.1','0.01']

if drag eq 0 then begin

umax=max(Cu,min=umin)
if abs(umax) gt abs(umin) then umin=-1.*umax else umax=-1.*umin
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cu[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='Pressure [bar]',min_value=umin,max_value=umax,/xstyle,xtickformat="(A1)",position=[0.16,0.5,0.99,0.9]
colorbar,range=[umin,umax],position=[0.23,0.55,0.93,0.58],format='(e8.1)',divisions=4,charsize=1.25,color=255
xyouts,10.,3.e-3,'a) d(UPE)/dt',color=255
umax=max(Ca,min=umin)
if abs(umax) gt abs(umin) then umin=-1.*umax else umax=-1.*umin
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Ca[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='Pressure [bar]',xtitle='longitude [degrees east of substellar]',min_value=umin,max_value=umax,/xstyle,position=[0.16,0.1,0.99,0.5],/noerase
colorbar,range=[umin,umax],position=[0.23,0.15,0.93,0.18],format='(e8.1)',divisions=4,charsize=1.25,color=255
xyouts,10.,3.e-3,'b) d(APE)/dt',color=255

endif else begin

umax=max(CuD,min=umin)
absmax=max([abs(umax),abs(umin)])
range=umax-umin
maxrange=2.*absmax  ;from -absmax to +absmax
rat=range/maxrange
cbottom=(absmax+umin)/maxrange*255.
nlevels=24.*rat+1.  ;so that delta color always = 255./24.
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CuD[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='Pressure [bar]',min_value=umin,max_value=umax,/xstyle,xtickformat="(A1)",position=[0.16,0.5,0.99,0.9],c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom
colorbar,range=[umin,umax],position=[0.23,0.55,0.93,0.58],format='(e8.1)',bottom=cbottom,ncolors=256*rat,divisions=4,charsize=1.25,color=255
xyouts,10.,3.e-3,'a) d(UPE)/dt',color=255
umax=max(CaD,min=umin)
absmax=max([abs(umax),abs(umin)])
range=umax-umin
maxrange=2.*absmax  ;from -absmax to +absmax
rat=range/maxrange
nlevels=24.*rat+1.  ;so that delta color always = 255./24.
cbottom=(absmax+umin)/maxrange*255.
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CaD[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='Pressure [bar]',xtitle='longitude [degrees east of substellar]',min_value=umin,max_value=umax,/xstyle,position=[0.16,0.1,0.99,0.5],/noerase,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom
colorbar,range=[umin,umax],position=[0.23,0.15,0.93,0.18],format='(e8.1)',bottom=cbottom,ncolors=256*rat,divisions=4,charsize=1.25,color=255
xyouts,10.,3.e-3,'b) d(APE)/dt',color=255
endelse

psclose

;        write_png,'Cua.png',tvrd(true=1)

endif else begin

loadct,4

;plots of qnf
psopen,'qnf_lons',/enc,/color
!P.MULTI=[0,2,2]
umax=max(qnf[*,0,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(qnf[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='qnf',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(qnf[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(qnf[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='qnf',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(qnf[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(qnf[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='qnf',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(qnf[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(qnf[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='qnf',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'qnf_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(qnf[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(qnf[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='qnf',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(qnf[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(qnf[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='qnf',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose
;plots of uT
;psopen,'uT_lons',/enc,/color
;!P.MULTI=[0,2,2]
;umax=max(uT[*,0,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(uT[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='uT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
;umax=max(uT[*,24,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(uT[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='uT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
;umax=max(uT[*,48,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(uT[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='uT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
;umax=max(uT[*,72,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(uT[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='uT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
;psclose
;psopen,'uT_lats',/enc,/color;,/landscape
;!P.MULTI=[0,1,2]
;umax=max(uT[12,*,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(uT[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='uT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
;umax=max(uT[24,*,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(uT[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='uT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
;psclose
;plots of aT
;psopen,'aT_lons',/enc,/color
;!P.MULTI=[0,2,2]
;umax=max(aT[*,0,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(aT[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='aT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
;umax=max(aT[*,24,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(aT[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='aT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
;umax=max(aT[*,48,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(aT[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='aT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
;umax=max(aT[*,72,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(aT[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='aT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
;psclose
;psopen,'aT_lats',/enc,/color;,/landscape
;!P.MULTI=[0,1,2]
;umax=max(aT[12,*,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(aT[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='aT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
;umax=max(aT[24,*,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(aT[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='aT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
;psclose
;plots of ke
psopen,'ke_lons',/enc,/color
!P.MULTI=[0,2,2]
umax=max(ke[*,0,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ke[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='ke',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(ke[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ke[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='ke',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(ke[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ke[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='ke',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(ke[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ke[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='ke',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'ke_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(ke[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ke[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='ke',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(ke[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ke[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='ke',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose

;plots of Cu
psopen,'Cu_lons',/enc,/color
!P.MULTI=[0,2,2]
umax=max(Cu[*,0,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cu[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='Cu',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(Cu[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cu[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='Cu',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(Cu[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cu[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='Cu',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(Cu[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cu[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='Cu',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'Cu_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(Cu[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cu[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='Cu',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(Cu[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cu[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='Cu',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose
;plots of Ca
psopen,'Ca_lons',/enc,/color
!P.MULTI=[0,2,2]
umax=max(Ca[*,0,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Ca[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='Ca',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(Ca[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Ca[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='Ca',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(Ca[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Ca[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='Ca',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(Ca[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Ca[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='Ca',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'Ca_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(Ca[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Ca[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='Ca',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(Ca[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Ca[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='Ca',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose
;plots of CuTaT
;psopen,'CuTaT_lons',/enc,/color
;!P.MULTI=[0,2,2]
;umax=max(CuTaT[*,0,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(CuTaT[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='CuTaT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
;umax=max(CuTaT[*,24,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(CuTaT[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='CuTaT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
;umax=max(CuTaT[*,48,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(CuTaT[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='CuTaT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
;umax=max(CuTaT[*,72,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(CuTaT[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='CuTaT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
;psclose
;psopen,'CuTaT_lats',/enc,/color;,/landscape
;!P.MULTI=[0,1,2]
;umax=max(CuTaT[12,*,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(CuTaT[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='CuTaT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
;umax=max(CuTaT[24,*,*],min=umin)
;nlevels=25
;step=(umax-umin)/nlevels
;mylevels=indgen(nlevels)*step+umin
;contour,reform(CuTaT[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='CuTaT',min_value=umin,max_value=umax,/xstyle
;colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
;psclose
;plots of Cake
psopen,'Cake_lons',/enc,/color
!P.MULTI=[0,2,2]
umax=max(Cake[*,0,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cake[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='Cake',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(Cake[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cake[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='Cake',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(Cake[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cake[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='Cake',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(Cake[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cake[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='Cake',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'Cake_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(Cake[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cake[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='Cake',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(Cake[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(Cake[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='Cake',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose
;plots of diss
if test[0] ne -1 then begin
   psopen,'ds_lons',/enc,/color
   !P.MULTI=[0,2,2]
   umax=max(ds[*,0,*],min=umin)
   nlevels=25
   step=(umax-umin)/nlevels
   mylevels=indgen(nlevels)*step+umin
   contour,reform(ds[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='ds',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(ds[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ds[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='ds',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(ds[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ds[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='ds',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(ds[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ds[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='ds',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'ds_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(ds[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ds[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='ds',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(ds[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(ds[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='ds',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose
;plots of CuD
psopen,'CuD_lons',/enc,/color
!P.MULTI=[0,2,2]
umax=max(CuD[*,0,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CuD[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='CuD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(CuD[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CuD[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='CuD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(CuD[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CuD[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='CuD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(CuD[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CuD[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='CuD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'CuD_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(CuD[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CuD[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='CuD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(CuD[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CuD[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='CuD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose
;plots of CaD
psopen,'CaD_lons',/enc,/color
!P.MULTI=[0,2,2]
umax=max(CaD[*,0,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CaD[*,0,*]),reform(lat[*,0,*]),reform(vert[*,0,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=0',title='CaD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.45,0.95],format='(e8.1)'
umax=max(CaD[*,24,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CaD[*,24,*]),reform(lat[*,24,*]),reform(vert[*,24,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=90',title='CaD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.93,0.95,0.95],format='(e8.1)'
umax=max(CaD[*,48,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CaD[*,48,*]),reform(lat[*,48,*]),reform(vert[*,48,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=180',title='CaD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.45,0.53],format='(e8.1)'
umax=max(CaD[*,72,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CaD[*,72,*]),reform(lat[*,72,*]),reform(vert[*,72,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lon=-90',title='CaD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.55,0.5,0.95,0.53],format='(e8.1)'
psclose
psopen,'CaD_lats',/enc,/color;,/landscape
!P.MULTI=[0,1,2]
umax=max(CaD[12,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CaD[12,*,*]),reform(lon[12,*,*]),reform(vert[12,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=45',title='CaD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.93,0.95,0.95],format='(e8.1)'
umax=max(CaD[24,*,*],min=umin)
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(CaD[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='lat=90',title='CaD',min_value=umin,max_value=umax,/xstyle
colorbar,range=[umin,umax],position=[0.05,0.5,0.95,0.53],format='(e8.1)'
psclose
endif

!P.MULTI=0

endelse

loadct,0


;stop

end
