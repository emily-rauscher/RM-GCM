pro fluxes,nlat,nlon,nlev,mbar=mbar,fbase=fbase,psfile=psfile,solc=solc,noheader=noheader

; plots the total LW (and SW) fluxes at each flux level (or maybe heating rate
; at full levels?), as well as what radiative transfer alone would have given
; reads fort.63
;
; only plots cirrad fluxes if fdiff keyword set
;
; pressures in bar, unless mbar keyword used (Jun 2011 change to
; fort.63 output format)

; format of fort.63 is:
; LATITUDE, LONGITUDE:  -81.6505907502981       0.000000000000000E+000
;           UPWARD FLUX (Wm-2)     DOWNWARD FLUX (Wm-2)     NET FLUX (Wm-2)           SW NET FLUX (Wm-2)
;         0.616      0.21033E+06           0.00000E+00           0.21033E+06           0.15074E+06
;       .... (nlev+1 entries total, last being surface) ........
;    101645.531      0.33613E+07           0.33578E+07           0.35000E+04           0.00000E+00

entry = create_struct( 'lat', 0.0, $
                       'lon', 0.0, $
                       'day', 0B, $
                       'pres', fltarr(nlev), $
                       'sflux', fltarr(nlev), $  ; SW flux
                       'flux', fltarr(nlev), $   ; actual IR flux
                       'cflux', fltarr(nlev), $   ; upward (cooling) flux
                       'rflux', fltarr(nlev))    ; cirrad flux
nent=(nlat-1)*nlon
rates=replicate(entry,nent)

read1=strarr(2)
data=fltarr(5,nlev+1)
read4=''
toss=fltarr(3*(nlev+1))
read5=''

dflux=fltarr(nlev)
nflux=fltarr(nlev)
dcflux=fltarr(nlev)
ncflux=fltarr(nlev)
dsflux=fltarr(nlev)

openr,17,'fort.63'
if not keyword_set(noheader) then begin
   header=strarr(2)
   readf,17,header
endif
for i=0,nent-1 do begin
   readf,17,read1
;   print,read1[0]
   rates[i].lat=float(strmid(read1[0],23,6))
;   print,rates[i].lat
   rates[i].lon=float(strmid(read1[0],47,6))
;   print,rates[i].lon
   readf,17,data
   rates[i].pres=reform(data[0,0:nlev-1])
   rates[i].rflux=reform(data[1,0:nlev-1]-data[2,0:nlev-1])
   rates[i].flux=reform(data[3,0:nlev-1])
   rates[i].sflux=reform(data[4,0:nlev-1])
   rates[i].cflux=rates[i].flux-rates[i].sflux
   readf,17,read4
;   print,read4
   readf,17,toss
   readf,17,read5
;   print,read5
   if ((rates[i].lon le 90.) or (rates[i].lon ge 270.)) then begin
      rates[i].day=1
      dflux+=(rates[i].flux*cos(rates[i].lat*!pi/180.))
      dcflux+=(rates[i].cflux*cos(rates[i].lat*!pi/180.))
      dsflux+=(rates[i].sflux*cos(rates[i].lat*!pi/180.))
   endif else begin
      nflux+=(rates[i].flux*cos(rates[i].lat*!pi/180.))
      ncflux+=(rates[i].cflux*cos(rates[i].lat*!pi/180.))
   endelse
endfor
close,17

dflux*=(!pi/2./nlon/nlat)
nflux*=(!pi/2./nlon/nlat)
dcflux*=(!pi/2./nlon/nlat)
ncflux*=(!pi/2./nlon/nlat)
dsflux*=(!pi/2./nlon/nlat)

if keyword_set(mbar) then begin
   rates.pres/=1.e3             ;convert to bar
endif

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=40,ysize=15;,bits_per_pixel=24
   !p.thick=5
   !x.thick=5
   !y.thick=5
   !p.font=0
endif else begin
  device,true_color=24,decomposed=0,retain=2
endelse
  !p.charsize=3

!p.multi=[0,3,1,1,0]

umin=min(rates.flux,max=umax)
pmin=min(rates.pres,max=pmax)
plot,rates[0].flux,rates[0].pres,yr=[pmax,pmin],/ylog,/ystyle,xr=[umin,umax],$
     xtit='Net LW flux [W/m^2]',ytit='Pressure [bar]',xmargin=[12,3];,/xlog
loadct,5,/silent
for i=0,nent-1 do begin
   if rates[i].day eq 0 then oplot,rates[i].flux,rates[i].pres,color=50 else $
      oplot,rates[i].flux,rates[i].pres,color=100
endfor
loadct,0,/silent
oplot,dflux,rates[0].pres,thick=8
oplot,nflux,rates[0].pres,thick=8
if keyword_set(fbase) then begin
   oplot,[fbase,fbase],[pmax,pmin],linestyle=2
endif

umin=min(rates.sflux,max=umax)
plot,-rates[0].sflux,rates[0].pres,yr=[pmax,pmin],/ylog,/ystyle,xr=[-umax,-umin],$
     xtit='Net SW flux [W/m^2]',ytit='Pressure [bar]',xmargin=[12,3];,/xlog
loadct,5,/silent
for i=0,nent-1 do begin
   if rates[i].day eq 0 then oplot,-rates[i].sflux,rates[i].pres,color=50 else $
      oplot,-rates[i].sflux,rates[i].pres,color=100
endfor
loadct,0,/silent
oplot,-dsflux,rates[0].pres,thick=8
if keyword_set(solc) then begin
   oplot,[-solc/2.,-solc/2.],[pmax,pmin],linestyle=1
endif

umin=min(rates.cflux,max=umax)
plot,rates[0].cflux,rates[0].pres,yr=[pmax,pmin],/ylog,/ystyle,xr=[umin,umax],$
     xtit='Net total flux [W/m^2]',ytit='Pressure [bar]',xmargin=[12,3];,/xlog
loadct,5,/silent
for i=0,nent-1 do begin
   if rates[i].day eq 0 then oplot,rates[i].cflux,rates[i].pres,color=50 else $
      oplot,rates[i].cflux,rates[i].pres,color=100
endfor
loadct,0,/silent
oplot,dcflux,rates[0].pres,thick=8
oplot,ncflux,rates[0].pres,thick=8
if keyword_set(fbase) then begin
   oplot,[fbase,fbase],[pmax,pmin],linestyle=2
endif
if keyword_set(solc) then begin
   oplot,[-solc/4.,-solc/4.],[pmax,pmin],linestyle=1
endif
topflux=(dcflux[0]+ncflux[0])/2.
oplot,[topflux,topflux],[pmin,pmin],psym=7,thick=5,symsize=2

print
print,'Dayside SW flux entering top boundary:',dsflux[0],' W/m^2'
print,'Total SW flux entering top boundary:',dsflux[0]/2.,' W/m^2'
print
print,'Dayside LW flux leaving top boundary:',dflux[0],' W/m^2'
print,'Nightside LW flux leaving top boundary:',nflux[0],' W/m^2'
print,'Total LW flux leaving top boundary:',(dflux[0]+nflux[0])/2.,' W/m^2'
print
print,'Dayside flux leaving top boundary:',dcflux[0],' W/m^2'
print,'Nightside flux leaving top boundary:',ncflux[0],' W/m^2'
print,'Total flux leaving top boundary:',topflux,' W/m^2'

!p.multi=0
if keyword_set(psfile) then begin
   !p.thick=1
   !x.thick=1
   !y.thick=1
   psclose
endif

;stop

end
