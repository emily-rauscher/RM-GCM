pro erates, p0, ww, nl, ndays

; reads heating rate from fort.53, energies from fort.52, fric diss
; from fort.55
;
;ww=2.06e-5 for C&S

p0*=1.e5 ; convert from bar to Pa
period=2.*!pi/ww  ; day in seconds

window,0

openr,17,'fort.52'
xy=fltarr(2,nl,ndays+1)         ;0: kinetic, 1: kinetic + cpT
readf,17,xy
close,17
xy*=p0

tke=total(reform(xy[0,*,*]),1)
tte=total(reform(xy[1,*,*]),1)

openr,18,'fort.53'
h=fltarr(3,ndays+1)
readf,18,h
close,18
h*=p0

ih=reform(h[2,*])*period

;what if ih negative?

openr,19,'fort.55'
fdiss=fltarr(ndays+1)
readf,19,fdiss
close,19
fdiss*=p0

ifdiss=fdiss*period


nonz=where([tke,tte,ih] ne 0.)
ymin=min(([tke,tte,ih])[nonz],max=ymax)

print,minmax(tte),minmax(tke),minmax(ih)

plot,tke,xtitle='day',ytitle='Energy (J)',yr=[ymin,ymax],$
     charsize=1.5,/xstyle,linestyle=1,/ylog
oplot,tte
oplot,tte-tke,linestyle=2
oplot,ih,linestyle=3
oplot,ifdiss,linestyle=4

;numerical dissipation effective friction timescale (in sec) by tke/h

;d(tte)/dt by (tte[i+1]-tte[i])/(day in sec)
; day in sec = 2pi/ww
dedt=fltarr(ndays)
for i=0,ndays-1 do dedt[i]=(tte[i+1]-tte[i])/period

stop

end
