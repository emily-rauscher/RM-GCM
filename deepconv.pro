pro deepconv,pres,tint=tint,akap=akap,kir=kir,ga=ga

; not sure of the point of this
;
; to get the profile for deep, adiabatic layers, assuming that the
; temperature gradient from radiative diffusion is the convective
; gradient

; transition from guillot to this as:  TFAC=(1.-(1.2)**(1.-TAU))

if not keyword_set(tint) then tint=500.
fint=5.67e-8*tint^4                      ;W/m^2
if not keyword_set(akap) then akap=0.286  ;R/c_p
if not keyword_set(kir) then kir=1.e-2    ;cm^2/g
if not keyword_set(ga) then ga=8.         ;m/s^2

tau=(kir/ga)*pres*1.e4

;cc=tint*(0.75*tau)^0.25*((0.25/akap)^0.25-1.)*20.
;ftemp=(0.75*tau)^0.25*tint+cc
ftemp=tint*(tau*3./16./akap)^0.25
tg=guillot(pres,1.,/night,tint=tint)

window,0

plot,tg,pres,yr=[max(pres),min(pres)],/ylog,ystyle=8;,linestyle=1
oplot,ftemp,pres,linestyle=2

;tfac=1.-1.2^(1.-tau)
;test=where(tau ge 1.)
;temp=tg
;temp[test]=tfac[test]*ftemp[test]+(1.-tfac[test])*tg[test]

;oplot,temp,pres

convt=tint*(3.*tau*akap/16.)^0.25
oplot,convt,pres,thick=3;,linestyle=1
convtf=tint*(3.*tau/akap/16.)^0.25
oplot,convtf,pres,thick=3,linestyle=2

axis,yaxis=1,/ylog,yr=[max(tau),min(tau)],/ystyle

window,1

tgakap=(pres/tg)*deriv(pres,tg)
ftakap=(pres/ftemp)*deriv(pres,ftemp)
;tempakap=(pres/temp)*deriv(pres,temp)
tgab=deriv(alog(pres),alog(tg))  ;THIS ONE GIVES CORRECT ANSWER FOR TG
ftab=deriv(alog(pres),alog(ftemp))
;nl=n_elements(pres)
;ftac=fltarr(nl)
;ftad=fltarr(nl)
;for i=0,nl-2 do begin
;   ftac[i]=(alog(ftemp[i])-alog(ftemp[i+1]))/(alog(pres[i])-alog(pres[i+1]))
;   ftad[i]=(pres[i]/ftemp[i])*(ftemp[i]-ftemp[i+1])/(pres[i]-pres[i+1])
;endfor
;ftae=(tau/ftemp)*deriv(tau,ftemp)
;ftaf=deriv(alog(tau),alog(ftemp))

prange=fltarr(2)
prange[0]=min([tgakap,ftakap])-akap
prange[1]=0.01

plot,tgakap-akap,pres,yr=[max(pres),min(pres)],/ylog,/ystyle,xr=prange;,linestyle=1
oplot,ftakap-akap,pres,linestyle=2
oplot,tgab-akap,pres,thick=3
oplot,ftab-akap,pres,linestyle=2,thick=3
;oplot,ftae-akap,pres,linestyle=1
;oplot,ftaf-akap,pres,linestyle=1,thick=3
;oplot,tempakap-akap,pres

;tgrad=(tau/tg)*0.25*tint*(3./4)^0.25/(2./3.+tau)^0.75
;tgrad=(tau/tg)*(3./16.)*tint^4/tg^3
tgrad=tau/4./(2./3.+tau)
;oplot,tgrad-akap,pres,thick=3;,linestyle=1

tfgrad=(tint/ftemp)^4*tau*3./16.
;oplot,tfgrad-akap,pres,thick=3,linestyle=2

;return,temp

end
