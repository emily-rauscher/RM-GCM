pro energy_diag, p0, nl, ndays, plotlev=plotlev, all=all

; plots kinetic and total energy at any level, throughout the run
; plotlev selects level to plot, otherwise can use all to plot all
; p0 (pressure at z=0) needs to be set to get mass normalization correct
; (input as bar, will convert to Pa)

  ;by default will plot surface values (nl-1), otherwise set with plotlev
  if not keyword_set(all) then begin
     if not keyword_set(plotlev) then plotlev=nl-1
  endif else plotlev=findgen(nl)

  window,0

  openr,17,'fort.52'
    xy=fltarr(2,nl,ndays+1)   ;0: kinetic, 1: kinetic + cpT
  readf,17,xy
  close,17

  xy*=p0*1.e5

  nonz=where(xy(*,plotlev,*) ne 0.)
  ymin=min((xy(*,plotlev,*))[nonz],max=ymax)
  print,ymin,ymax
  nonz=where(xy(0,plotlev,*) ne 0.)
  print,min((xy[0,plotlev,*])[nonz]),max((xy[0,plotlev,*])[nonz])

  if (size(plotlev))[0] eq 0 then begin
     plot,xy(0,plotlev,*),xtitle='day',ytitle='Energy (J)',yr=[ymin,ymax],$
          linestyle=1,charsize=1.5,/xstyle,/ylog
     oplot,xy(1,plotlev,*)
     print,'L=',plotlev+1
  endif else begin
     plot,findgen(ndays+1),replicate(0.,ndays+1),$
          xtitle='day',ytitle='Energy (J)',yr=[ymin,ymax],charsize=1.5,/xstyle,/ylog
     for i=0,nl-1 do begin
        oplot,xy(0,i,*),linestyle=1
        oplot,xy(1,i,*)
     endfor
  endelse

end
