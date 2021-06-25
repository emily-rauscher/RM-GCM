pro diagnostics, nl, ndays, plotlev=plotlev, dayrange=dayrange,uz=uz,p0=p0,oom=oom,psfile=psfile

; reads wind diagnostics from fort.27
; may also read temperature diagnostics from fort.28 (unless keyword
;    uz set)
; requires: get_sigma.pro, psopen.pro, psclose.pro, colorbar.pro
;
; nl: number of vertical levels
; ndays: length of simulation
; plotlev: level to plot (0:surface is the default)
; dayrange: range of days to plot, 0:ndays is the default
; uz: if set, will plot the zonal averaged u at the equator as
;     function of time, for each level, and no temperature data
; p0: surface pressure in bar, default: 100
; oom: orders of magnitude in pressure domain, default:5
; psfile: if set will make eps file instead of plotting in a window
;         (can set as string for filename, otherwise file is plot.eps)

if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=30,ysize=20;,bits_per_pixel=24
  !p.font=0
  !p.charsize=1.75
endelse

;by default will plot surface values (nl-1), otherwise set with plotlev
if not keyword_set(plotlev) then plotlev=0
if not keyword_set(dayrange) then dayrange=[0.,ndays]

openr,17,'fort.27'
xy=fltarr(3,nl,ndays+1)         ;0: max u, 1: min u, 2: [u]
readf,17,xy
close,17

if keyword_set(uz) then begin

   if not keyword_set(p0) then p0=1.e2
   if not keyword_set(oom) then oom=5
   days=findgen(dayrange[1]-dayrange[0]+1)+dayrange[0]
   levs=p0*get_sigma(oom,nl)
   uzn=min(xy[2,*,*],max=uzx)
   range=uzx-uzn
   nlevels=45
   step=(uzx-uzn)/(nlevels-1)
   mylevels=findgen(nlevels)*step+uzn
   loadct,10,/silent
   contour,transpose(reform(xy[2,*,dayrange[0]:dayrange[1]])),days,levs,/cell_fill,$
           levels=mylevels,yr=[max(levs),min(levs)],min_value=uzn,max_value=uzx,/ylog,$
           /xstyle,/ystyle,xmargin=[9,8],ymargin=[4,2],$
           xtitle='Planet Day',ytitle='Pressure [bar]';,title='[u] at equator'
   contour,transpose(reform(xy[2,*,dayrange[0]:dayrange[1]])),days,levs,levels=[0],$
           color=255,/over
   colorbar,range=[uzn,uzx],position=[0.91,0.17,0.94,0.87],$
            format='(i6)',/vertical,/right,charsize=1.25,$
            divisions=4
   loadct,0,/silent
endif else begin

   openr,18,'fort.28'
   ab=fltarr(6,nl,ndays+1)
   readf,18,ab
   close,18
;  window,0
  ymax=max(xy(*,plotlev,*))
  ymin=min(xy(*,plotlev,*))
  plot,xy(0,plotlev,dayrange[0]:dayrange[1]),xtitle='days',ytitle='u (m/s)',yr=[ymin,ymax],$
       linestyle=1,charsize=1.5,/xstyle
  oplot,xy(1,plotlev,dayrange[0]:dayrange[1]),linestyle=2
  oplot,xy(2,plotlev,dayrange[0]:dayrange[1])

;  window,2
  ymax=max(ab(*,plotlev,*))
  ymin=min(ab(*,plotlev,*))
  plot,ab(0,plotlev,dayrange[0]:dayrange[1]),xtitle='days',ytitle='T (K)',$
       yr=[ymin,ymax],charsize=1.5,/xstyle
  oplot,ab(1,plotlev,dayrange[0]:dayrange[1])
  oplot,ab(2,plotlev,dayrange[0]:dayrange[1]),linestyle=1
  oplot,ab(3,plotlev,dayrange[0]:dayrange[1]),linestyle=1
  oplot,ab(4,plotlev,dayrange[0]:dayrange[1]),linestyle=2
  oplot,ab(5,plotlev,dayrange[0]:dayrange[1]),linestyle=2

  print,strcompress('L'+string(plotlev+1),/remove_all)

  print,'Relative variation of [u]',$
        stddev(xy[2,plotlev,dayrange[0]:dayrange[1]])/mean(xy[2,plotlev,dayrange[0]:dayrange[1]])
  print,'Relative change of [u]',$
        (xy[2,plotlev,dayrange[0]]-xy[2,plotlev,dayrange[1]])/mean(xy[2,plotlev,dayrange[0]:dayrange[1]])
  print,'Relative variation of T (north pole)',$
        stddev(ab[4,plotlev,dayrange[0]:dayrange[1]])/mean(ab[4,plotlev,dayrange[0]:dayrange[1]])
  print,'Relative change:',$
        (ab[4,plotlev,dayrange[0]]-ab[4,plotlev,dayrange[1]])/mean(ab[4,plotlev,dayrange[0]:dayrange[1]])
  print,'Relative variation of T (east limb)',$
        stddev(ab[2,plotlev,dayrange[0]:dayrange[1]])/mean(ab[2,plotlev,dayrange[0]:dayrange[1]])
  print,'Relative change:',$
        (ab[2,plotlev,dayrange[0]]-ab[2,plotlev,dayrange[1]])/mean(ab[2,plotlev,dayrange[0]:dayrange[1]])
  print,'Relative variation of T (substellar)',$
        stddev(ab[0,plotlev,dayrange[0]:dayrange[1]])/mean(ab[0,plotlev,dayrange[0]:dayrange[1]])
  print,'Relative change:',$
        (ab[0,plotlev,dayrange[0]]-ab[0,plotlev,dayrange[1]])/mean(ab[0,plotlev,dayrange[0]:dayrange[1]])
endelse

if keyword_set(psfile) then begin
   !p.font=-1
   !p.charsize=1
   psclose
   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
   spawn,'convert plot.png eps3:plot.eps'
endif

;stop

end
