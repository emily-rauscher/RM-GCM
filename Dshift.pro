pro Dshift,radea,ww,psfile=psfile,urange=urange,average=average

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

    lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
    u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])

;terme=where(lon[0,*,0] eq 90.)
;termw=where(lon[0,*,0] eq 270.)
terme=nlon*0.25
termw=nlon*0.75
if keyword_set(average) then begin
   terme=[terme-1,terme,terme+1]
   termw=[termw-1,termw,termw+1]
endif

  OpenR, 18, 'vertical.txt'
  ReadF, 18, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  ab=fltarr(3,nlat,nlon,nlev)   
  ReadF, 18, ab
  Close, 18

  sig=reform(ab[0,*,*,*])
  pres=reform(ab[1,*,*,*])
  z=reform(ab[2,*,*,*])

  vels=fltarr(4,nlat*nlev*2)    ; 0: r, 1: theta, 2: losv, 3: w/rot
  count=0

if keyword_set(average) then begin
   for i=0,nlat-1 do begin
      for j=0,nlev-1 do begin
         vels[0,count]= mean(z[i,terme,j])+radea
         vels[0,count+1]= mean(z[i,termw,j])+radea
         vels[1,count]= (!pi/180.)*(180.-mean(lat[i,terme,j]))
         vels[1,count+1]= (!pi/180.)*mean(lat[i,termw,j])
         vels[2,count]= mean(u[i,terme,j]*sin((!pi/180.)*lon[i,terme,j])*cos((!pi/180.)*lat[i,terme,j]) $
                            +v[i,terme,j]*cos((!pi/180.)*lon[i,terme,j])*sin((!pi/180.)*lat[i,terme,j]))
         vels[2,count+1]= mean(u[i,termw,j]*sin((!pi/180.)*lon[i,termw,j])*cos((!pi/180.)*lat[i,termw,j]) $
                              +v[i,termw,j]*cos((!pi/180.)*lon[i,termw,j])*sin((!pi/180.)*lat[i,termw,j]))
         vels[3,count]=vels[2,count]+ ww*vels[0,count] $
                       *mean(cos((!pi/180.)*lat[i,terme,j])*sin((!pi/180.)*lon[i,terme,j]))
         vels[3,count+1]=vels[2,count+1]+ ww*vels[0,count+1] $
                       *mean(cos((!pi/180.)*lat[i,termw,j])*sin((!pi/180.)*lon[i,termw,j]))
         vels[0,count:count+1]=vels[0,count:count+1]/radea
         count=count+2
      endfor
   endfor
endif else begin
   for i=0,nlat-1 do begin
      for j=0,nlev-1 do begin
         vels[0,count]= z[i,terme,j]+radea
         vels[0,count+1]= z[i,termw,j]+radea
         vels[1,count]= (!pi/180.)*(180.-lat[i,terme,j])
         vels[1,count+1]= (!pi/180.)*lat[i,termw,j]
         vels[2,count]= cos((!pi/180.)*lat[i,terme,j])*u[i,terme,j]
         vels[2,count+1]= -1.*cos((!pi/180.)*lat[i,termw,j])*u[i,termw,j]
         vels[3,count]=vels[2,count]+ ww*vels[0,count]*cos((!pi/180.)*lat[i,terme,j])
         vels[3,count+1]=vels[2,count+1]- ww*vels[0,count+1]*cos((!pi/180.)*lat[i,termw,j])
         vels[0,count:count+1]=vels[0,count:count+1]/radea
         count=count+2
      endfor
   endfor
endelse

;switch velocities so that negative it towards us (blue-shifted):
vels[2:3,*]*=-1.

if not keyword_set(urange) then begin
   umax=max(vels[2:3,*],min=umin) 
endif else begin
   umin=urange[0]
   umax=urange[1]
endelse

absmax=max([abs(umax),abs(umin)])
range=umax-umin
maxrange=2.*absmax  ;from -absmax to +absmax
rat=range/maxrange
cbottom=(absmax+umin)/maxrange*255.
nlevels=44.*rat+1.  ;so that delta color always = 255./44.
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
edges=max(alog10(vels[0,*]),min=bcirc)

plevels=[10.,1.,0.1,1.e-2,1.e-3,1.e-4]
pstring=['10 bar','1 bar','0.1 bar','10 mbar','1 mbar','0.1 mbar']
if keyword_set(average) then begin
   tinde=terme[1]
   tindw=termw[1]
endif else begin
   tinde=terme
   tindw=termw
endelse
temppe=reform(pres[*,tinde,*])
temple=(!pi/180.)*(180.-reform(lat[*,tinde,*]))
tempze=reform(z[*,tinde,*])
temppw=reform(pres[*,tindw,*])
templw=(!pi/180.)*reform(lat[*,tindw,*])
tempzw=reform(z[*,tindw,*])


if keyword_set(psfile) then begin
   psopen,psfile,/enc,/color,x=15,y=30
   tpp=[0.1,0.45,0.9,0.85]
   bpp=[0.1,0.05,0.9,0.45]
   cbp=[0.15,0.87,0.85,0.9]
   tptit=''
   bptit=''
endif else begin
  device,true_color=24,decomposed=0,retain=2
  window, 0, xsize=800, ysize=800, retain=2
   tpp=[0.1,0.1,0.9,0.9]
   bpp=[0.1,0.1,0.9,0.9]
   cbp=[0.15,0.91,0.85,0.95]
   tptit='in the frame of the planet'
   bptit='including the effect of rotation'
endelse
  loadct,0,/silent

polar_contour,reform(vels[2,*]),reform(vels[1,*]),alog10(reform(vels[0,*])),/cell_fill,min_value=umin,max_value=umax,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels,/closed,/irregular,/isotropic,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=tpp,subtitle=tptit
loadct,25,/silent
polar_contour,reform(vels[2,*]),reform(vels[1,*]),alog10(reform(vels[0,*])),/cell_fill,min_value=umin,max_value=umax,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels,/closed,/irregular,/isotropic,/overplot,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=tpp
colorbar,range=[umin,umax],position=cbp,format='(e8.1)',bottom=cbottom,ncolors=256*rat,/top,divisions=4,charsize=1.15
loadct,0,/silent
tvcircle,bcirc,0.,0.,/data,/fill,color='black'

for i=0,n_elements(plevels)-1 do begin
   near=min(abs(temppe-plevels[i]),index)
   near=where((temppe gt 0.9*temppe(index)) and (temppe lt 1.1*temppe(index)))
   plot,/polar,alog10((tempze(near)+radea)/radea),temple[near],/noerase,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=tpp;,psym=1
   near=min(abs(temppw-plevels[i]),index)
   near=where((temppw gt 0.9*temppw(index)) and (temppw lt 1.1*temppw(index)))
   plot,/polar,alog10((tempzw(near)+radea)/radea),templw[near],/noerase,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=tpp;,psym=1
   xyouts,0.,0.-max(alog10((tempzw(near)+radea)/radea)),pstring[i],align=0.5
endfor

if not keyword_set(psfile) then begin
  window, 2, xsize=800, ysize=800, retain=2
  polar_contour,reform(vels[3,*]),reform(vels[1,*]),alog10(reform(vels[0,*])),/cell_fill,min_value=umin,max_value=umax,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels,/closed,/irregular,/isotropic,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=bpp,subtitle=bptit
endif else begin
   polar_contour,reform(vels[3,*]),reform(vels[1,*]),alog10(reform(vels[0,*])),/cell_fill,min_value=umin,max_value=umax,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels,/closed,/irregular,/isotropic,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=bpp,subtitle=bptit,/noerase
endelse
loadct,25,/silent
polar_contour,reform(vels[3,*]),reform(vels[1,*]),alog10(reform(vels[0,*])),/cell_fill,min_value=umin,max_value=umax,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels,/closed,/irregular,/isotropic,/overplot,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=bpp
colorbar,range=[umin,umax],position=cbp,format='(e8.1)',bottom=cbottom,ncolors=256*rat,/top,divisions=4,charsize=1.15
loadct,0,/silent
tvcircle,bcirc,0.,0.,/data,/fill,color='black'

for i=0,n_elements(plevels)-1 do begin
   near=min(abs(temppe-plevels[i]),index)
   near=where((temppe gt 0.9*temppe(index)) and (temppe lt 1.1*temppe(index)))
   plot,/polar,alog10((tempze(near)+radea)/radea),temple[near],/noerase,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=bpp;,psym=1
   near=min(abs(temppw-plevels[i]),index)
   near=where((temppw gt 0.9*temppw(index)) and (temppw lt 1.1*temppw(index)))
   plot,/polar,alog10((tempzw(near)+radea)/radea),templw[near],/noerase,xr=[-1.*edges,edges],yr=[-1.*edges,edges],/xstyle,/ystyle,position=bpp;,psym=1
   xyouts,0.,0.-max(alog10((tempzw(near)+radea)/radea)),pstring[i],align=0.5
endfor

if keyword_set(psfile) then psclose

loadct,0,/silent

if keyword_set(average) then begin
   elx=max(lon[*,terme,*],min=eln)
   wlx=max(lon[*,termw,*],min=wln)
   print,'averaged over: ',eln,elx,' and ',wln,wlx
endif else begin
   print,'at: ',lon[0,terme,0],lon[0,termw,0]
endelse

;stop

end
