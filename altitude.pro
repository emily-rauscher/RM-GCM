pro altitude,R,g,TGR,p0=p0,oom=oom,base=base,z=z

;purpose: to calculate altitudes at each point in 3D model space
;         (sigma, lat, lon) and print file with: sigma, pressure, altitude
;
;reads fort.26 for T,u,v; fort.50 for p_surf

;if not keyword_set(R) then R=3500.
;if not keyword_set(g) then g=10.
;TGR is the "base" temperature for p=p0
;
;if not integrating from p0,TGR, can set to integrate from base=[pressure,temperature]

if not keyword_set(p0) then p0=1.

;get surface presssures
openr,17,'fort.50'
readf,17,nlat,nlon
xy=fltarr(3,nlat,nlon)
readf,17,xy
close,17

lon=reform(xy[0,*,*]) & lat=reform(xy[1,*,*]) & sp=reform(xy[2,*,*])
sp=(sp+1.)*p0

;get T
openr,18,'fort.26'
readf,18,nlat2,nlon2,nlev
if (nlat2 ne nlat) or (nlon2 ne nlon) then print,'ERROR: dimensions of fort.26 do not match fort.50'
ab=fltarr(6,nlat,nlon,nlev)
readf,18,ab
close,18

temp=reform(ab[5,*,*,*])

;set sigma
sigma=fltarr(nlev)
if keyword_set(oom) then begin
   stp=-1.*oom/nlev
   sigma[nlev-1]=10.^(stp/2.)
   for i=nlev-2,0,-1 do sigma[i]=sigma[i+1]*10.^(stp)
endif else begin
   sigma=(findgen(nlev)+0.5)/nlev
endelse


if not keyword_set(base) then begin
;create array to hold heights (z=0 at p=p0)
   z=fltarr(nlat,nlon,nlev)

;set altitude of first level (up from base=p0, where T=TGR)
   z[*,*,nlev-1]=(R/g)*0.5*(temp[*,*,nlev-1]+TGR)*alog(p0/sp[*,*]/sigma[nlev-1])

;integrate hydrostatic to solve for higher levels
   for i=nlev-2,0,-1 do begin
      z[*,*,i]=z[*,*,i+1]+(R/g)*0.5*(temp[*,*,i]+temp[*,*,i+1])*alog(sigma[i+1]/sigma[i])
   endfor

   openw,19,'vertical.txt'
   printf,19,nlat,nlon,nlev
;printf,19,z
   for i=0,nlev-1 do begin
      for j=0,nlon-1 do begin
         for k=0,nlat-1 do begin
            printf,19,sigma[i],sigma[i]*sp[k,j],z[k,j,i]
         endfor
      endfor
   endfor
   if not keyword_set(oom) then oom=0
   printf,19
   printf,19,'R,g,TGR,p0,OOM:',R,g,TGR,p0,oom
   close,19
endif else begin

;find first pressure/sigma level up from base
   test=where(sigma*p0 le base[0])
   if test[0] eq -1 then begin
      print,"input error"
      stop
   endif else begin
      bnlev=n_elements(test)
      z=fltarr(nlat,nlon,bnlev)
      z[*,*,bnlev-1]=(R/g)*0.5*(temp[*,*,test[bnlev-1]]+base[1])*alog(base[0]/sp[*,*]/sigma[test[bnlev-1]])
      for i=bnlev-2,0,-1 do begin
         z[*,*,i]=z[*,*,i+1]+(R/g)*0.5*(temp[*,*,test[i]]+temp[*,*,test[i+1]])*alog(sigma[test[i+1]]/sigma[test[i]])
      endfor

      openw,19,'vertical.txt'
      printf,19,nlat,nlon,bnlev
      for i=0,bnlev-1 do begin
         for j=0,nlon-1 do begin
            for k=0,nlat-1 do begin
               printf,19,sigma[i],sigma[i]*sp[k,j],z[k,j,i]
            endfor
         endfor
      endfor
      if not keyword_set(oom) then oom=0
      printf,19
      printf,19,'R,g,p0,OOM,base:',R,g,p0,oom,base
      close,19
   endelse

endelse

end
