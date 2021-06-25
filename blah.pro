pro blah

mg=96.
jg=24.
oom=5.34363
p0=220.16*1.e5
radea=9.44e7
ga=9.42
gascon=4593.
akap=0.321
cp=gascon/akap

openr,17,'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)
readf,17,xy
close,17

lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*])
u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])

openr,18,'fort.50'
readf,18,nlat,nlon
ab=fltarr(3,nlat,nlon)
readf,18,ab
close,18

sp=reform(ab[2,*,*])
sp=(sp+1.)*p0

sigma=fltarr(nlev)
stp=-1.*oom/nlev
sigma[nlev-1]=10.^(stp/2.)
for i=nlev-2,0,-1 do sigma[i]=sigma[i+1]*10.^(stp)
dsigma=fltarr(nlev)
dsigma[0]=((sigma[0]+sigma[1])/2.-0.)
dsigma[nlev-1]=(1.-(sigma[nlev-2]+sigma[nlev-1])/2.)
for i=1,nlev-2 do dsigma[i]=(sigma[i+1]+sigma[i])/2.-(sigma[i]+sigma[i-1])/2.

ke=fltarr(nlat,nlon,nlev)
lke=fltarr(nlev)
te=fltarr(nlat,nlon,nlev)
lte=fltarr(nlev)

for l=0,nlev-1 do begin
   dmassl=radea^2*(2.*!pi/mg)*(!pi/2./jg)/ga*dsigma[l]
   for j=0,nlat-1 do begin
      dmass=dmassl*cos(lat[j,0.,l]*!pi/180.)
      ke[j,*,l]=dmass*sp[j,*]*0.5*(u[j,*,l]^2+v[j,*,l]^2)
      te[j,*,l]=ke+dmass*sp[j,*]*cp*temp[j,*,l]
   endfor
   lke[l]=total(ke[*,*,l])
   lte[l]=total(te[*,*,l])
endfor

plot,lte,sigma,yr=[max(sigma),min(sigma)],/ylog,/xlog
oplot,lke,sigma

stop

end
