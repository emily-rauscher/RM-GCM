pro brsm_plot

psopen,'test',/enc,x=30,y=25,/color

!p.font=0

night=[670.,678.,687.,698.,711.,725.,742.,760.,782.,806.,$
833.,864.,899.,938.,983.,1036.,1094.,1159.,1240.,1345.,$
1428.,1488.,1531.,1546.,1524.,1453.,1401.,1366.,1341.,1321.,$
1307.,1295.,1287.,1280.,1275.,1271.,1269.,1269.,1269.,1270.,$
1270.,1272.,1274.,1277.,1283.]

day0=night+[402.,401.,401.,401.,402.,403.,406.,409.,413.,418.,$
423.,426.,427.,428.,426.,422.,415.,399.,361.,267.,$
228.,241.,294.,393.,539.,723.,855.,941.,999.,1042.,$
1071.,1092.,1106.,1115.,1121.,1125.,1126.,1127.,1127.,1126.,$
1126.,1125.,1123.,1120.,1115.]

day3=night+[1085.,1088.,1094.,1101.,1110.,1120.,1131.,1142.,1152.,1159.,$
1162.,1157.,1140.,1109.,1059.,972.,801.,606.,480.,316.,$
194.,121.,97.,147.,283.,480.,623.,720.,788.,839.,$
875.,899.,915.,926.,933.,936.,938.,938.,938.,938.,$
937.,936.,934.,931.,926.]

trad=[0.0023,0.0025,0.0028,0.0032,0.0036,0.0043,0.0051,0.0060,0.0074,0.0085,$
0.0099,0.011,0.013,0.015,0.017,0.019,0.022,0.025,0.031,0.039,$
0.048,0.059,0.065,0.072,0.079,0.096,0.13,0.18,0.25,0.34,$
0.46,0.59,0.89,2.2,5.5,13.,26.,57.,520.,2100.,$
8200.,33000.,82000.,210000.,820000.]*2.*!pi/2.06e-5

print,minmax(trad)

p=100.*get_sigma(7,45)

plot,p,night,/xlog,yr=[500,2500],ystyle=8,linestyle=2,xmargin=[8,8],$
     xtit='Pressure [bar]',ytit='Temperature [K]',$
     thick=7,ythick=7,xthick=7,charsize=2,xr=[1.e-5,1.e2],/xstyle,$
     xtickname=['','10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','10!U0!N','10!U1!N']
oplot,p,day0,linestyle=1,thick=7
oplot,p,day3,linestyle=3,thick=7
xyouts,5,1150,'night',charsize=2
xyouts,2.e-5,2100,'day, w/absorber',charsize=2
xyouts,2.e-5,1400,'day, w/o',charsize=2
xyouts,2.e-5,1300,'absorber',charsize=2
loadct,5
plot,p,trad,/noerase,ystyle=4,/ylog,yrange=[100,1.e12],xr=[1.e-5,1.e2],$
     thick=7,/xlog,xmargin=[8,8],xthick=7,ythick=7,charsize=2,/xstyle,$
     xtickname=['','10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','10!U0!N','10!U1!N']
oplot,p,trad,thick=7,color=100
axis,yaxis=1,/ylog,yrange=[100,1.e12],ythick=7,$
     ytit='Relaxation time [s]',charsize=2
xyouts,1.e-2,2.e3,'relaxation time',color=100,charsize=2
loadct,0
;LABELS, COLOR(S)

psclose

;plot,night,p,/ylog,yr=[100.,1.e-5],/ystyle,xr=[600,2500],$
;     ytit='Pressure [bar]',xtit='Temperature [K]',thick=5,xthick=5,ythick=5,$
;     charsize=2
;oplot,day0,p,linestyle=1,thick=5
;oplot,day3,p,linestyle=2,thick=5

end
