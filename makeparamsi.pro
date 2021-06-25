PRO makeparamsi,Tnum=Tnum,Lnum=Lnum

;PROGRAM TO GENERATE ABBREVIATED PARAMS.I FILE FOR EASE OF QUICKLY CHANGING RESOLUTIONS;
;This was written to assist with optimization ~ mtroman 8/7/15

;;;;;;;;;INPUTS: 
;;;;;;;;; Tnum -- The horizontal resolution parameter as specified in params.i
;;;;;;;;; L -- The number of vertical layers
;;;;;;;;; Note that all other parameters were taken from params.i file as written


;CODE

	  randomfnum=3
	  filename='params.i'
	  ;pick some random file number, close it and open a new one
	  filnum = randomfnum
      close, filnum
      ;pick some 
      openw, filnum, filename
      print,'opened file:',filnum

 if (N_params() LT 1)then begin
     print,'Syntax - ' + $
             'makeparamsi, T-resolution", # vertical layers (in quotes)'
     goto, YOUPICK
     endif 
     
     if (N_params() EQ 1) then begin
     print,'Syntax - ' + $
            'makeparamsi, T-resolution, # vertical layers'
     printf,filnum,'Assuming L30'
     T=strtrim(string(T),2)
     L='L30'
     endif 
     
    if (N_params() EQ 2)then begin
	T=strtrim(string(Tnum),2)
	L=strtrim(string(L),2)
	endif
    if (N_params() GT 2)then goto, QUIT
    
	  
    
    
If Tnum eq 'T2'  then begin

   if L eq 'L20' then begin
   printf,filnum,'C T2 L20  full sphere'
   printf,filnum,'   PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=2 0,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
   printf,filnum,'    +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
   goto,quit
   endif
   
   if L eq 'L30' then begin
printf,filnum,'C T2 L30  full sphere'
printf,filnum,'     PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=30,MOCT=1,MG=8,JG=2 ,NWJ2=2'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	
    if L eq 'L40' then begin
printf,filnum,'C T2 L40  full sphere'
printf,filnum,'     PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=40,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit 
	endif
	
   if L eq 'L50' then begin
printf,filnum,'C T2 L50  full sphere'
printf,filnum,'     PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=50,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit 
	endif
	
	   if L eq 'L100' then begin
printf,filnum,'C T2 L100  full sphere'
printf,filnum,'     PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=100,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	
	   if L eq 'L200' then begin
printf,filnum,'C T2 L200  full sphere'
printf,filnum,'     PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=2 00,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	
endif


If Tnum eq 'T5' then begin

   if L eq 'L20' then begin
printf,filnum,'C T5 L20  full sphere'
printf,filnum,'     PARAMETER(NN=5,MM=5,NHEM=2 ,NL=2 0,MOCT=1,MG=16,JG=4,NWJ2=9        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
   endif
endif
   
If Tnum eq 'T10' then begin
	if L eq 'L5' then begin
printf,filnum,'C T10 L5  full sphere'
printf,filnum,'     PARAMETER(NN=10,MM=10,NHEM=2 ,NL=5,MOCT=1,MG=32,JG=8,NWJ2=30'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit 
	endif

   if L eq 'L10' then begin
printf,filnum,'C T10 L10  full sphere'
printf,filnum,'     PARAMETER(NN=10,MM=10,NHEM=2 ,NL=10,MOCT=1,MG=32,JG=8,NWJ2=30'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit 
	endif
	if L eq 'L20' then begin
printf,filnum,'C T10 L20  full sphere'
printf,filnum,'     PARAMETER(NN=10,MM=10,NHEM=2 ,NL=2 0,MOCT=1,MG=32,JG=8,NWJ2=30'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit 
	endif
	if L eq 'L30' then begin
printf,filnum,'C T10 L30  full sphere'
printf,filnum,'     PARAMETER(NN=10,MM=10,NHEM=2 ,NL=30,MOCT=1,MG=32,JG=8,NWJ2=30'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'
goto,quit  
	endif
	if L eq 'L36' then begin
printf,filnum,'C T10 L36  full sphere'
printf,filnum,'     PARAMETER(NN=10,MM=10,NHEM=2 ,NL=36,MOCT=1,MG=32,JG=8,NWJ2=30'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	if L eq 'L50' then begin
printf,filnum,'C T10 L50  full sphere'
printf,filnum,'     PARAMETER(NN=10,MM=10,NHEM=2 ,NL=50,MOCT=1,MG=32,JG=8,NWJ2=30'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	if L eq 'L100' then begin
printf,filnum,'C T10 L100  full sphere'
printf,filnum,'     PARAMETER(NN=10,MM=10,NHEM=2 ,NL=100,MOCT=1,MG=32,JG=8,NWJ2=30'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit 
	endif
endif
	
;;;############ T21 #########
If Tnum eq 'T20' then begin
	if L eq 'L5' then begin
printf,filnum,'C T21 L5  full sphere'
printf,filnum,'     PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121'       
printf,filnum,'   &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif                       
	if L eq 'L10' then begin
printf,filnum,'C T21 L10  full sphere'
printf,filnum,'PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=10,MOCT=1,MG=64,JG=16,NWJ2=121'       
printf,filnum,'&         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit                       
		endif                       
	if L eq 'L20' then begin
printf,filnum,'C T21 L20  full sphere'
printf,filnum,'     PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=2 0,MOCT=1,MG=64,JG=16,NWJ2=121'       
printf,filnum,'   &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit                       
	endif                       
	if L eq 'L30' then begin
printf,filnum,'C T21 L30  full sphere'
printf,filnum,'     PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=30,MOCT=1,MG=64,JG=16,NWJ2=121'       
printf,filnum,'   &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit                        
	endif                       
	if L eq 'L40' then begin
printf,filnum,'C T21 L40  full sphere'
printf,filnum,'     PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=40,MOCT=1,MG=64,JG=16,NWJ2=121'       
printf,filnum,'   &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
    endif
endif

                       
;;;############ T31 #########
If Tnum eq 'T31' then begin
	if L eq 'L10' then begin
printf,filnum,'C T31 L10  full sphere'
printf,filnum,'     PARAMETER(NN=31,MM=31,NHEM=2 ,NL=10,MOCT=1,MG=96,JG=24,NWJ2=256'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	if L eq 'L24' then begin
printf,filnum,'C T31 L24  full sphere'
printf,filnum,'     PARAMETER(NN=31,MM=31,NHEM=2 ,NL=24,MOCT=1,MG=96,JG=24,NWJ2=256'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit 
	endif
	if L eq 'L30' then begin
printf,filnum,'C T31 L30  full sphere'
printf,filnum,'       PARAMETER(NN=31,MM=31,NHEM=2 ,NL=30,MOCT=1,MG=96,JG=24,NWJ2=256'        
printf,filnum,'     +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	if L eq 'L45' then begin
printf,filnum,'C T31 L45  full sphere'
printf,filnum,'     PARAMETER(NN=31,MM=31,NHEM=2 ,NL=45,MOCT=1,MG=96,JG=24,NWJ2=256'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
endif

;;;############ T42 #########
If Tnum eq 'T42'  then begin
	if L eq 'L15' then begin

printf,filnum,'C T42 L15  full sphere'
printf,filnum,'     PARAMETER(NN=42,MM=42,NHEM=2 ,NL=15,MOCT=1,MG=128,JG=32,NWJ2=462'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	if L eq 'L24' then begin
printf,filnum,'C T42 L24  full sphere'
printf,filnum,'     PARAMETER(NN=42,MM=42,NHEM=2 ,NL=24,MOCT=1,MG=128,JG=32,NWJ2=462'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	if L eq 'L30' then begin
printf,filnum,'C T42 L30  full sphere'
printf,filnum,'     PARAMETER(NN=42,MM=42,NHEM=2 ,NL=30,MOCT=1,MG=128,JG=32,NWJ2=462'        
printf,filnum,'   +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)' 
goto,quit
	endif
endif
 
;;;############ T63 #########
If Tnum eq 'T63' then begin
	if L eq 'L27' then begin

printf,filnum,'C T63 L27  full sphere'
printf,filnum,'     PARAMETER(NN=63,MM=63,NHEM=2 ,NL=2 7,MOCT=1,MG=192,JG=48'        
printf,filnum,'   +,NWJ2=1024,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
endif

;;;############ T85 #########
If Tnum eq 'T85' then begin
	if L eq 'L27' then begin

printf,filnum,'C T85 L27  full sphere'
printf,filnum,'     PARAMETER(NN=85,MM=85,NHEM=2 ,NL=2 7,MOCT=1,MG=256,JG=64,NWJ2=1849'      
printf,filnum,'   +,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif

	if L eq 'L20' then begin
printf,filnum,'C T85 L20  full sphere'
printf,filnum,'     PARAMETER(NN=85,MM=85,NHEM=2 ,NL=2 0,MOCT=1,MG=256,JG=64,NWJ2=1849 '     
printf,filnum,'   +,NCRAY=1,JGL=JG)'
goto,quit
	endif
endif

;;;############ T170 #########
If Tnum eq 'T170' then begin
	if L eq 'L27' then begin

printf,filnum,'C T170 L27  full sphere'
printf,filnum,'     PARAMETER(NN=170,MM=170,NHEM=2 ,NL=2 7,MOCT=1,MG=512,JG=128'
printf,filnum,'   +,NWJ2=7310,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'  
goto,quit
	endif
	if L eq 'L20' then begin
printf,filnum,'C T170 L20  full sphere'
printf,filnum,'     PARAMETER(NN=170,MM=170,NHEM=2 ,NL=2 0,MOCT=1,MG=512,JG=128,'
printf,filnum,'   & NWJ2=7310,NCRAY=1,JGL=JG)'
goto,quit
	endif
endif

;;;############ T213 #########

If Tnum eq 'T213' then begin
	if L eq 'L30' then begin

printf,filnum,'C T213 L30  full sphere'
printf,filnum,'     PARAMETER(NN=213,MM=213,NHEM=2 ,NL=30,MOCT=1,MG=1024,JG=160'
printf,filnum,'   +,NWJ2=11449,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'
goto,quit
	endif
	if L eq 'L60' then begin
printf,filnum,'C T213 L60  full sphere'
printf,filnum,'     PARAMETER(NN=213,MM=213,NHEM=2 ,NL=60,MOCT=1,MG=1024,JG=160,NWJ2=11449'  
printf,filnum,'   +,NCRAY=1,JGL=JG)'  
goto,quit
	endif
endif

if Tnum ne 'T2'  and Tnum ne 'T5' and Tnum ne 'T10' and Tnum ne 'T20' and Tnum ne 'T31' and Tnum ne 'T42'  $
and Tnum ne 'T63' and Tnum ne 'T85' and Tnum ne 'T170' and Tnum ne 'T213' then goto, YOUPICK

YOUPICK:

print,"Hit Y to choose, any other key for next choice:"
line1='T2 L20  full sphere'
line2='   PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=2 0,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T2 L30  full sphere'
line2= '   PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=30,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T2 L40  full sphere'
line2='   PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=40,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T2 L50  full sphere'
line2='   PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=50,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T2 L100  full sphere'
line2='   PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=100,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T2 L200  full sphere'
line2='   PARAMETER(NN=2 ,MM=2 ,NHEM=2 ,NL=2 00,MOCT=1,MG=8,JG=2 ,NWJ2=2'         
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T5 L20  full sphere'
line2='   PARAMETER(NN=5,MM=5,NHEM=2 ,NL=2 0,MOCT=1,MG=16,JG=4,NWJ2=9'       
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L5  full sphere'
line2='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=5,MOCT=1,MG=32,JG=8,NWJ2=30'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L10  full sphere'
line2='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=10,MOCT=1,MG=32,JG=8,NWJ2=30'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L20  full sphere'
line2='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=2 0,MOCT=1,MG=32,JG=8,NWJ2=30'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L30  full sphere'
line2='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=30,MOCT=1,MG=32,JG=8,NWJ2=30'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L36  full sphere'
line2='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=36,MOCT=1,MG=32,JG=8,NWJ2=30'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L50  full sphere'
line2='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=50,MOCT=1,MG=32,JG=8,NWJ2=30'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L100  full sphere'
line2='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=100,MOCT=1,MG=32,JG=8,NWJ2=30'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L5  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121'       
line3='  &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'                         
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L10  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=10,MOCT=1,MG=64,JG=16,NWJ2=121'       
line3='  &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'                         
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L20  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=2 0,MOCT=1,MG=64,JG=16,NWJ2=121'       
line3='  &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'                         
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L30  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=30,MOCT=1,MG=64,JG=16,NWJ2=121'       
line3='  &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'                         
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L40  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=40,MOCT=1,MG=64,JG=16,NWJ2=121'       
line3='  &         ,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'                         
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T31 L10  full sphere'
line2='   PARAMETER(NN=31,MM=31,NHEM=2 ,NL=10,MOCT=1,MG=96,JG=24,NWJ2=256'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T31 L24  full sphere'
line2='   PARAMETER(NN=31,MM=31,NHEM=2 ,NL=24,MOCT=1,MG=96,JG=24,NWJ2=256'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T31 L30  full sphere'
line2='  PARAMETER(NN=31,MM=31,NHEM=2 ,NL=30,MOCT=1,MG=96,JG=24,NWJ2=256'        
line3='    +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T31 L45  full sphere'
line2='   PARAMETER(NN=31,MM=31,NHEM=2 ,NL=45,MOCT=1,MG=96,JG=24,NWJ2=256'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T42 L15  full sphere'
line2='   PARAMETER(NN=42,MM=42,NHEM=2 ,NL=15,MOCT=1,MG=128,JG=32,NWJ2=462'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T42 L24  full sphere'
line2='   PARAMETER(NN=42,MM=42,NHEM=2 ,NL=24,MOCT=1,MG=128,JG=32,NWJ2=462'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T42 L30  full sphere'
line2='   PARAMETER(NN=42,MM=42,NHEM=2 ,NL=30,MOCT=1,MG=128,JG=32,NWJ2=462'        
line3='  +,NCRAY=8,JGL=JG,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T63 L27  full sphere'
line2='   PARAMETER(NN=63,MM=63,NHEM=2 ,NL=2 7,MOCT=1,MG=192,JG=48'        
line3='  +,NWJ2=1024,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T85 L27  full sphere'
line2='   PARAMETER(NN=85,MM=85,NHEM=2 ,NL=2 7,MOCT=1,MG=256,JG=64,NWJ2=1849'      
line3='  +,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'  

print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T170 L27  full sphere'
line2='   PARAMETER(NN=170,MM=170,NHEM=2 ,NL=2 7,MOCT=1,MG=512,JG=128'
line3='  +,NWJ2=7310,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T213 L30  full sphere'
line2='   PARAMETER(NN=213,MM=213,NHEM=2 ,NL=30,MOCT=1,MG=1024,JG=160'
line3='  +,NWJ2=11449,NCRAY=8,JGL=1,NTRAC=1,NLEVRF=1)'

print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T10 L5  full sphere'
line='   PARAMETER(NN=10,MM=10,NHEM=2 ,NL=5,MOCT=1,MG=32,JG=8,NWJ2=30'        
line='  +,NCRAY=1,JGL=JG)'  

print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

;;;############ T21 #########
line1='T21 L5  1/2 hemisphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=1,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121'        
line3='  +,NCRAY=1,JGL=JG)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L5  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121'        
line3='  +,NCRAY=1,JGL=JG)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L10  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=10,MOCT=1,MG=64,JG=16,NWJ2=121'        
line3='  +,NCRAY=1,JGL=JG)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T21 L1  full sphere'
line2='   PARAMETER(NN=2 1,MM=2 1,NHEM=2 ,NL=1,MOCT=1,MG=64,JG=16,NWJ2=121'        
line3='  +,NCRAY=1,JGL=JG)'  

print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T31 L10  full sphere'
line2='   PARAMETER(NN=31,MM=31,NHEM=2 ,NL=10,MOCT=1,MG=96,JG=24,NWJ2=256'        
line3='  +,NCRAY=1,JGL=JG)'  


;;;############ T42 #########

print,line1
print,line2
print,line3

kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T42 L15  full sphere'
line2='   PARAMETER(NN=42,MM=42,NHEM=2 ,NL=15,MOCT=1,MG=128,JG=32,NWJ2=462'        
line3='  +,NCRAY=1,JGL=JG)'  

print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T42 L40  full sphere'
line2='   PARAMETER(NN=42,MM=42,NHEM=2 ,NL=40,MOCT=1,MG=128,JG=32,NWJ2=462'        
line3='  +,NCRAY=1,JGL=JG)'  


;;;############ T63 #########
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T63 L20  full sphere'
line2='   PARAMETER(NN=63,MM=63,NHEM=2 ,NL=2 0,MOCT=1,MG=192,JG=48'        
line3='  +,NWJ2=1024,NCRAY=1,JGL=JG)'  

print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit


line1='T63 L40  full sphere' ; (seg faults on beehive with OMP, even with very large stacks)
line2='   PARAMETER(NN=63,MM=63,NHEM=2 ,NL=40,MOCT=1,MG=192,JG=48        
line3='  +,NWJ2=1024,NCRAY=1,JGL=JG)'  


print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T85 L20  full sphere'
line2='   PARAMETER(NN=85,MM=85,NHEM=2 ,NL=2 0,MOCT=1,MG=256,JG=64,NWJ2=1849'      
line3='  +,NCRAY=1,JGL=JG)'



print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

line1='T170 L20  full sphere'
line2='   PARAMETER(NN=170,MM=170,NHEM=2 ,NL=2 0,MOCT=1,MG=512,JG=128,'
line3='  & NWJ2=7310,NCRAY=1,JGL=JG)'


print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit


line1='T106 L30  full sphere'
line2='   PARAMETER(NN=106,MM=106,NHEM=2 ,NL=30,MOCT=1,MG=512,JG=80,NWJ2=2862'      
line3='  +,NCRAY=1,JGL=JG)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

;;;############ T213 #########

line1='T213 L60  full sphere'
line2='   PARAMETER(NN=213,MM=213,NHEM=2 ,NL=60,MOCT=1,MG=1024,JG=160,NWJ2=11449'   
line3='  +,NCRAY=1,JGL=JG)'  
print,line1
print,line2
print,line3
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then goto, chosen
if kybrd eq 'q' or kybrd eq 'Q' then goto, quit

CHOSEN:
print,''
print,''
print,'You chose: '
print,line1
print,line2
print,line3
print,''
print,''
print,'write to params.i file?'
kybrd=get_kbrd(/key_name)
if kybrd eq 'y' or kybrd eq 'Y' then begin
printf,filnum,'C     '+strtrim(string(line1),2)
printf,filnum,'      '+strtrim(string(line2),2)
printf,filnum,'     '+strtrim(string(line3),2)
print,'wrote file: params.i'
endif

QUIT:
close,filnum
print,'(closed file:',strtrim(filnum,2),')'

end
