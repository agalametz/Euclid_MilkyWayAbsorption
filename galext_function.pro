

PRO galext_function;,filter

reddeninglaw,Lambda,kL
readcol,'/Users/audreygalametz/WORK/EUCLID/OU-PHZ/4-Phosphoros/GUI_Phosphoros/RedCurve/calzetti.dat',Lambdacalz,calzinit

;sed = ['/Users/audreygalametz/WORK/EUCLID/OU-PHZ/1-PAST_ALGORITHMS/LePhare/lephare_dev/sed/GAL/COSMOS_SED/'+['Ell1_A_0.sed','Ell2_A_0.sed','Ell3_A_0.sed','Ell4_A_0.sed','Ell5_A_0.sed','Ell6_A_0.sed','Ell7_A_0.sed','S0_A_0.sed','Sa_A_0.sed','Sa_A_1.sed','Sb_A_0.sed','Sb_A_1.sed','Sc_A_0.sed','Sc_A_1.sed','Sc_A_2.sed','Sd_A_0.sed','Sd_A_1.sed','Sd_A_2.sed','Sdm_A_0.sed','SB0_A_0.sed','SB1_A_0.sed','SB2_A_0.sed','SB3_A_0.sed','SB4_A_0.sed','SB5_A_0.sed','SB6_A_0.sed','SB7_A_0.sed','SB8_A_0.sed','SB9_A_0.sed','SB10_A_0.sed','SB11_A_0.sed']]

sed = '/Users/audreygalametz/WORK/EUCLID/OU-PHZ/1-PAST_ALGORITHMS/LePhare/lephare_dev/sed/GAL/COSMOS_SED/Ell1_A_0.sed'
;sed = '/Users/audreygalametz/WORK/EUCLID/OU-PHZ/1-PAST_ALGORITHMS/LePhare/lephare_dev/sed/GAL/COSMOS_SED/SB11_A_0.sed'
readcol,sed,w,f
;Intebv = 0
;readcol,'/Users/audreygalametz/WORK/EUCLID/OU-PHZ/4-Phosphoros/GUI_Phosphoros/RedCurve/calzetti.dat',Lambdacalz,calzinit
;fc = interpol(Lambdacalz,calzinit,w,SPLINE=spline)
;fr = f * 10^(-0.4*fc*Intebv)
fr = f

z = 0.1
wr = w * (1+z)

fd = interpol(kL,Lambda,wr,SPLINE=spline)

;;;;;;;;;;;;;;;;;;;;;;;;

;e.g.,filter = ['VIS','DESu','Hnir_WFC3f160w']
path = '/Users/audreygalametz/WORK/EUCLID/OU-PHZ/3-TFA_TESTS/FILTERS/'
filter = 'LSSTg'
filtername = 'LSST g-band'
;if filter eq 'Y' then filter = 'Ynir_WFC3f105w'
;if filter eq 'J' then filter = 'Jnir_WFC3f125w'
;if filter eq 'H' then filter = 'Hnir_WFC3f160w'
readcol,path+filter+'.txt',wf,ff
wf = wf * 10. ;(if filter eq 'LSSTg')

col = mbplotpars('galext_function.ps',/ps,/big)

multiplot,[1,1]

size = 31
Ebv = findgen(size)/100.
magcorr = fltarr((size(Ebv))[1])

integfilter,wr,fr,wf,ff,mag

for ii = 0,(size(Ebv))[1]-1 do begin
   fext = fr * 10^(-0.4*fd*Ebv[ii])
   integfilter,wr,fext,wf,ff,mage
   magcorr[ii] = mage - mag
endfor

;plot,Ebv,magcorr
a=(magcorr[size-1]-magcorr[0])/(Ebv[size-1]-Ebv[0])
;a=(magcorr[10]-magcorr[0])/(Ebv[10]-Ebv[0])
plot,Ebv,a*Ebv,ytitle='A!ig!n',charsize=1.2,position=[0.15,0.35,0.95,0.8];xtitle='pD'
xyouts,0.02,magcorr[size-1]*7./8.,filtername,charsize=1.2

;magcorr2 = fltarr((size(Ebv))[1],10)
;for ii = 0,(size(Ebv))[1]-1 do begin
;   x = -0.4 * alog(10) * Ebv[ii] * fd
;   fext = fr * (1 + x)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,0] = mage - mag
;   fext = fr * (1 + x + x^2/2.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,1] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,2] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6. + x^4/24.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,3] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6. + x^4/24. + x^5/120.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,4] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6. + x^4/24. + x^5/120. + x^6/720.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,5] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6. + x^4/24. + x^5/120.+ x^6/720. + x^7/5040.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,6] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6. + x^4/24. + x^5/120.+ x^6/720. + x^7/5040. + x^8/40320.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,7] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6. + x^4/24. + x^5/120.+ x^6/720. + x^7/5040. + x^8/40320. + x^9/362880.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,8] = mage - mag
;   fext = fr * (1 + x + x^2/2. + x^3/6. + x^4/24. + x^5/120.+ x^6/720. + x^7/5040. + x^8/40320. + x^9/362880. + x^10/3628800.)
;   integfilter,wr,fext,wf,ff,mage
;   magcorr2[ii,9] = mage - mag
;endfor
;for ii=0,9 do oplot,Ebv,magcorr2[*,ii],linestyle=1

;oplot,Ebv,magcorr,color=fsc_color('red')

;if fil eq (size(filter))[1]-1 then begin
;plot,Ebv,magcorr[*,0],xtitle='E(B-V)',ytitle='!9D!3 mag',xr=[0,5],xs=1,yr=[0,max(magcorr[50,*])+0.1],ys=1
;endif else begin
;plot,Ebv,magcorr[*,0],ytitle='!9D!3 mag',xr=[0,5],xs=1,yr=[0,max(magcorr[50,*])+0.1],ys=1
;endelse
;xyouts,0.05,(max(!Y.CRANGE)*3.3/4.),filter[fil]
;xyouts,0.05,(max(!Y.CRANGE)*2.3/4.),'z = [0,6]'

x = 0.4 * alog(10) * fd
integfilter,wr,fr,wf,ff,m0,ord0
integfilter,wr,fr*x,wf,ff,m,ord1
integfilter,wr,fr*(x^2)/2.,wf,ff,m,ord2
integfilter,wr,fr*(x^3)/factorial(3),wf,ff,m,ord3
integfilter,wr,fr*(x^4)/factorial(4),wf,ff,m,ord4
integfilter,wr,fr*(x^5)/factorial(5),wf,ff,m,ord5
integfilter,wr,fr*(x^6)/factorial(6),wf,ff,m,ord6
integfilter,wr,fr*(x^7)/factorial(7),wf,ff,m,ord7
integfilter,wr,fr*(x^8)/factorial(8),wf,ff,m,ord8
integfilter,wr,fr*(x^9)/factorial(9),wf,ff,m,ord9
integfilter,wr,fr*(x^10)/factorial(10),wf,ff,m,ord10

magcorr3 = fltarr((size(Ebv))[1],10)
for ii = 0,(size(Ebv))[1]-1 do begin
;   magcorr3[ii,0] = -2.5*alog10(ord0) - m0
   magcorr3[ii,0] = -2.5*alog10(ord0 - ord1*Ebv[ii]) - m0
   magcorr3[ii,1] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2) - m0
   magcorr3[ii,2] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3) - m0
   magcorr3[ii,3] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3 + ord4*Ebv[ii]^4) - m0
   magcorr3[ii,4] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3 + ord4*Ebv[ii]^4 - ord5*Ebv[ii]^5) - m0
   magcorr3[ii,5] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3 + ord4*Ebv[ii]^4 - ord5*Ebv[ii]^5 + ord6*Ebv[ii]^6) - m0
   magcorr3[ii,6] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3 + ord4*Ebv[ii]^4 - ord5*Ebv[ii]^5 + ord6*Ebv[ii]^6 - ord7*Ebv[ii]^7) - m0
   magcorr3[ii,7] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3 + ord4*Ebv[ii]^4 - ord5*Ebv[ii]^5 + ord6*Ebv[ii]^6 - ord7*Ebv[ii]^7 + ord8*Ebv[ii]^8) - m0
   magcorr3[ii,8] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3 + ord4*Ebv[ii]^4 - ord5*Ebv[ii]^5 + ord6*Ebv[ii]^6 - ord7*Ebv[ii]^7 + ord8*Ebv[ii]^8 - ord9*Ebv[ii]^9) - m0
   magcorr3[ii,9] = -2.5*alog10(ord0 - ord1*Ebv[ii] + ord2*Ebv[ii]^2 - ord3*Ebv[ii]^3 + ord4*Ebv[ii]^4 - ord5*Ebv[ii]^5 + ord6*Ebv[ii]^6 - ord7*Ebv[ii]^7 + ord8*Ebv[ii]^8 - ord9*Ebv[ii]^9 + ord10*Ebv[ii]^10) - m0
endfor
for ii = 0,4 do oplot,Ebv,magcorr3[*,ii],linestyle=1

oplot,Ebv,magcorr,color=fsc_color('red')

oplot,0.05*[1,1],[0.165883,0.187566],color=fsc_color('dark gray')
oplot,0.1*[1,1],[0.331432,0.374752],color=fsc_color('dark gray')
oplot,0.15*[1,1],[0.496639,0.561560],color=fsc_color('dark gray')
oplot,0.2*[1,1],[0.661497,0.747988],color=fsc_color('dark gray')
oplot,0.25*[1,1],[0.825995,0.934035],color=fsc_color('dark gray')

xyouts,0.18,1.1,'(1)',charsize=1
xyouts,0.28,0.675,'(2)',charsize=1
xyouts,0.27,1.1,'(3)',charsize=1
xyouts,0.28,0.92,'(4)',charsize=1
xyouts,0.285,1.052,'(5)',charsize=1

;p=[0.1,0.1]
;std=[0.1,0.1]
;yfit = mpcurvefit(Ebv,magcorr,std,p,sig,FUNCTION_NAME='YeqAX',/autoderivative,/quiet,CHISQ=chisq)
;print,p
;oplot,Ebv,p[0]*Ebv+p[1]*Ebv*Ebv,color=fsc_color('green')

plot,Ebv,magcorr - (a*Ebv),position=[0.15,0.1,0.95,0.3],charsize=1.2,yr=[-0.02,0.02],xtitle='pD',ytitle='A!ig!n real - model';,ytickname=['-1e!e-2!n','-5e!e-3!n','0','5e!e-3!n','1e!e-2!n']
oplot,[0,5],0*[1,1],color=fsc_color('red')
for ii = 0,4 do oplot,Ebv,magcorr - magcorr3[*,ii],linestyle=1
oplot,0.05*[1,1],magcorr[5]-[0.165883,0.187566],color=fsc_color('dark gray')
oplot,0.1*[1,1],magcorr[10]-[0.331432,0.374752],color=fsc_color('dark gray')
oplot,0.15*[1,1],magcorr[15]-[0.496639,0.561560],color=fsc_color('dark gray')
oplot,0.2*[1,1],magcorr[20]-[0.661497,0.747988],color=fsc_color('dark gray')
oplot,0.25*[1,1],magcorr[25]-[0.825995,0.934035],color=fsc_color('dark gray')

xyouts,0.055,-0.018,'(1)',charsize=1
xyouts,0.102,0.014,'(2)',charsize=1
xyouts,0.212,-0.018,'(3)',charsize=1
xyouts,0.263,0.014,'(4)',charsize=1
xyouts,0.285,-0.008,'(5)',charsize=1

device,/close

END

PRO integfilter,wi,fi,wf,ff,mout,fout

linterp, wf,ff,wi,fff,missing=0.
i1=min(where(fff gt 0.05))
i2=max(where(fff gt 0.05))
;ic=floor((i1+i2)/2.)
scale = 1.;abs(fi[ic])
nu = 2.99792456211d18/wi
fnu = wi*0.+3631*1.d-23
fnu2fl,nu,fnu,wz,flam
k1 = scale*trap_int(wi,(wi/1.d3)*(fi/scale)*fff)
k2 = 1.d-6*trap_int(wi,(wi/1.d3)*(flam*1.d9)*fff)
;fout = -2.5*alog10(K1/K2)-7.5
mout = -2.5*alog10(K1/K2)-7.5
fout = K1/K2*1e3

;ffn = interpol(ff,wf,wi,SPLINE=spline)
;fout = fi * ffn
;plot,wi,fi,xr=[1e3,1e5],xs=1,/xlog
;oplot,wf,ff
;oplot,wi,fout
;fout = int_tabulated(wi,fout,/sort)

END

PRO reddeninglaw,Lambda,kL
  
; Reddening curve k(Lambda) = A(Lambda)/E(B-V)
; Fitzpatrick 1999 (MW) - Rv-dependent                                                                                                                                           
x = findgen(1000)/100. + 0.01
kL = fltarr(1000)
Lambda = 1./(x/1e4)

; For Lambda < 2700A
a=where(1./(x/1e4) le 2700)
Rv = 3.1
C2 = -0.824 + 4.717/Rv
C1 = 2.030 - 3.007*C2
xo = 4.596 ;(um^-1)
gamma = 0.99 ;(um^-1)
C3 = 3.23
C4 = 0.41
D = x[a]^2 / ( (x[a]^2 - xo^2)^2 + (x[a]*gamma)^2)
F = fltarr((size(a))[1])
b=where(x[a] ge 5.9)
F[b] = (0.5392 * (x[a[b]]-5.9)^2) + (0.05644 * (x[a[b]]-5.9)^3)
kL[a] = C1 + C2*x[a] + C3*D + C4*F + Rv

; Interpolation at Lambda > 2700A
a = where(1./(x/1e4) gt 2700.)
wi = [0,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846]
kLi = [0,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591]
kL[a] = interpol(kLi,wi,x[a],/spline)

;plot,1./(x/1e4),kL,xtitle = '1/Lambda (um^-1)',ytitle='k(Lambda)',xr=[0,1e4],xs=1
;oplot,1./[0,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846]*1e4,$
;[0,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591],psym=1 
;col = mbplotpars('Fitzpatrick1999.ps',/ps,/big)
;plot,Lambda,kL,xr=[1e3,1e4],xs=1,/xlog,xtitle = 'Wavelenght (A)',ytitle = 'A / E(B-V)',charsize=1,position=[0.1,0.1,0.95,0.7]
;readcol,'/Users/audreygalametz/WORK/EUCLID/OU-PHZ/4-Phosphoros/1-DC1/FILTERS/MER/Rext_ACSf606w.txt',wf,ff
;oplot,wf,ff*20.,color=fsc_color('gray')
;int = trap_int(wf,wf*ff)
;int2 = trap_int(wf,ff/wf)
;pivot = sqrt(int/int2)
;kpivot = interpol(kL,Lambda,pivot,SPLINE=spline)
;plotsym,0,1,/fill
;oplot,pivot*[1,1],kpivot*[1,1],psym=8
;print,kpivot
;device,/close

END

;PRO YeqAX, X , P , YMOD, dP
;     YMOD = P[0]*X + P[1]*X*X
;END
