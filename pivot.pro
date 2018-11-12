
PRO pivot

reddeninglaw,Lambda,kL

;readcol,'/Users/audreygalametz/WORK/EUCLID/OU-PHZ-3-205/3-TFA_TESTS/FILTERS/'+filter+'.txt',wf,ff
;if filter eq 'VIS' then wf = wf *10.
;if filter eq 'LSSTu' then wf = wf*10.
;if filter eq 'LSSTz' then wf = wf*10.
;readcol,'/Users/audreygalametz/Phosphoros/AuxiliaryData/Filters/PANSTARRS/filter.ps1_gband_2.724',wf,ff
;readcol,'/Users/audreygalametz/Phosphoros/AuxiliaryData/Filters/PANSTARRS/filter.ps1_rband_2.724',wf,ff
;readcol,'/Users/audreygalametz/Phosphoros/AuxiliaryData/Filters/PANSTARRS/filter.ps1_iband_2.724',wf,ff
;readcol,'/Users/audreygalametz/Phosphoros/AuxiliaryData/Filters/PANSTARRS/filter.ps1_zband_2.724',wf,ff
;readcol,'/Users/audreygalametz/Phosphoros/AuxiliaryData/Filters/PANSTARRS/filter.ps1_yband_2.725',wf,ff
;readcol,'/Users/audreygalametz/WORK/EUCLID/OU-PHZ-3-205/3-TFA_TESTS/FILTERS/ThroughputVIS.txt',wf,ff - wf=wf*10.

readcol,'/Users/audreygalametz/Phosphoros/AuxiliaryData/Filters/PANSTARRS/filter.ps1_gband_2.724',wf,ff

int = trap_int(wf,wf*ff)
int2 = trap_int(wf,ff/wf)
pivot = sqrt(int/int2)
kpivot = interpol(kL,Lambda,pivot,SPLINE=spline)
print,'pivot: ',pivot
print,'kpivot: ',kpivot
print,'Dm at Ebv=0.1: ',kpivot*0.1   
print,'Dm at Ebv=0.3: ',kpivot*0.3
print,'Dm at Ebv=0.02: ',kpivot*0.02

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Definition of Fitzpatrick 99
; Reddening curve k(Lambda) = A(Lambda)/E(B-V)
; Rv-dependent 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO reddeninglaw,Lambda,kL                                                                                                                                                                   

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

END
