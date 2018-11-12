
PRO ebv_schlegel,ra,dec;,filter,mag

if ra lt 0 or ra gt 360 or dec lt -90 or dec gt 90 then begin
   print,'Wrong Coordinates - R.A. = [0,360]; Dec. = [-90,90]'
   stop
endif

; Schlegel maps are in galactic coordinates - 
; First transformation from R.A./Dec. (J2000) to galactic - 
radectolb,ra,dec,l,b

; Convert to the pixel values in the dust maps - 
if b ge 0 then begin
dust_wcs2xy,'n',l,b,x,y ;North
;x2 = 2048 * sqrt(1 - sin(b*180./!pi)) * cos(l*180./!pi) + 2047.5
;y2 = -2048 * sqrt(1 - sin(b*180./!pi)) * sin(l*180./!pi) + 2047.5
;print,'Formula XY: ',x2,y2
map = 'n'
endif else begin
dust_wcs2xy,'s',l,b,x,y ;South
;x2 = 2048 * sqrt(1 + sin(b*180./!pi)) * cos(l*180./!pi) + 2047.5
;y2 = 2048 * sqrt(1 + sin(b*180./!pi)) * sin(l*180./!pi) + 2047.5
;print,'Formula XY: ',x2,y2
map = 's'
endelse
x = (x > 0) < 4095 ; Force pixel locations to fall within the image bounds [0,4095]
y = (y > 0) < 4095 ; Force pixel locations to fall within the image bounds [0,4095]

; Extract one value in the dust maps
path = '/Users/audreygalametz/SOFT/IDL/SCHLEGEL/SFD_4096/'
if map eq 'n' then begin
   map = 'SFD_dust_4096_ngp.fits'
endif else begin
   map = 'SFD_dust_4096_sgp.fits'
endelse
hdr=headfits(path+map)
x1 = (x-1) > 0
x2 = (x+1) < 4095
y1 = (y-1) > 0
y2 = (y+1) < 4095
fxread, path+map, value, hdr, x1, x2, y1, y2
value = median(value)
print,'Schlegel E(B-V) from fits: ',value

map = mrdfits('SCHLEGEL/lambda_sfd_ebv.fits',1,/silent)
radectolb,ra,dec,l,b
phi = l / (180./!pi)
theta = (90 - b) / (180./!pi)
ang2pix_nest,512,theta,phi,id
ebv = map[id].TEMPERATURE
print,'Schlegel E(B-V) from Healpix: ',ebv

;!!!!

; Multiplying by the A/E(B-V) coefficients
;readcol,'dust_schlegel_coefipac.txt',band,wave,aebvSF,aebvSc,format='A,F,F,F'
;wave = wave *1e4
;;SandF: Schlafly and Finkbeiner 2011 (ApJ 737, 103)
;;SFD: Schlegel et al. 1998 (ApJ 500, 525)
;mag = value * interpol(aebvSc,wave,filter)
;print,'Mag. Corr. at '+strtrim(string(filter),1)+'A: ',mag

END

PRO radectolb,ra,dec,l,b

; Definition from http://terraformers.info/files/galactic.pdf
  radeg = 180.0d/!DPI           ; Conversion to radians
  ras = ra / radeg
  decs = dec /radeg
  rapol = (12d0 + 51.4/60.0d0) * 15.0d0 ; RA of North Pole
  rapol = rapol / radeg
  decpol = 27.0d0 + 77.0d0/600.0d0 ; Dec of North Pole
  decpol = decpol / radeg
  sdecs = sin(decs)
  cdecs = cos(decs)
  sdecpol = sin(decpol)
  cdecpol = cos(decpol)
  b1 = sdecs * sdecpol + cdecs*cdecpol*cos(ras-rapol)
  b = asin(b1) * radeg
  cb = cos(b/radeg)
  racen = (17d0 + 45.6/60.0d0) * 15.0d0 ; RA of Galactic Center
  racen = racen / radeg
  deccen = atan(-cos(racen-rapol)/tan(decpol)) ;* radeg ; Dec Galactic Center = -28.929 656 ...
;  deccen = deccen / radeg
  sdeccen = sin(deccen)
  j = ((sdecs*cdecpol)-(cdecs*sdecpol*cos(ras-rapol))) / cb * radeg
  k = asin((cdecs*sin(ras-rapol)) / cb) * radeg
  q = acos(sdeccen/cdecpol) * radeg
  if j lt 0 then begin
     l = q + k - 180.0d0
  endif else begin
     l = q - k
  endelse
  if l lt 0 then l = l + 360

print,'Galactic coordinates: ',l,b

END

PRO dust_wcs2xy,flag,l,b,xpix,ypix

;path = '/Users/audreygalametz/SOFT/IDL/SCHLEGEL/SFD_4096/'
  crpix = 2048.5
  ;lonpole = 180.
  dradeg = 180.d0 / !dpi
  
  if flag eq 'n' then begin
;mapNorth = 'SFD_dust_4096_ngp.fits' - crval = 90.
     theta = b
     phi = 270. + l             ; l + 180 + lonpole - crval
     cd1_1 = -0.039564682
     cd2_2 = 0.039564682
  endif
  if flag eq 's' then begin
;mapSouth = 'SFD_dust_4096_sgp.fits' - crval = -90
     theta = -b
     phi = 90. - l             ; l + 180 + lonpole - crval
     cd1_1 = 0.039564682
     cd2_2 = -0.039564682
  endif
  
  phi = (phi MOD 360) +360 * (phi LT 0) ; phi between [0,360]
  rtheta = 2 * dradeg * sin((90 - theta)/2. / dradeg)
  xtemp = rtheta * sin(phi / dradeg)
  ytemp = -rtheta * cos(phi / dradeg)
;  denom = cd1_1 * cd2_2                        ;- cd1_2 * cd2_1
  xpix = (xtemp/cd1_1) + (crpix - 1.0)  ;(cd2_2*xtemp - cd1_2*ytemp) / denom + (crpix - 1.0)
  ypix = (ytemp/cd2_2) + (crpix - 1.0)  ;(cd1_1*ytemp - cd2_1*xtemp) / denom + (crpix - 1.0)
  xpix = fix(xpix + 0.5)                       ;Round to nearest pixel
  ypix = fix(ypix + 0.5)                       ;Round to nearest pixel
  
;  if flag eq 'n' then print,'X,Y (North Map): ', xpix,ypix
;  if flag eq 's' then print,'X,Y (South Map): ', xpix,ypix
  
END
