



PRO ebv_planck;,ra,dec

ra = [0, 90,180, 30,60,90,120,150,60, 192.84,  266.40]
dec= [-90,-60,-30,  0,10,20, 30, 60,90,27.1283, -28.929]
;l=0
;b=0
for ii=0,(size(ra))[1]-1 do begin
;   radectolb,ra[ii],dec[ii],l,b
;   print,l,b
;endfor
;stop

;read_fits_map,'HFI_CompMap_ThermalDustModel_2048_R1.20.fits',ebv
;ebv = reform(ebv[*,2])
;writefits,'PlanckEbv.fits',ebv
;read_fits_map,'HFI_CompMap_ThermalDustModel_2048_R1.20.fits',tau353
;tau353 = reform(tau353[*,2])
;writefits,'PlanckTau353.fits',tau353
;read_fits_map,'HFI_CompMap_ThermalDustModel_2048_R1.20.fits',radiance
;radiance = reform(radiance[*,2])
;writefits,'PlanckRadiance.fits',radiance

l=0
b=0
ebv = readfits('PLANCK/PlanckEbv.fits',/silent)
radectolb,ra[ii],dec[ii],l,b
phi = l / (180./!pi)
theta = (90 - b) / (180./!pi)
ang2pix_nest,2048,theta,phi,id
;print,'Planck E(B-V): ',ebv[id]
print,ebv[id]

endfor

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

;  print,'Galactic coordinates: ',l,b

END
