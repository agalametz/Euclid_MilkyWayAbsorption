
PRO ebv_schlafly,ra,dec

map = mrdfits('SCHLAFLY/ps1-ebv-4.5kpc.fits',1,/silent)

radectolb,ra,dec,l,b
phi = l / (180./!pi)
theta = (90 - b) / (180./!pi)
ang2pix_ring,512,theta,phi,id
ebv = map[id].ebv
print,'Schlafly E(B-V): ',ebv

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
