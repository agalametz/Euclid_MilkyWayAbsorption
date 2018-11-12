
PRO euclidFoV_ebv

flag = 'plot'

if flag eq 'disc' then begin

col = mbplotpars('euclidFoV_ebv_disc.ps',/ps,/big)
loadct,0
multiplot,[1,1]

; North
schlegel = readfits('DUSTMAPS/SCHLEGEL/SFD_dust_4096_ngp_euclid.fits')
;schlegel = schlegel /0.92 - factor in Planck 2014
planck = readfits('DUSTMAPS/PLANCK/Planck_ngp_euclid.fits')
;plot,planck,schlegel

im = fltarr(120,120)

for ii = 0,4095 do begin
   for jj = 0,4095 do begin
      if planck[ii,jj] ne 0 and schlegel[ii,jj] ne 0 then begin
         x = floor(planck[ii,jj]*200.)
         y = floor(schlegel[ii,jj]*200.)
         if x lt 119 and y lt 119 then im[x,y] = im[x,y]+1
      endif
   endfor
endfor
plot,0*[1,1],0*[1,1],/isotropic,xr=[0,0.6],yr=[0,0.6],xs=1,ys=1,position=[0.1,0.2,0.5,0.6],charsize=1.2,xtitle='Planck E!iB-V!n',ytitle='SFD E!iB-V!n' 
for ii = 0,119 do begin
   for jj = 0,119 do begin
      if im[ii,jj] ne 0 then begin
         if im[ii,jj] lt 5000000 then begin
            polyfill,[ii,ii+1,ii+1,ii,ii]/200.,[jj,jj,jj+1,jj+1,jj]/200.,color=200 - (200./alog10(500000)*alog10(im[ii,jj]))
         endif else begin
            polyfill,[ii,ii+1,ii+1,ii,ii]/200.,[jj,jj,jj+1,jj+1,jj]/200.,color=10
         endelse
      endif
   endfor
endfor
oplot,[0,1],[0,1],linestyle=1
xyouts,0.03,0.54,'North',charsize=1.2
polyfill,[3,4,4,3,3]/200.,[2,2,3,3,2]/200.,color=10
polyfill,[4,5,5,4,4]/200.,[3,3,4,4,3]/200.,color=10
; legend >> plot,0*[1,1],0*[1,1],position=[0.65,0.1,0.68,0.6],xr=[0,1],yr=[0,200]

schlegel = readfits('DUSTMAPS/SCHLEGEL/SFD_dust_4096_sgp_euclid.fits')
planck = readfits('DUSTMAPS/PLANCK/Planck_sgp_euclid.fits')
im = fltarr(120,120)

for ii = 0,4095 do begin
   for jj = 0,4095 do begin
      if planck[ii,jj] ne 0 and schlegel[ii,jj] ne 0 then begin
         x = floor(planck[ii,jj]*200.)
         y = floor(schlegel[ii,jj]*200.)
         if x lt 119 and y lt 119 then im[x,y] = im[x,y]+1
      endif
   endfor
endfor
;print,max(im)

plot,0*[1,1],0*[1,1],/isotropic,xr=[0,0.6],yr=[0,0.6],xs=1,ys=1,position=[0.55,0.2,0.95,0.6],charsize=1.2,xtitle='Planck E!iB-V!n'
for ii = 0,119 do begin
   for jj = 0,119 do begin
      if im[ii,jj] ne 0 then begin
         if im[ii,jj] lt 500000 then begin
            polyfill,[ii,ii+1,ii+1,ii,ii]/200.,[jj,jj,jj+1,jj+1,jj]/200.,color=200 - ((200)/alog10(500000)*alog10(im[ii,jj]))
         endif else begin
            polyfill,[ii,ii+1,ii+1,ii,ii]/200.,[jj,jj,jj+1,jj+1,jj]/200.,color=10
         endelse
      endif
   endfor
endfor
oplot,[0,1],[0,1],linestyle=1
xyouts,0.03,0.54,'South',charsize=1.2

plot,0*[1,1],0*[1,1],position=[0.2,0.65,0.85,0.67],xr=[0,alog10(500000)],yr=[0,1],xticks=5,yticks=1,ytickname=[' ',' '],xtickname=['5e0-','5e1','5e2','5e3','5e4','5e5+'],xs=1,ys=1,charsize=0.7,xtitle='# pixels in SFD map'
for ii = 0,113 do begin
   polyfill,[ii,ii+1,ii+1,ii,ii]/20.,[0,0,1,1,0],color=200 - (200./alog10(500000)*ii/20.)
endfor

device,/close

endif

if flag eq 'plot' then begin

col = mbplotpars('euclidFoV_ebv.ps',/ps,/big)

multiplot,[1,1]

loadct,0
read_jpeg,'DUSTMAPS/PLANCK/Planck_ngp_euclid_line.jpg',map
plotimage,bytscli(map[0,*,*]),position=[0.01,0.41,0.35,0.75],/preserve,/noerase,xs=4,ys=4
read_jpeg,'DUSTMAPS/PLANCK/Planck_sgp_euclid_line.jpg',map
plotimage,bytscli(map[0,*,*]),position=[0.01,0.06,0.35,0.4],/preserve,/noerase,xs=4,ys=4

readcol,'DUSTMAPS/euclidFoV_ebv_schlegelhealpix_histo.cat',bin,hn,hs
hnsg = hn * 100.
hssg = hs * 100.
;readcol,'DUSTMAPS/euclidFoV_ebv_schlaflyhealpix_histo.cat',bin,hn,hs
;hnsf = hn * 100.
readcol,'DUSTMAPS/euclidFoV_ebv_planckhealpix_histo.cat',bin,hn,hs
hnp = hn * 100.
hsp = hs * 100.

; North

plot,bin,hnp,psym=3,xr=[0,2],xs=1,ys=1,ytitle='% of HEALpix pixels',position=[0.48,0.47,0.97,0.73],charsize=1.5,yr=[1e-4,1e2],/ylog,ytickname=['1e-4','1e-3','1e-2','1e-1','1e0','1e1','1e2'],xtitle='E!iB-V!n'
xyouts,0.8,1e1,'North (b>30; |Dec.|>15)',charsize=1.2
oplot,[1.4,1.5],1e0*[1,1]
xyouts,1.55,1e0,'Planck',charsize=1.2
oplot,[1.4,1.5],1.3e-1*[1,1],color=fsc_color('slate gray')
xyouts,1.55,1.3e-1,'SFD',charsize=1.2

for ii = 0 , n_elements(bin)-1 do oplot,bin[ii]+[-0.005,0.005],hnp[ii]*[1,1]
for ii = 1 , n_elements(bin)-1 do oplot,bin[ii]-0.005*[1,1],[hnp[ii-1],hnp[ii]]
for ii = 0 , n_elements(bin)-2 do oplot,bin[ii]+0.005*[1,1],[hnp[ii],hnp[ii+1]]
a2 = where(bin ne 0)
a1 = where(bin gt 0.1)
print,"Planck North Map - % above Ebv=0.1: ",total(hnp[a1]) * 100. / total(hnp[a2])
;n_elements(a1)*1./n_elements(a2)*100.                         
a1 = where(bin gt 0.3)
print,"Planck North Map - % above Ebv=0.3: ",total(hnp[a1]) * 100. / total(hnp[a2])
;n_elements(a1)*1./n_elements(a2)*100.                         

for ii = 0 , n_elements(bin)-1 do oplot,bin[ii]+[-0.005,0.005],hnsg[ii]*[1,1],color=fsc_color('slate gray')
for ii = 1 , n_elements(bin)-1 do oplot,bin[ii]-0.005*[1,1],[hnsg[ii-1],hnsg[ii]],color=fsc_color('slate gray')
for ii = 0 , n_elements(bin)-2 do oplot,bin[ii]+0.005*[1,1],[hnsg[ii],hnsg[ii+1]],color=fsc_color('slate gray')
a2 = where(bin ne 0)
a1 = where(bin gt 0.1)
print,"Schlegel North Map - % above Ebv=0.1: ",total(hnsg[a1]) * 100. / total(hnsg[a2])
a1 = where(bin gt 0.3)
print,"Schlegel North Map - % above Ebv=0.3: ",total(hnsg[a1]) * 100. / total(hnsg[a2])

;for ii = 0 , n_elements(bin)-1 do oplot,bin[ii]+[-0.005,0.005],hnsf[ii]*[1,1],color=fsc_color('gray')
;for ii = 1 , n_elements(bin)-1 do oplot,bin[ii]-0.005*[1,1],[hnsf[ii-1],hnsf[ii]],color=fsc_color('gray')
;for ii = 0 , n_elements(bin)-2 do oplot,bin[ii]+0.005*[1,1],[hnsf[ii],hnsf[ii+1]],color=fsc_color('gray')
;a2 = where(bin ne 0)
;a1 = where(bin gt 0.1)
;print,"Schlafly North Map - % above Ebv=0.1: ",total(hnsf[a1]) * 100. / total(hnsf[a2])
;a1 = where(bin gt 0.3)
;print,"Schlafly North Map - % above Ebv=0.3: ",total(hnsf[a1]) * 100. / total(hnsf[a2])

; South

plot,bin,hsp,psym=3,xr=[0,2],xs=1,ys=1,ytitle='% of HEALpix pixels',position=[0.48,0.12,0.97,0.38],charsize=1.5,yr=[1e-4,1e2],/ylog,ytickname=['1e-4','1e-3','1e-2','1e-1','1e0','1e1','1e2'],xtitle='E!iB-V!n'
xyouts,0.8,1e1,'South (b<-30; |Dec.|>15)',charsize=1.2

for ii = 0 , n_elements(bin)-1 do oplot,bin[ii]+[-0.005,0.005],hsp[ii]*[1,1]
for ii = 1 , n_elements(bin)-1 do oplot,bin[ii]-0.005*[1,1],[hsp[ii-1],hsp[ii]]
for ii = 0 , n_elements(bin)-2 do oplot,bin[ii]+0.005*[1,1],[hsp[ii],hsp[ii+1]]
a2 = where(bin ne 0)
a1 = where(bin gt 0.1)
print,"Planck South Map - % above Ebv=0.1: ",total(hsp[a1]) * 100. / total(hsp[a2])
a1 = where(bin gt 0.3)
print,"Planck South Map - % above Ebv=0.3: ",total(hsp[a1]) * 100. / total(hsp[a2])

for ii = 0 , n_elements(bin)-1 do oplot,bin[ii]+[-0.005,0.005],hssg[ii]*[1,1],color=fsc_color('slate gray')
for ii = 1 , n_elements(bin)-1 do oplot,bin[ii]-0.005*[1,1],[hssg[ii-1],hssg[ii]],color=fsc_color('slate gray')
for ii = 0 , n_elements(bin)-2 do oplot,bin[ii]+0.005*[1,1],[hssg[ii],hssg[ii+1]],color=fsc_color('slate gray')
a2 = where(bin ne 0)
a1 = where(bin gt 0.1)
print,"Schlegel South Map - % above Ebv=0.1: ",total(hssg[a1]) * 100. / total(hssg[a2])
a1 = where(bin gt 0.3)
print,"Schlegel South Map - % above Ebv=0.3: ",total(hssg[a1]) * 100. / total(hssg[a2])


device,/close

endif 

if flag eq 'histo' then begin

;readcol, 'DUSTMAPS/euclidFoV_ebv_schlegelhealpix.cat',l,b,field,ebv,format='F,F,A,F'
;openw,1,'DUSTMAPS/euclidFoV_ebv_schlegelhealpix_histo.cat'
;nbcell = (size(l))[1]*1.

;readcol,'DUSTMAPS/euclidFoV_ebv_schlaflyhealpix.cat',l,b,field,ebv,format='F,F,A,F'
;openw,1,'DUSTMAPS/euclidFoV_ebv_schlaflyhealpix_histo.cat'
;nbcell = (size(where(ebv gt 0)))[1]*1.

readcol,'DUSTMAPS/euclidFoV_ebv_planckhealpix.cat',l,b,field,ebv,format='F,F,A,F'
openw,1,'DUSTMAPS/euclidFoV_ebv_planckhealpix_histo.cat'
nbcell = (size(where(ebv gt 0)))[1]*1.

n = where(field eq 'north' and ebv ge 0)
hn = histo(ebv[n],bin,binsize=0.01,min=0,max=2)
hn = hn / nbcell
s = where(field eq 'south' and ebv ge 0)
hs = histo(ebv[s],bin,binsize=0.01,min=0,max=2)
hs = hs /  nbcell
for ii = 0,(size(bin))[1]-1 do printf,1,bin[ii],hn[ii],hs[ii],format='(F,F,F)'
close,1

endif

if flag eq 'prep' then begin

; Schlegel previous
;lb = mrdfits('DUSTMAPS/SCHLEGEL/pixel_coords_map_nested_galactic_res9.fits',1)
;ebv = fltarr(3145728)
;openw,1,'euclidFoV_ebv_schlegelhealpix.cat'
;schlegel = mrdfits('DUSTMAPS/SCHLEGEL/lambda_sfd_ebv.fits',1,/silent)
;for ii = 0d0,3145727 do begin
;   if abs(lb[ii].LATITUDE) gt 30 then begin
;      glactc,ra,dec,2000,lb[ii].LONGITUDE,lb[ii].LATITUDE,2
;      if abs(dec) gt 15 then begin
;         ebv[ii] = schlegel[ii].TEMPERATURE
;         if lb[ii].LATITUDE gt 0 then begin
;            printf,1,lb[ii].LONGITUDE,lb[ii].LATITUDE,' north ',ebv[ii],format='(F,F,A,F)'
;         endif else begin
;            printf,1,lb[ii].LONGITUDE,lb[ii].LATITUDE,' south ',ebv[ii],format='(F,F,A,F)'
;         endelse      
;      endif      
;   endif 
;endfor
 
; SCHLEGEL
;   ebv = fltarr(3145728)
;   schlegel = mrdfits('DUSTMAPS/SCHLEGEL/lambda_sfd_ebv.fits',1,/silent)
;   openw,1,'DUSTMAPS/euclidFoV_ebv_schlegelhealpix.cat'
;   for ii = 0d0,3145727 do begin
;      pix2ang_nest,512,long64(ii),theta,phi

; PLANCK
;   ebv = fltarr(50331648)
;   planck = mrdfits('DUSTMAPS/PLANCK/PlanckEbv.fits')
;   openw,1,'DUSTMAPS/euclidFoV_ebv_planckhealpix.cat'
;   for ii = 0d0,50331647 do begin
;      pix2ang_nest,2048,long64(ii),theta,phi

; SCHLAFLY 
   ebv = fltarr(3145728) 
   schlafly= mrdfits('DUSTMAPS/SCHLAFLY/ps1-ebv-4.5kpc.fits',1,/silent)
   openw,1,'DUSTMAPS/euclidFoV_ebv_schlaflyhealpix.cat'
   for ii = 0d0,3145727 do begin
      pix2ang_nest,512,long64(ii),theta,phi
      
      l = phi * (180/!pi)
      b = 90 - (theta * (180/!pi))
      if abs(b) gt 30 then begin
         glactc,ra,dec,2000,l,b,2
         if abs(dec) gt 15 then begin
;            ebv[ii] = schlegel[ii].TEMPERATURE
;            ebv[ii] = planck[ii]
            ebv[ii] = schlafly[ii].ebv
            if b gt 0 then begin
               printf,1,l,b,' north ',ebv[ii],format='(F,F,A,F)'
            endif else begin
               printf,1,l,b,' south ',ebv[ii],format='(F,F,A,F)'
            endelse
         endif
      endif
   endfor
   
endif

; prev

;im = readfits('DUSTMAPS/PLANCK/Planck_ngp_euclid.fits',hdr)
;a=where(im ne 0)
;h=histo(im,bin,binsize=0.01,min=0,max=0.71)
;plot,bin,h,psym=3,xr=[0,2],xs=1,/ylog,yr=[1e0,1.5e7],ys=1,ytitle='# of pixels in E(B-V) map',position=[0.45,0.47,0.95,0.73],charsize=1.,xtitle='E(B-V)'
;for ii = 0 , n_elements(bin)-1 do oplot,bin[ii]+[-0.005,0.005],h[ii]*[1,1]
;for ii = 1 , n_elements(bin)-1 do oplot,bin[ii]-0.005*[1,1],[h[ii-1],h[ii]]
;for ii = 0 , n_elements(bin)-2 do oplot,bin[ii]+0.005*[1,1],[h[ii],h[ii+1]]
;xyouts,1.6,1e6,'North'
;a1 = where(im gt 0.1)
;a2 = where(im ne 0)
;print,"Planck North Map - % above Ebv=0.1: ",n_elements(a1)*1./n_elements(a2)*100.
;a1 = where(im gt 0.3)
;print,"Planck North Map - % above Ebv=0.3: ",n_elements(a1)*1./n_elements(a2)*100.


END

