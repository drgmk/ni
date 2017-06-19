;; $Id$

;
; NAME:
;         LBTI_ZODI_LIMITS
;
; PURPOSE:
;         Calculate simple Zodi limits for an LBTI observation
;
; CATEGORY:
;         LBTI
;
; EXPLANATION:
;         Calculate Zodi limits using the simple disk model, given an LBTI leak or upper
;         limit. Result is a series of plots
;
; SYNTAX:
;         LBTI_ZODI_LIMITS,name,leak,e_leak,dist,lstar,fstar,dec,hamin,hamax,
;                          pa=pa,inc=inc,no=no,nr=nr,zs=zs
;
; INPUTS:
;            name: name of source for file names
;            leak: observed leak
;          e_leak: 1sigma uncertainty on leak
;            dist: distance to star in pc
;           lstar: stellar luminosity in L_sun
;           fstar: stellar flux density in Jy
;             dec: Declination of star (degrees)
;           hamin: hour angle at start of obs, can be an array for many obs (in hours)
;           hamax: hour angle at end of obs, can be an array for many obs (in hours)
;
; KEYWORD PARAMETERS:
;
;              pa: Position angle of disk E of N on sky (degrees)
;             inc: Disk inclination (from face-on, degrees)
;              no: number of random orientations
;              nr: number of radii to use in disk model
;             r11: photometric flux ratio at 11um
;
; OUTPUTS:
;              zs: array of zodi limits
;
; EXAMPLE:
;
;         Make plots for the first LBTI detection of an exo-Zodi around eta Corvi
;
;         lbti_zodi_limits,'eta Crv',0.044,0.0035,18.28,5.2,1.59,-16.1958,[-0.28,0.52,1.97],[0.02,0.92,2.17],inc=46.8,pa=116.3,zs=zs,nr=100,no=5000

;
; MODIFICATION HISTORY:
;         Made in 2014 by GMK
;
;-

pro lbti_zodi_limits,name,leak,e_leak,dist,lstar,fstar,dec,hamin,hamax,pa=pa,inc=inc,r11=r11,no=no,nr=nr,zs=zs

  ;; sanitise object name
  fn = STRJOIN( STRSPLIT(name,' ',/EXTRACT), '-' )
  thk = 6
  csz = 1.4
  !P.FONT=0
  set_plot,'ps'

  ;; u,v coverage for the observation
  sz = 10
  ntx = lbti_observe(0,0,0,0,hamin,hamax,dec,uvu=u,uvv=v,/UVONLY,has=has)
  device,file=fn+'-uv.eps',/ENC,/TIMES,XSIZE=8.2,YSIZE=8,/INCHES
  plot,u,v,psym=1,/iso,xr=[sz,-sz],yr=[-sz,sz],xtitle='u (m)',ytitle='v (m)',$
       xthick=thk,ythick=thk,charthick=thk,charsize=csz+0.6,position=[0.14,0.12,0.99,0.99]
  a0 = [-1,-1]*9
  al = 5
  arrow,a0[0],a0[1],a0[0]+al,a0[1],/data,thick=thk,hsize=500
  arrow,a0[0],a0[1],a0[0],a0[1]+al,/data,thick=thk,hsize=500
  xyouts,a0[0],a0[1]+al+0.2,'N',alignment=0.5,charsize=csz
  xyouts,a0[0]+al+0.1,a0[1]-0.2,'E',alignment=1,charsize=csz
  if N_ELEMENTS(pa) ne 0 then begin
     oplot,[sz,-sz],[sz,-sz]/TAN(pa*!dtor),linestyle=2,thick=thk ; disk PA
     xyouts,0,sz/50.,'Outer disk PA',orientation=pa-90,alignment=0.5,charsize=csz+0.6,charthick=thk
;     al_legend,['Outer disk PA'],/bot,/lef,box=0,charsize=csz+0.6,charthick=thk,thick=thk,linestyle=2,margin=0
  endif
  for i=0,N_ELEMENTS(hamin) -1 do push,lstr,'HA: '+STRING(hamin[i],FOR='(F4.1)')+' to '+STRING(hamax[i],FOR='(F4.1)')+'h'
  al_legend,/top,/lef,box=0,[name,lstr],charsize=csz+0.6,margin=0,charthick=thk
  device,/CLOSE

  ;; now model with many orientations
  if N_ELEMENTS(no) eq 0 then no = 1000                ; no of orientations
  if N_ELEMENTS(nr) eq 0 then nr = 2000l               ; no of radii
  fd = lbti_disk(rarcs=rd,dist=dist,lstar=lstar,nr=nr) ; disk model fluxes and radii

  ;; make plot showing how model leaks vary with HA
  device,file=fn+'-ha.eps',/ENC,/TIMES,XSIZE=8,YSIZE=6,/INCHES
  xr = MINMAX([hamin,hamax])
  plot,[0],[0],xrange=xr,yrange=100*[0,2*leak],xtitle='HA (h)',ytitle='Calibrated leak (%)',xthick=thk,ythick=thk,charthick=thk,charsize=csz,position=[0.12,0.12,0.99,0.99],/xst,/yst

  txs = FLTARR(no)
  incs = FLTARR(no)
  pas = FLTARR(no)
  txsnoha = txs                 ; also do it for no movement in HA (at mean HA)
  for i=0,no-1 do begin
     incs[i] = acos(randomu(sd))/!dtor ; more high inclinations than low
     pas[i]  = randomu(sd)*180.0       ; could be 90deg due to symmetry

     ;; transmission at this orientation for true range of HA
     txs[i] = lbti_observe(fd,rd,pas[i],incs[i],hamin,hamax,dec,has=has,hanulltx=hatx)

     ;; add these lines to plot
     if i lt 200 then begin
        oplot,has,100*hatx*leak/MEAN(hatx),color=100,psym=6,symsize=0.1,thick=thk ; show how null tx varies with HA
        oplot,has,100*hatx*leak/MEAN(hatx),color=100                              ; show how null tx varies with HA
     endif

     ;; also calculate at mean HA, so can show improvement due to HA coverage
     txsnoha[i] = lbti_observe(fd,rd,pas[i],incs[i],mean(has)-0.01,mean(has)+0.01,dec,nha=2)
     counter,i+1,no,'Orientations done: '                                             
  endfor                                                                              
  zs = leak / (txs/FLOAT(fstar))         ; limits in zodi units
  zsnoha = leak / (txsnoha/FLOAT(fstar)) ; limits in zodi units

  ;; now if given PA/Inc from elsewhere
  if N_ELEMENTS(pa) eq 1 then begin
     tx1 = lbti_observe(fd,rd,pa,inc,hamin,hamax,dec,has=has,hanulltx=hatx)
     z1 = leak/(tx1/FLOAT(fstar))
     pstr = 'At PA='+STRING(pa,FOR='(I3)')+', I='+STRING(inc,FOR='(I2)')+': '+STRING(z1,FOR='(I4)')+textoidl('\pm')+STRING(e_leak/leak*z1,FOR='(I3)')+'z'
     oplot,has,100*hatx*leak/MEAN(hatx),thick=thk,color=FSC_COLOR('dark red')
  endif else pstr = ''

  ;; add measurements
  oplot,xr,100*[1,1]*leak,linestyle=2,thick=thk
  oplot,xr,100*[1,1]*(leak-e_leak),linestyle=1,thick=thk
  oplot,xr,100*[1,1]*(leak+e_leak),linestyle=1,thick=thk
  device,/CLOSE

  ;; figure range to plot for histogram
  zs1 = zs[SORT(zs)]
  maxin = zs1[no*0.97]
  binsz = maxin/20.0

  device,file=fn+'-zlim.eps',/ENC,/TIMES,XSIZE=8,YSIZE=6,/INCHES
  histoplot,zs,/FILLPOLYGON,/FREQUENCY,thick=thk,xthick=thk,ythick=thk,charthick=thk,charsize=csz,xtitle='Zodi',ytitle='Fraction of orientations',maxinput=maxin,binsize=binsz,yrange=[0,1],position=[0.12,0.12,0.99,0.99],ytickformat='(F3.1)'
  
;  histoplot,zsnoha,/oplot,/freq,binsize=binsz ; show histogram for no range of HA

  ;; figure stats
  zf = [0.05,0.16,0.5,0.84,0.95] ; fractions we care about
  zq = zs1[no*zf]                 ; CI

  ezstat = e_leak/leak*zq[2]                      ; stats error on meadian
  zmin = zq[2] - sqrt((zq[2]-zq[1])^2 + ezstat^2) ; quad added orientation and stat errors
  zmax = zq[2] + sqrt((zq[3]-zq[2])^2 + ezstat^2) ; quad added orientation and stat errors

  ;; show cumulative distribution
  cumnum,zs,/GTX,x=x,y=y
  y /= FLOAT(no)                ; turn into fraction
  oplot,x,y,thick=thk,color=150
  for i=0,N_ELEMENTS(zf)-1 do begin
     oplot,[zq[i]],[1-zf[i]],thick=thk,psym=1
     xyouts,zq[i],(1-zf[i])+0.01,STRN(100*(1-zf[i]),FORMAT='(I)')+'%',alignment=0.5
  endfor

  al_legend,/top,/rig,box=0,charsize=csz,[name,'',lstr,'','Calibrated null: '+STRING(leak*100,FOR='(F4.1)')+textoidl('\pm')+STRING(e_leak*100,FOR='(F4.2)')+'%','Zodi range: '+STRING(zmin,FOR='(I4)')+'<'+STRING(zq[2],FOR='(I4)')+'<'+STRING(zmax,FOR='(I4)')+'z',pstr]
  
  ;; add median and uncertainty
  b1 = 0.91*!Y.CRANGE[1]
  ba = 0.93*!Y.CRANGE[1]
  b2 = 0.95*!Y.CRANGE[1]
  oplot,[zq[2]],[ba],thick=thk,psym=1,symsize=1.2
  oplot,[1,1]*(leak-e_leak)*zq[2]/leak,[b1,b2],thick=thk
  oplot,[1,1]*(leak+e_leak)*zq[2]/leak,[b1,b2],thick=thk
  oplot,[leak-e_leak,leak+e_leak]*zq[2]/leak,[1,1]*ba,thick=thk
  xyouts,(leak+e_leak)*zq[2]/leak,(b1+ba)/2.,textoidl('median \pm null uncertainty'),alignment=-0.05,charsize=csz;-0.3

  device,/CLOSE
 
  set_plot,'x'

  save,file=fn+'.xdr'
;stop
end

  ;; figure median zodi limit for varying disk width about EEID
;;   ndr = 5
;;   drs = linarrw(min=0.1,max=2*r0,n=ndr); to twice r0 for this star
;;   zqm = FLTARR(3,ndr)
;;   for j=0,ndr-1 do begin
;;      fd = lbti_disk(dr=drs[j],nr=nr,rarcs=rd,dist=dist,lstar=lstar)
;;      txs = FLTARR(no)
;;      for i=0,no-1 do begin
;;         incr = acos(randomu(sd))/!dtor                        ; more high inclinations than low
;;         par  = randomu(sd)*180.0                              ; could be 90deg due to symmetry
;;         txs[i] = lbti_observe(fd,rd,par,incr,hamin,hamax,dec) ; 
;;      endfor                                                   ;
;;      zstmp = leak / (txs/fstar)                               ; limits in zodi units
;;      zstmp = zstmp[SORT(zstmp)]
;;      zqm[*,j] = zstmp[no*[0.25,0.5,0.75]]
;;      counter,j+1,ndr,'Radial widths done: ' ;
;;   endfor

;;   device,file=fn+'-zlim-dr.eps',/ENC,/TIMES,XSIZE=8,YSIZE=6,/INCHES
;;   plot,drs,zqm[1,*],thick=thk,charthick=thk,xthick=thk,ythick=thk,charsize=csz,xtitle='Radial width (AU)',ytitle='Median and IQR Zodi',/xst,/yl
;;   oplot,drs,zqm[0,*],linestyle=2,thick=thk
;;   oplot,drs,zqm[2,*],linestyle=2,thick=thk
;;   al_legend,/top,/rig,box=0,[name,'','Disk centered on EEID at '+STRING(r0,FOR='(F3.1)')+'AU','',lstr],charsize=csz,charthick=thk
;;   device,/CLOSE

;;   set_plot,'x'

  ;; plot stats, box/whiskers styles
  ;; b1 = 0.92*!Y.CRANGE[1]
  ;; ba = 0.94*!Y.CRANGE[1]
  ;; b2 = 0.96*!Y.CRANGE[1]
  ;; oplot,[1,1]*zq[0],[b1,b2],thick=thk,color=00
  ;; oplot,[1,1]*zq[2],[b1,b2],thick=thk,color=00
  ;; oplot,[1,1]*zq[1],[b1,b2],thick=thk,color=00
  ;; oplot,[zq[1],zq[3]],[1,1]*b1,thick=thk,color=00
  ;; oplot,[zq[1],zq[3]],[1,1]*b2,thick=thk,color=00
  ;; oplot,[1,1]*zq[3],[b1,b2],thick=thk,color=00
  ;; oplot,[1,1]*zq[4],[b1,b2],thick=thk,color=00
  ;; oplot,[zq[0],zq[1]],[ba,ba],thick=thk,color=00
  ;; oplot,[zq[3],zq[4]],[ba,ba],thick=thk,color=00
  ;; xyouts,zq[4],(b1+ba)/2.,'z>5,16,50,84,95%',alignment=-0.05,charsize=csz-0.3

;; end
