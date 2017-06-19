;; $Id$

;
; NAME:
;         LBTI_TRANSX
;
; PURPOSE:
;         Calculate transmission through LBTI
;
; CATEGORY:
;         LBTI
;
; EXPLANATION:
;         Calculate transmission of an azimuthally integrated radial disk profile though
;         LBTI at null and peak. The returned value is the transmission at null
;
;         TOTAL(fnudisk*transxnull/peak)
;
;         the transmission at the peak is also calculated as both are needed to calculate
;         the leak (see Millan-Gabet et al 2012). The trig terms used here are a bit
;         different to Millan-Gabet et al (2012) but are equivalent.
;
;         The position angle is the disk position angle relative to the LBTI projected
;         finges, so a disk with PA=0deg will be parallel to the transmission pattern and
;         if highly inclined may not be detected.
;
; SYNTAX:
;         result = LBTI_TRANSX(fnudisk,rarcs,pa,inc,wav=wav,
;                              peaktx=peaktx,rtxnull=rtxnull,rtxpeak=rtxpeak,KECK=keck)
;
; INPUTS:
;         fnudisk: Disk radial profile, azimuthally integrated (i.e. SB * 2pi x r X dr)
;           rarcs: Radial bin locations (centers) of disk profile in arc seconds
;              pa: Position angle of disk relative to LBTI fringe pattern
;             inc: Disk inclination (from face-on)
;
; KEYWORD PARAMETERS:
;             wav: Wavelength to compute transmission at in micron, default=10um
;            KECK: set baseline and transmission for Keck Nuller
;
; OUTPUTS:
;          peaktx: Total transmission at peak
;         rtxnull: Average radial transmission at null
;         rtxpeak: Average radial transmission at peak
;
; MODIFICATION HISTORY:
;         Made by GMK in 2013
;         Added PSF term, 15 July 2014 (took it out, 17 July...)
;
;-

function lbti_transx,fnudisk,rarcs,pa,inc,wav=wav,peaktx=peaktx,rtxnull=rtxnull,rtxpeak=rtxpeak,KECK=keck

  ;; transmission function
  if N_ELEMENTS(wav) ne 1 then wav = 11.0 else wav = FLOAT(wav) ;
  blbti = 14.4                                                  ; mirror spacing (m)
  if KEYWORD_SET(keck) then begin                               ;
     blbti = 85.0                                               ; for Keck Nuller
     dlbti = 4.0                                                ; size for PSF attenuation (~half Keck aperture)
  endif                                                         ;
  null1as = (648000.0/!pi)*0.5*wav*1e-6/blbti                   ; Location of the first null in arcsec

  ;; compute transmission for this orientation
  nth = 1000
  theta = findgen(nth)*2*!pi/FLOAT(nth) ; 100 points around azimuth
  nr = N_ELEMENTS(rarcs)
  rtxnull = FLTARR(nr)
  rtxpeak = FLTARR(nr)

  ;; projection of a unit radius point in the disk onto the sky plane
  xrt = sin(pa*!dtor)*cos(theta) + cos(pa*!dtor)*sin(theta)*cos(inc*!dtor)

  ;; transmitted flux for a disk 'wedge' is (flux/no wedges)*trans and total is
  ;; Sum((flux[i]/no)*trans[i]), which is the same as Sum(flux[i]*(trans[i]/no))
;  for i=0,nr-1 do rtxnull[i] = mean(sin(0.5*!pi*rarcs[i]*xrt/null1as)^2)

  ;; or make a matrix of r,theta to avoid for loop
  xrt1 = REBIN(xrt,nth,nr)                    ; xrt along first dimension
  rarcs1 = TRANSPOSE(REBIN(rarcs,nr,nth))     ; radius along second dimension
  rtxall = sin(0.5*!pi*rarcs1*xrt1/null1as)^2 ; "fringe" term
  rtxnull = MEAN(rtxall,DIMENSION=1)          ; average over thetas (first dimension)

  ;; if we're actually doing Keck we'll include this...
  if KEYWORD_SET(keck) then begin
     u = !pi*dlbti*(!pi*rarcs/648000.0)/(wav*1e-6) ;
     psf = ( 2*BESELJ(u,1)/u )^2                   ; PSF term for transmission
     rtxnull *= psf
  endif

  if ARG_PRESENT(peaktx) then begin
;     for i=0,nr-1 do rtxpeak[i] = mean(cos(0.5*!pi*rarcs[i]*xrt/null1as)^2)
     rtxall = cos(0.5*!pi*rarcs1*xrt1/null1as)^2
     rtxpeak = MEAN(rtxall,DIMENSION=1)
     if KEYWORD_SET(keck) then rtxpeak *= psf
     peaktx = TOTAL(fnudisk*rtxpeak)
  endif

  return,TOTAL(fnudisk*rtxnull)

end
