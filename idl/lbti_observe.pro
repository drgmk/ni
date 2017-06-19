;; $Id$

;
; NAME:
;         LBTI_OBSERVE
;
; PURPOSE:
;         Calculate transmission through LBTI for a real observation
;
; CATEGORY:
;         LBTI
;
; EXPLANATION:
;         Calculate transmission of an azimuthally integrated radial disk profile though
;         LBTI at null and peak for a real observation, specifically accounting for the
;         rotation of the disk on the sky during the integration.
;
;         Uses some equations grabbed from http://www.ucolick.org/~sla/deimos/swpdr/va.html
;
;         Good summary for fixed telescopes (e.g. Keck Nuller) here
;         http://www.vlti.org/events/assets/1/proceedings/1.4_Segransan.pdf
;
; SYNTAX:
;         result = LBTI_OBSERVE(fnudisk,rarcs,pa,inc,hamin,hamax,dec,wav=wav
;                              peaktx=peaktx,rtxnull=rtxnull,rtxpeak=rtxpeak)
;
; INPUTS:
;         fnudisk: Disk radial profile, azimuthally integrated (i.e. SB * 2pi x r X dr)
;           rarcs: Radial bin locations (centers) of disk profile in arc seconds
;              pa: Position angle of disk E of N on sky (degrees)
;             inc: Disk inclination (from face-on, degrees)
;           hamin: hour angle at start of obs, can be an array for many obs (in hours)
;           hamax: hour angle at end of obs, can be an array for many obs (in hours)
;             dec: declination of source (degrees)
;
; KEYWORD PARAMETERS:
;
;             wav: Wavelength to compute transmission at in micron, default=11um
;           hadeg: set to give hour angle in degrees
;             nha: number of hour angle steps to use (default=100)
;
; OUTPUTS:
;             alt: source altitude at hour angles (degrees)
;             has: array of hour angles (hours)
;              va: PA E of N of local vertical (degrees)
;          peaktx: Total transmission at peak
;        nulltxha: null transmission at hour angle positions
;        peaktxha: peak transmission at hour angle positions
;
; MODIFICATION HISTORY:
;         Made in 2013 by GMK
;
;-

function lbti_observe,fnudisk,rarcs,pa,inc,hamin,hamax,dec,wav=wav,hadeg=hadeg,nha=nha,$
                      alt=alt,has=hah,va=va,peaktx=peaktx,hanulltx=nulltxha,hapeaktx=peaktxha,$
                      uvu=u,uvv=v,uvonly=uvonly,elev=h

  ;; get star declination in radians
  delta = dec * !dtor

  ;; LBTI details that we won't change until the next nuller comes along, at which
  ;; point we will probably change both.
  phi = 32.7 * !dtor            ; latitude
  dlbti = 14.4                  ; baseline in m

  ;; sort out arrays of hour angles
  if N_ELEMENTS(hamin) gt 1 then begin
     nobs = N_ELEMENTS(hamin)
     if nobs ne N_ELEMENTS(hamax) then STOP,'LBTI_OBSERVE: hamin/max not same size'
     dha = hamax - hamin        ; observation lengths
     delvarx,alt,hah,va,nulltx,peaktx,nulltxha,peaktxha,u,v,h
     for i=0,nobs-1 do begin

        if ARG_PRESENT(peaktx) then begin
           nulltx1 = lbti_observe(fnudisk,rarcs,pa,inc,hamin[i],hamax[i],dec,wav=wav,hadeg=hadeg,nha=nha,$
                                  alt=alt1,has=hah1,va=va1,peaktx=peaktx1,hanulltx=nulltxha1,$
                                  hapeaktx=peaktxha1,uvu=u1,uvv=v1,uvonly=uvonly,elev=h1)
        endif else begin
           nulltx1 = lbti_observe(fnudisk,rarcs,pa,inc,hamin[i],hamax[i],dec,wav=wav,hadeg=hadeg,nha=nha,$
                                  alt=alt1,has=hah1,va=va1,hanulltx=nulltxha1,$
                                  uvu=u1,uvv=v1,uvonly=uvonly,elev=h1)
        endelse
        push,alt,alt1
        push,hah,hah1
        push,va,va1
        if N_ELEMENTS(nulltx1) ne 0 then push,nulltx,nulltx1
        if N_ELEMENTS(peaktx1) ne 0 then push,peaktx,peaktx1
        if N_ELEMENTS(nulltxha1) ne 0 then push,nulltxha,nulltxha1
        if N_ELEMENTS(peaktxha1) ne 0 then push,peaktxha,peaktxha1
        push,u,u1
        push,v,v1
        push,h,h1
     endfor
     ;; average according to time spent along each HA track
     if N_ELEMENTS(peaktx) ne 0 then peaktx = TOTAL(dha*peaktx)/TOTAL(dha)
     if N_ELEMENTS(nulltx) ne 0 then return,TOTAL(dha*nulltx)/TOTAL(dha) else return,!VALUES.F_NAN
  endif

  ;; make an array of HAs to observe at, these are evenly spaced in time
  if N_ELEMENTS(nha) eq 0 then nha = 20 else nha = FIX(nha)
  has = hamin + (hamax-hamin)*findgen(nha)/(nha-1)
  if ~ KEYWORD_SET(hadeg) then has *= 360/24.0 ; convert to degrees unless given in degrees
  has *= !dtor                                 ; convert to radians
  hah = has * 24/360.0/!dtor                   ; HA in hours for returning

  ;; figure the position angles of these, LBTIs transmission pattern is always constant
  ;; in the local vertical (i.e. HorPA=0)
  h = asin(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos(has))                  ; altitude in radians
  cosva = (sin(phi)*cos(delta) - cos(phi)*sin(delta)*cos(has))/sqrt(1-sin(h)^2) ; cos(VA)
  sinva = sin(has) * cos(phi) / sqrt(1 - sin(h)^2)                              ; sin(VA)
  va = atan(sinva,cosva)/!dtor                                                  ; VA in degrees
  alt = h / !dtor                                                               ; altitude

  ;; figure u-v coverage, which depends on baseline and vertical angle. positive u is to
  ;; the East
  u = dlbti/2. * sin( (90+va)*!dtor) ; LBTI baseline is perpendicular to the vertical angle
  v = dlbti/2. * cos( (90+va)*!dtor) ; ditto
  u = [-u,u]
  v = [-v,v]
  if KEYWORD_SET(uvonly) then return,!VALUES.F_NAN

  ;; now get the transmission at each of these times
  nulltxha = FLTARR(nha)
  peaktxha = FLTARR(nha)
  for i=0,nha-1 do begin
     palbti = va[i] - pa        ; PA of disk relative to LBTI fringes (see LBTI_TRANSX)
     if ARG_PRESENT(peaktx) then begin
        nulltxha[i] = LBTI_TRANSX(fnudisk,rarcs,palbti,inc,wav=wav,peaktx=peaktx)
        peaktxha[i] = peaktx
     endif else begin
        nulltxha[i] = LBTI_TRANSX(fnudisk,rarcs,palbti,inc,wav=wav) ; save some time if peaktx not needed
     endelse
  endfor

  ;; and average them
  nulltx = MEAN(nulltxha)
  if ARG_PRESENT(peaktx) then peaktx = MEAN(peaktxha)

  return,nulltx

end

;; http://www.ucolick.org/~sla/deimos/swpdr/va.html

;; Finding the Vertical Angle (VA)

;; The declination of an object is delta, with north being positive. The latitude of an
;; observer is phi, with north being positive. The Hour Angle of an object is HA; zero on
;; the meridian, negative to the east, and positive to the west. In the equatorial
;; (RA,Dec) coordinate system the position angle of a directed vector on the sky is
;; measured from equatorial North toward equatorial East; this is EquPA. In the Horizon
;; (Alt-Az) coordinate system the position angle of a directed vector on the sky is
;; measured from Up toward an upright observer's Left; this is HorPA. The elevation of an
;; object is h, and cos(h) is always positive. The Vertical Angle (VA) is defined as the
;; PA of the local vertical as expressed in the equatorial system.

;; VA = EquPA - HorPA

;; (Note that in many applications this angle is mis-named the ``Parallactic
;; Angle''. Lately this has been promoted by terminology found in the extremely useful
;; SLALIB software by Pat Wallace, but the misnomer can be traced back at least as far as
;; Spherical Astronomy by Smart. The correct usage of Parallactic Angle is in the context
;; of the Equatorial and Ecliptic coordinate systems. In that context it describes the
;; orientation of the ellipse that stellar coordinates traverse due to annual
;; parallax. The term ``Parallactic Angle'' should not be applied in the context of the
;; Horizon system.)

;; Spherical trigonometry gives

;; sin(h) = sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(HA) 
;; cos(VA) = (sin(phi)*cos(delta) - cos(phi)*sin(delta)*cos(HA)) / sqrt(1 - sin(h)**2) 
;; sin(VA) = sin(HA) * cos(phi) / sqrt(1 -sin(h)**2)
