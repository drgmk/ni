;; $Id: bnuw.pro 114 2012-08-29 16:12:00Z gkennedy $

;+
;
; NAME:
;         BNUW
;
; PURPOSE:
;         A function that returns Bnu in Jy/sr for a given temperature
;         and at one or more given wavelengths (in um)
;
; CATEGORY:
;         Optprops
;
; SYNTAX:
;         Return = bnuw(wavelengths,temperature)
;
; INPUTS:
;       wav: array of wavelengths in microns
;         t: temperature
;
; MODIFICATION HISTORY:
;         Made by Mark sometime...
;         Modified by Grant in March 2010 to make floating overflow values zero.
;
;-

function bnuw,wav,t

  if (n_params() eq 0) then begin
     print,'Correct usage is ... = bnuw(wav,t),'
     print,'  where wav is in um, t in K and the'
     print,'  output is in Jy/sr'
  endif

  ;; Give the constants for this function, where:
  ;; h = 6.6260755*10^{-34} Js is Planck's constant,
  ;; k = 1.380658*10^{-23} J/K is Boltzmann's constant,
  ;; c = 2.997924580*10^8 m/s is the speed of light.
  
  k1 = 3.9728949e19             ;**** = 2hc in Jy*mum^3, i.e., = (2hc)*1e(18+26)
  k2 = 14387.69                 ;**** = hc/k in mum*K, i.e., = (hc/k)*1e6

  wav = float(wav)
  fact1 = k1/wav^3
  fact2 = k2/(wav*t)

  ;; Avoid the message "Floating Overflow" if fact2 > 88.0
  ret = fact1/(exp(fact2 < 88.0) - 1.0)

  ;; now set values where fact2 would have been > 88 to zero. not doing this can lead to
  ;; problems if you integrate the curve, since extra flux turns up where bnuw should
  ;; have been zero but isn't if fact2 <= 88
  tmp = where(fact2 gt 88.0)
  if tmp[0] ne -1 then ret(tmp) = 0.0

  return,ret

end
