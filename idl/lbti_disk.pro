;; $Id$

;
; NAME:
;         LBTI_DISK
;
; PURPOSE:
;         Simple parameterised disk model
;
; CATEGORY:
;         LBTI
;
; EXPLANATION:
;         Generate radial flux and surface brightness profiles for a simple power-law
;         axisymmetric optically-thin disk. The default model parameters are inspired by
;         the Kelsall COBE/DIRBE model, using their surface density and radial profile.
;
; SYNTAX:
;         result = LBTI_DISK(sp=sp,lstar=lstar,tstar=tstar,dist=dist,
;                            r0=r0,rin=rin,rout=rout,alpha=alpha,z=z,nr=nr
;                            wav=wav,
;                            bigsig=bigsig,tbb=tbb,rmas=rmas,fnustar=fnustar,snudisk=snudisk)
;
; INPUTS:
;         spty: spectral type, used to get stellar parameters via Schmidt-Kaler table
;        lstar: Stellar luminosity in L_sun, default=1 (overrides spty)
;        tstar: Stellar effective temperature, default=5800 (overrides spty)
;         dist: Distance to star in pc, default=10
;
;           r0: Disk reference radius in AU, default=278.3K (e.g. propto SQRT(Lstar))
;          rin: Disk innner radius, default=1500K
;         rout: Disk outer radius, default=10AU x SQRT(Lstar)
;           dr: Disk width when centered on rcen (alternative to using rin,rout)
;         rcen: Disk center (alternative to using rin,rout)
;        alpha: Disk surface density power-law index, default=0.34 (i.e. Kelsall)
;        zodis: Number of zodis in disk, simple multiplier
;           nr: Number of radial locations in disk
;
;          wav: Wavelength of interest, default=11um
;
; OUTPUTS:
;         Returns an array with the disk flux at each radial location r. The flux is that
;         integrated around in azimuth (i.e the surface brightness x 2pi x rarcs x dr)
;
; KEYWORD PARAMETERS:
;
;       bigsig: Disk surface density profile at rarcs (in AU^2/AU^2)
;      snudisk: Disk surface brightness (in Jy/arcsec^2)
;          tbb: Disk temperature at r
;        rarcs: Disk radial locations (bin centers) in arcseconds
;      fnustar: Stellar flux density at wavelength wav (using blackbody model)
;
; MODIFICATION HISTORY:
;         Made in 2013 by GMK
;
;-

function lbti_disk, $

   ;; star input
   spty=spty,lstar=lstar,tstar=tstar,dist=dist, $
   ;; disk input
   r0=r0,rin=rin,rout=rout,dr=dr,rcen=rcen,alpha=alpha,zodis=z,nr=nr, $
   ;; Some observation input
   wav=wav, $
   ;; some knobs
   zodipictemp=zodipictemp, $
   ;; output
   bigsig=bigsig,tbb=tbb,rarcs=rarcs,fnustar=fnustar,snudisk=snudisk

  if n_elements(spty) eq 1 then starparams,spty,ls=lstar,te=tstar         ; use SpType to get L*,T* if given
  if n_elements(lstar) ne 1 then lstar = 1.0    else lstar = FLOAT(lstar) ; in Lsun
  if n_elements(tstar) ne 1 then tstar = 5800.0 else tstar = FLOAT(tstar) ; in K
  if n_elements(dist) ne 1 then   dist = 10.0   else dist = FLOAT(dist)   ; in pc
  if n_elements(alpha) ne 1 then alpha = 0.34   else alpha = FLOAT(alpha) ;
  if n_elements(z) ne 1 then         z = 1.0    else z = FLOAT(z)         ; in zodis at r0
  if n_elements(wav) ne 1 then     wav = 11.0   else wav = FLOAT(wav)     ; in micron
  
  ;; star
  fnustar = 1.77*bnuw(wav,tstar)*lstar*tstar^(-4)/FLOAT(dist)^2

  ;; Define inner edge at 1500K
  if N_ELEMENTS(rin) ne 1 then rin = (278.3/1500.0)^2*sqrt(lstar) else rin = FLOAT(rin) ; in AU

  ;; disk T and r and dr, accounting for stellar luminosity
  if N_ELEMENTS(rout) eq 0 then rout = 10*SQRT(lstar) else rout = FLOAT(rout) ; in AU

  ;; Define a zodi, doing with temperature means that r0 varies with sqrt(Lstar)
  if N_ELEMENTS(r0) eq 0 then r0 = sqrt(lstar) else r0 = FLOAT(r0)
  
  ;; or use dr if given
  if N_ELEMENTS(dr) ne 0 then begin
     if N_ELEMENTS(rcen) eq 0 then rcen = r0
     rin = rcen - dr/2.0
     rout = rcen + dr/2.0
  endif

  ;; set up disk
  if (n_elements(nr) ne 1) then nr = 2000l
  nrb     = nr + 1
  rb      = rin + (rout-rin)*findgen(nrb)/(nrb-1)
  rm      = 0.5*(rb(1:nr)+rb(0:nr-1))
  drm     = rb(1:nr) - rb(0:nr-1)
  rarcs   = rm/dist
  tbb     = 278.3*lstar^(0.25)/sqrt(rm)
  if KEYWORD_SET(zodipictemp) then tbb = 286.0 * rm^(-0.467)* (lstar^(0.467/2.0))

  ;; 1 zodi surface density by definition after integrating Kelsall model g(xi) over
  ;; column to get multiplier of 0.629991 for 1.13e-7 volume density
  sigzodi = 7.11889e-8

  ;; disk flux density
  bigsig   = sigzodi*z*(rm/r0)^(-alpha)                   ; AU^2/AU^2 in each of the bins
  aream    = 2*!pi*rm*drm                                 ; Area in AU^2 (of the bin)
  sig      = aream*bigsig                                 ; Total dust area in each of the rm bins in AU^2
  fnudisk  = 2.95e-10*bnuw(wav,tbb)*(0.25*sig/!pi)/dist^2 ; Jy in each of the bins
  snudisk  = fnudisk/(aream/dist^2)                       ; surface brightness (Jy/arcsec^2)
  
  return,fnudisk
  
end
