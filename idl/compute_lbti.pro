;; rmg oct 2009
;; 
;; turn kin.7.pro into a function to use with KIN data fitting code
;;
;; INPUTS: 
;;
;; zodipic map of zodi cloud ONLY (no central star) (Jy).
;; wavelength (um)
;; pixel scale (mas/pixel)
;; dec (deg)
;; HA (hrs)
;; d (pc) 
;;
;; OUTPUT: the KIN leakage, specifically: the transmitted flux
;; throught the total KIN pattern with long baseline at [peak,null]. 
;; other factors to compute the measured KIN null are applied outside this function.

;; rmg jan 2013
;; adapted for LBTI

function compute_lbti, fnu, pixel_scale, wavelength, dec, ha, u, v, plot, d, lbti_null

s = pixel_scale
lambda = wavelength

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; external parameters:
;; -------------------
;; none ....?

;; read image
image = fnu
n = (size(image))[1] ; assume nx=ny

;; initiallize some arrays:
;; -----------------------
;; long baseline beam pattern at peak and null: depends on input uv
;; 
lbti_peak = dblarr(n,n)
lbti_null = dblarr(n,n)

;; image sky xy axis coordinates (mas) for display.
;; i also need them in units of radians for the computation.
x = s*(findgen(n)+0.5-n/2) ; mas
y = x
xrad = x/206264807.
yrad = y/206264807.
;plot,x,y,psym=6

;; Compute the images:
;; --------------------
;; changed the signs of x for east=left.
for i=0,n-1 do begin
   for j=0,n-1 do begin

      lbti_null[i,j] = 1./2. * (1 - cos(2.*!pi*(-xrad[i]*u + yrad[j]*v)/lambda/1.e-6))
      lbti_peak[i,j] = 1./2. * (1 + cos(2.*!pi*(-xrad[i]*u + yrad[j]*v)/lambda/1.e-6))

   endfor
endfor

;; compute KIN transmitted flux:
;; ----------------------------
;; all these images now have consistent orientation.
peak_signal = total(image*lbti_peak)
null_signal = total(image*lbti_null)
print,'LBTI 1-zodi transmitted flux at null: (Jy)',null_signal
print,'LBTI 1-zodi transmitted flux at peak: (Jy)',peak_signal
print,'==================='

if (plot) then begin

;cleanplot,/silent
;window,1,title='KIN transmission'

;white = FSC_Color("White", !D.Table_Size-4)

x2 = reverse(x,1)

str = string(format='("HA = ",f4.1)',ha)
str0 =  'Input image (East = left, North = up)'
str3 = strcompress(str + ' - LBTI long baseline fringes at Null')
str4 = strcompress(str + ' - LBTI total transmission pattern (at Null)')
str5 = strcompress(str + ' - LBTI transmission * Image (at Null)')
str6 = strcompress(str + ' - LBTI long baseline fringes at Peak')
str7 = strcompress(str + ' - LBTI total transmission pattern (at Peak)')
str8 = strcompress(str + ' - LBTI transmission * Image (at Peak)')

!p.multi = [0,3,2]

;; input image and fixed KIN patterns
;;
image_cont,image;,xval=x2,yval=y,xtit='X(mas)',ytit='Y(mas)',tit=str0,/aspect,/noc

;; patterns at null
;;
image_cont,lbti_null;,xval=x2,yval=y,xtit='X(mas)',ytit='Y(mas)',tit=str3,/aspect,/noc
image_cont,image*lbti_null;,xval=x2,yval=y,xtit='X(mas)',ytit='Y(mas)',tit=str5,/aspect,/noc
;;
;; patterns at peak
;;

;; zoom-in
;;
sub_size = 64
imin = (size(fnu))[1]/2 - sub_size/2
imax = (size(fnu))[1]/2 + sub_size/2
image_cont,image[imin:imax,imin:imax];,xval=x2[imin:imax],yval=y[imin:imax],xtit='X(mas)',ytit='Y(mas)',tit=str0,/aspect,/noc
foo = image*lbti_null
image_cont,foo[imin:imax,imin:imax];,xval=x2[imin:imax],yval=y[imin:imax],xtit='X(mas)',ytit='Y(mas)',tit=str5,/aspect,/noc
foo = image*lbti_peak
image_cont,foo[imin:imax,imin:imax];,xval=x2[imin:imax],yval=y[imin:imax],xtit='X(mas)',ytit='Y(mas)',tit=str8,/aspect,/noc
;;

!p.multi = 0

endif

return,[peak_signal,null_signal]

end
