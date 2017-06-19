;; $Id: cumnum.pro 256 2010-03-26 18:02:53Z gkennedy $

;+
;
; NAME:
;         Find the cumulative number of something
;
; CATEGORY:
;         Utilities
;
; SYNTAX:
;         cumnum,sizes,x=x,y=y,/GTX
;         plot,x,y,/xlog,/ylog
;
; PARAMETERS:
;         gtx   find cumulative number greater than x instead of less than
;
; MODIFICATION HISTORY:
;         Made by Grant in March 2010
;
;-

pro cumnum,arr,x=x,y=y,GTX=gtx,REVERSE_INDICES=srt

  ;; make x array
  srt = SORT(arr)
  x = arr[srt]

  ;; make y array
  n = N_ELEMENTS(arr)
  if KEYWORD_SET(gtx) then y = REVERSE(LINDGEN(n)+1) else y = LINDGEN(n)+1

end
