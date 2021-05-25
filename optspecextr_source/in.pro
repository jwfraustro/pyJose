;+
; NAME:
;      IN
;
;
; PURPOSE:
;      Checks to see if an element is in an array
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;            - Miscellaneous
;
;
; CALLING SEQUENCE:
;      Result = IN(el, ar)
;
;
; INPUTS:
;      EL: an element
;      AR: an array
;
;
; OUTPUTS:
;      True or false
;
;
; MODIFICATION HISTORY:
;      Written by John Dermody 08 Jul 2004
; 
;      $Log: in.pro,v $
;      Revision 3.0  2004/07/08 20:07:51  jfd28
;      Moving up to 3.0 (reflects this version anyway)
;
;      Revision 2.0  2004/07/08 20:05:01  jfd28
;      Tagging all files to 2.0
;
;      Revision 1.1  2004/07/08 19:22:43  jfd28
;      initial CVS
;
;
;-

function in, el, ar
nel = n_elements(el)
case nel of
  0: return, false
  1: return, not (where(el eq ar) eq [-1])
  else: begin
    for i = 0, nel-1 do begin
      if where(el[i] eq ar) eq [-1] then return, false
    endfor
    return, true
  endelse
endcase

end

