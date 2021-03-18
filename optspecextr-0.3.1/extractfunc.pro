;+
; NAME:
;      EXTRACTFUNC
;
;
; PURPOSE:
;      A function which extracts the optimal spectrum from a vector
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector fitting functions
;
;
; CALLING SEQUENCE:
;      Result = EXTRACTFUNC (XVALS, DATAV, VARV, PROFV, EVAL, COEFFV, OPVAR)
;
;
; INPUTS:
;      XVALS:  The x values for the data (not used)
;      DATAV:  The sky subtracted data vector
;      VARV:   The variance vector
;      PROFV:  The profile vector.  0 at bad pixel locations.
;      EVAL:   (Not used)
;      COEFFV: (Not used)
;
;
; OUTPUTS:
;      OPVAR: The variance of the optimal spectrum
;      Returns the optimal extraction
;
;
; RESTRICTIONS:
;      All vectors must be the same length
;
;
; EXAMPLE:
;      nx = 100
;      xvals = findgen(nx)
;      profv = ((-(xvals/nx*2 - 1)^2 + 1)*.75)/nx*2
;      specl = 10000
;      gain = 10
;      basey = profv*specl
;      datav = basey + basey*randomn(seed, nx)/gain
;      varv  = basey/gain
;      badloc = randomu(seed, 4)*50
;      datav[badloc] = randomu(seed, 4)*max(datav)*2
;      profv[badloc] = 0
;      goodloc = where(profv ne 0)
;      optspec = EXTRACTFUNC (xvals, datav, varv, profv, foo, foo, opvar)
;      print, "Actual, Standard with Bad, Standard no Bad, Optimal"
;      print, specl, total(datav), total(datav[goodloc]), optspec
;
;
; MODIFICATION HISTORY:
;      Written by:    John Dermody, 22 Sep 2003
;
;      $Log: extractfunc.pro,v $
;      Revision 3.0  2004/07/08 20:07:51  jfd28
;      Moving up to 3.0 (reflects this version anyway)
;
;      Revision 2.0  2004/07/08 20:05:01  jfd28
;      Tagging all files to 2.0
;
;      Revision 1.3  2004/07/08 18:45:58  jfd28
;      added comments
;      added CVS
;
;-

function extractfunc, xvals, datav, varv, profv, eval, coeffv, opvar

; CHECK INPUTS

nx = n_elements(datav)

str = ""
str += nx ne n_elements(profv) ? "profv " : ""
str += nx ne n_elements(varv)  ? "varv "  : ""
if str ne "" then message, "Poorly formed inputs: " + str


; ALWAYS EXTRACT

gl = where(profv ne 0.)
if (gl eq [-1]) then begin
  print, "No good pixels in datav"
  return, 0
endif
denom = total( (profv[gl] * profv[gl]) / varv[gl]) ; avoid recalc
opt   = total( (profv[gl] * datav[gl]) / varv[gl]) / denom
opvar = total(  profv[gl] ) / denom

return, opt
end
