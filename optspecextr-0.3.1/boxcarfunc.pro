;+
; NAME:
;      BOXCARFUNC
;
;
; PURPOSE:
;      Function for fitting the data with a boxcar average or median fit.
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Program
;        - Vector fitting functions
;
;
; CALLING SEQUENCE:
;      Result = BOXCARFUNC (XVALS, DATAV, VARV, SPECV, COEFFV, BOXCARHW)
;
;
; INPUTS:
;      XVALS:    The x values for the data
;      DATAV:    The vector with cosmic rays to be fitted
;      VARV:     The variance vector
;      SPECV:    The spectrum vector.  Bad pixels have a specv of 0
;      EVAL:     Set to 1 to evaluate the inputs
;      COEFFV:   A filler in this routine
;      BOXCARHW: The halfwidth of the median or smooth filter
;
;
; OUTPUTS:
;      Returns an estimation of the data by returning datav ran
;      through a median filter.  If eval is set the routine smooths
;      over the good pixels, then interpolates the bad pixels. 
;
;
; RESTRICTIONS:
;      All vectors must be the same length.  Boxcarhw must be an int >= 1
;
;
; PROCEDURE:
;      If eval eq 0 then do a median filter over the datav / specv.  Else,
;      smooth over the good pixels and interpolate over the bad pixels.  
;
;
; EXAMPLE:
;      ; Creates demo vectors using simple functions and randomness
;      noise = 2
;      nx = 50
;      xvals = findgen(nx)
;      basey = (-(xvals - nx/2)^2/(nx/2) + nx)
;      smult = 10
;      specv = sin(xvals/50*6*!pi)*2 + smult
;      datav = basey*specv/smult + randomn(seed, nx)*noise
;      varv  = noise^2 + randomn(seed, nx)/4
;      badloc = randomu(seed, 4)*50
;      datav[badloc] = randomu(seed, 4)*100
;      eval = 0
;      boxcarhw = 3
;      ; Run without knowledge of bad pixels
;      est = boxcarfunc(xvals, datav, varv, specv, eval, coeffv, boxcarhw)
;      plot, datav / specv
;      oplot, est
;      specv[badloc] = 0
;      eval = 1
;      ; Evaluate with knowledge of bad pixels
;      est = boxcarfunc(xvals, datav,  varv, specv, eval, coeffv, boxcarhw)
;      plot, datav / specv
;      oplot, est
;
;
; MODIFICATION HISTORY:
;       Written by: John Dermody, 10 Oct 2003
;
;       Thu Dec 18 22:53:47 2003, John Dermody
;		cleaned up and documented.  
;
;       $Log: boxcarfunc.pro,v $
;       Revision 3.0  2004/07/08 20:07:51  jfd28
;       Moving up to 3.0 (reflects this version anyway)
;
;       Revision 2.0  2004/07/08 20:05:01  jfd28
;       Tagging all files to 2.0
;
;       Revision 1.6  2004/07/08 18:43:37  jfd28
;       added comments
;
;       Revision 1.5  2004/07/02 19:03:55  jfd28
;       added an example to show how the function works
;       check for properly formed inputs
;
;       Revision 1.4  2004/06/03 16:33:57  jfd28
;       changed print in error catch to a message command
;
;       Revision 1.3  2004/06/03 15:32:38  jfd28
;       Cleaned up header. Set category to vector fitting functions
;       Added extensive example (will be used with all fitting functions
;       Removed divide by zero
;
;-

function BOXCARFUNC, xvals, datav, varv, specv, eval, coeffv, boxcarhw


; CHECK INPUTS
nx = n_elements(xvals)

str = ""
str += nx ne n_elements(datav) ? "datav " : ""
str += nx ne n_elements(varv ) ? "varv "  : ""
str += nx ne n_elements(specv) $
     && 1 ne n_elements(specv) ? "specv " : "" 
str += boxcarhw*2+1 ge nx      ? "boxcarhw to large for datav " : ""
str += boxcarhw lt 0           ? "boxcarhw less than 0 " : ""
if str ne "" then message, "Poorly formed inputs: " + str

coeffv = [mean(datav/specv)]

; FIT DATA
if (eval eq 0) then begin
  return, median(datav/specv, boxcarhw*2+1) ; median seems to work the best
endif

; EVALUATE
gv = where(specv ne 0, complement = bv) ; good and bad pixel locations
if (gv eq [-1]) then return, dblarr(nx) ; a row of zeros of right length
goodx = xvals[gv]
estg = smooth(datav[gv]/specv[gv], boxcarhw*2+1, /edge_truncate)
fiteval = fltarr(nx)
if (bv ne [-1]) then begin      ; interpolate the bad from good pixels
  badx = xvals[bv] 
  estb = interpol(estg, goodx, badx)
  fiteval[badx ] = estb
endif
fiteval[goodx] = estg

return, fiteval

end
