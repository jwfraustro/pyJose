;+
; NAME:
;      POLYFUNC
;
;
; PURPOSE:
;      Function for fiting the data using a polynomial fit
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector fitting functions
;
;
; CALLING SEQUENCE:
;      Result = POLYFUNC (XVals, Datav, Varv, Specv, Eval, Coeffv, Deg)
;
;
; INPUTS:
;      XVals:  The indicies of the datav
;      Datav:  The data vector to fit, must be same size as xvals
;      Varv:   The variance vector
;      Specv:  The spectrum vector
;      Deg:    The degree of fit
;
;
; OUTPUTS:
;      The estimated spectrum
;
;
; OPTIONAL OUTPUTS:
;      Coeffv:  Will contain the coefficients of the fit
;
;
; PROCEDURE:
;      If eval then evaluate the coefficents passed in. Else, fit the
;      data / spectrum using the variances as weights using a
;      polynomial function.  If deg is 0 then use mean, if deg is 1
;      then use linfit.
;
;
; EXAMPLE:
;      nx = 100
;      xvals = findgen(nx)
;      profv = ((-(xvals/nx*2 - 1)^2 + 1)*.75)/nx*2
;      sbump = 10
;      specv = sin(xvals/nx * !pi * 6) + sbump
;      gain = 10
;      basey = profv*specv
;      datav = basey + basey*randomn(seed, nx)/gain
;      varv  = basey/gain
;      badloc = randomu(seed, 4)*50
;      datav[badloc] = randomu(seed, 4)*max(datav)*2
;      profv[badloc] = 0
;      goodloc = where(profv ne 0)
;      eval = 0
;      deg = 2
;      est = POLYFUNC (xvals, datav, varv, specv, eval, coeffv, deg)
;      plot, est & oplot, datav/specv
;      eval = 1 
;      est = POLYFUNC (xvals, datav, varv, specv, eval, coeffv, deg)
;      plot, est
;
;
; MODIFICATION HISTORY:
;
;       Thu Dec 18 22:58:02 2003, John Dermody
;
;		If eval is set then evaluate Coeff.  Cleaned and Documented.
;
;       Sun Oct 19 17:26:22 2003, John Dermody
;
;		Created.  Check for locations inwhich variance equals
;		0, then use polyfit to evaluate other points.
;
;       $Log: polyfunc.pro,v $
;       Revision 3.0  2004/07/08 20:07:51  jfd28
;       Moving up to 3.0 (reflects this version anyway)
;
;       Revision 2.0  2004/07/08 20:05:01  jfd28
;       Tagging all files to 2.0
;
;       Revision 1.9  2004/07/08 18:43:24  jfd28
;       added comments
;
;       Revision 1.8  2004/07/02 18:34:18  jfd28
;       added an example to show how the function works
;       check inputs for proper formation
;
;       Revision 1.7  2004/06/04 19:06:44  jfd28
;       cleaned up header
;       added example
;       corrected evaluation techique when variance is 0.
;
;       Revision 1.6  2004/05/27 20:50:24  jfd28
;       ok done testing now
;
;       Revision 1.5  2004/05/27 20:48:17  jfd28
;       typo 2
;
;       Revision 1.4  2004/05/27 20:46:00  jfd28
;       fixed typo
;
;       Revision 1.3  2004/05/27 20:44:46  jfd28
;       deg eq 0 does a mean over all values, deg eq 1 does a linfit
;
;-

function POLYFUNC, xvals, datav, varv, specv, eval, coeffv, deg

; CHECK INPUTS
nx = n_elements(xvals)
deg = n_elements(deg) ne 0 ? deg : 2

str = ""
str += nx ne n_elements(datav) ? "datav " : ""
str += nx ne n_elements(varv ) ? "varv "  : ""
str += nx ne n_elements(specv) ? "specv " : ""
str += deg lt 0 || nx le deg   ? "deg "   : ""
if str ne "" then message, "Poorly formed inputs: " + str


; EVALUATE COEFFICIENTS
if (eval eq 1) then begin       ; evaluate the given coefficients
  fiteval = polyeval(coeffv, xvals)
  zl = where(varv eq 0.)        ; use actual data at varv 0 locations
  if (zl ne [-1]) then fiteval[zl] = (datav[zl] / specv[zl])
  return, fiteval
endif


; FIT DATA
nz = where(varv ne 0.)          ; locations where variance is zero
est = datav / specv             ; initial estimate
if (nz ne [-1]) then begin
  if (deg eq 0) then begin      ; use a average over the column
    mn = mean(datav/specv)
    est[nz] = mn
    coeffv = [mn, 0]            ; correct for right length
  endif else begin
    merrors = varv/specv^2 > 1e-8 ; use all pixels for estimation
    if (deg eq 1) then begin
      coeffv = linfit(xvals, datav/specv, $
                      measure_errors = merrors, yfit = estz)
    endif else begin
      coeffv = poly_fit(xvals, datav/specv, deg, /double, $
                        yfit = estz, yband = yband, measure_errors = merrors)
    endelse
    est[nz] = estz[nz]
  endelse 
endif

return, est
end
