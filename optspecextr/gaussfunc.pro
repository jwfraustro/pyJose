; NAME: GAUSSEVAL
; PURPOSE: Small helper routine
;
; simple function that evaluates a Gaussian of coeffecents a over
; values x, and returns in f

pro gausseval, x, a, f, p

z = (x - a[1]) / a[2]

f = a[0]*exp(-z^2 / 2)

if n_params(0) le 3 then return

p = fltarr(n_elements(x), 3)

p[*, 0] = f/a[0]
p[*, 1] = f * z / a[2]
p[*, 2] = p[*, 1] * z

end

;+
; NAME:
;      GAUSSFUNC
;
;
; PURPOSE:
;      Function for estimating a Gaussian curve to the data
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector fitting functions
;
;
; CALLING SEQUENCE:
;      Result = GAUSSFUNC (XVALS, DATAV, VARV, SPECV, COEFF, REEST)
;
;
; INPUTS:
;      XVALS:  The x values for the data
;      DATAV:  The vector with cosmic rays to be fitted
;      VARV:   The variance vector, pixels not to fit have a 0
;      SPECV:  The spectrum at that point, bad pixels have a 0
;      EVAL:   Do not fit, instead, evaluate based on coeff
;      COEFF:  Will contain the coefficients of the fit
;      REEST:  Set to reestimate the spectrum after bad pixel rejection
;
;
; KEYWORDS:
;      None
;
;
; OUTPUTS
;      Returns an estimation of the data by fitting a Gaussian.  If a
;      good fit is not found quickly to the data, the Gaussian is
;      fitted to the data smoothed instead.  If eval is set it simply
;      evaluates the Gaussian coeffeciants passed in at all points.
;
;
; RESTRICTIONS:
;      All vectors must be the same length
;
;
; EXAMPLE:
;      mult = 2
;      nx = 100*mult
;      xvals = findgen(nx)
;      profv = ((-(xvals/nx*mult*2 - mult)^2 + 1)*.75)/nx*mult*2
;      profv[0:nx/2 - nx/mult/2] = 0 & profv[nx/2+nx/mult/2:nx-1] = 0
;      specl = 1000
;      gain = 10
;      rn = 1
;      basey = profv*specl
;      datav = basey + basey*randomn(seed, nx)/gain + randomn(seed, nx)*rn
;      varv  = basey/gain + rn^2
;      badloc = randomu(seed, 4)*nx
;      datav[badloc] = randomu(seed, 4)*max(datav)*2
;      profv[badloc] = 0
;      goodloc = where(profv ne 0)
;      eval = 0
;      reest = 0
;      est = GAUSSFUNC (xvals, datav, varv, specl, eval, coeffv, reest)
;      plot, est & oplot, datav/specl
;      eval = 1 
;      est = GAUSSFUNC (xvals, datav, varv, specv, eval, coeffv, reest)
;      plot, est
;
;
; MODIFICATION HISTORY:
;      Written by:    John Dermody, 22 Sep 2003
;
;      Thu Dec 18 22:55:25 EST 2003, John Dermody
;            Cleaned up and docummented
;
;      $Log: gaussfunc.pro,v $
;      Revision 3.0  2004/07/08 20:07:51  jfd28
;      Moving up to 3.0 (reflects this version anyway)
;
;      Revision 2.0  2004/07/08 20:05:01  jfd28
;      Tagging all files to 2.0
;
;      Revision 1.4  2004/07/08 18:44:21  jfd28
;      merged with gausseval
;      gausseval now returns partial derivates
;      fit a median if the first fit fails insteadof smoothed data
;
;      Revision 1.3  2004/06/04 19:40:38  jfd28
;      added example
;      cleaned documentation
;      removed normalization from this routine
;      check inputs for proper form
;
;-

function GAUSSFUNC, xvals, datav, varv, specv, eval, coeff, reest

; CHECK INPUTS
nx = n_elements(xvals)

str = ""
str += nx ne n_elements(datav) ? "datav " : ""
str += nx ne n_elements(varv ) ? "varv "  : ""
str += nx ne n_elements(specv) $
     && 1 ne n_elements(specv) ? "specv " : "" 
str += nx le 4                 ? "datav has less than 4 elements" : ""
if str ne "" then message, "Poorly formed inputs: " + str


; EVALUATE COEFFICIENTS 
if (eval eq 1) then begin
  gausseval, xvals, coeff, ret
  return, ret
endif


; PREPARE FITTING

if keyword_set(reest) then begin
  specv = int_tabulated(xvals, datav) ; reestimate spectrum
endif

hbw = 4 < nx/2-1                ; half box width for estimation
sdata = smooth(datav / specv, 2*hbw+1, /edge_truncate) ; est Gauss
shi = max(sdata, chi)
slo = min(sdata, cli)
ci = abs(shi) gt abs(slo) ? chi : cli
hi = (datav/specv > sdata)[ci]  ; estimated Gaussian height
center = xvals[ci]              ; estimated Gaussian center
top = where(abs(datav/specv) gt abs(hi/exp(1)))
fwhm = max(xvals[top]) - min(xvals[top]) + 1 
                                ; estimated Gaussian FullWidthHalfMax
base = median(sdata)            ; estimated base level
coeff = [hi, center, fwhm/2]    ; coefficiance for estimate

nz = where(varv ne 0)           ; locations to fit

est = datav / specv
origcoeff = coeff

; FIT DATA

;est[nz] = gaussfit(xvals[nz], datav[nz]/specv, coeff, $
;                     measure_errors = sqrt(varv[nz]), nterms = 4, $
;                     estimates = origcoeff, chisq = chisq)

merrors = specv ^ 2 / varv > 1e-8
estz = curvefit(xvals, datav/specv, merrors, $
                   coeff, iter = fititer, function_name = 'gausseval', $
                   tol = 1e-6, /double, status = status, $
                   itmax = 10, chisq = chisq)
est[nz] = estz[nz]

if (status ne 0) then begin     ; estimation wrong so use median instead
  mdata = median(datav/specv, 2*hbw+1)
  estz = curvefit(xvals, mdata, merrors, $
                     origcoeff, iter=fititer, function_name='gausseval', $
                     /double, status = status) 
  est[nz] = estz[nz]
endif 

return, est
end
