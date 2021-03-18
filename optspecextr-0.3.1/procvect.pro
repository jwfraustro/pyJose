;+ NAME:
;      PROCVECT
;
;
; PURPOSE:
;      Fit using a user defined function to the input vector,
;      interating until no more cosmic rays are found using either
;      sigma or absolute comparison of the residuals rejection
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector fitting driver routine
;
;
; CALLING SEQUENCE:
;      Result = PROCVECT (datav, absthresh = absthresh, bgv = bgv, $
;                         bpct = bpct, coeffv = coeffv, crv = crv, errflag = $
;                         errflag, extract = extract, func = func, maskv = $
;                         maskv, multv = multv, noupdate = noupdate, opvar = $
;                         opvar, parm = parm, plottype = plottype, q = q, $
;                         skyvarv = skyvarv, thresh = thresh, v0 = v0, varv = $
;                         varv, vectnum = vectnum, verbose = verbose, xvals = $
;                         xvals)
;
;
; INPUTS:
;      DATAV:   The vector with cosmic rays to be fitted
;
;
; OPTIONAL KEYWORDS:
;      General:
;      VARV:      The variance vector
;      MULTV:     An array of same length as datav containing values to
;                 be multiplied by the fitted data in the cosmic ray rejection
;      MASKV:     cosmic ray mask (0 for cosmic ray, otherwise 1)
;      EXTRACT:   Instead of fitting the data, extract the optimal
;                 spectrum
;
;      Variance calculation:
;      BGV:       An array of wavelength dimension that contains the
;                 background values at datav
;      Q:         The gain of the array
;      V0:        The readnoise of the array
;      SKYVARV:   The sky variance vector
;      NOUPDATE:  Don't update the variance vector
;
;      Bad pixel rejection:
;      CRV:       Will contain the residual of the data fit 
;      THRESH:    sigma threshold for rejection as cosmic ray, default 5
;      BPCT:      The percentage of bad pixels that will end fitting
;      ABSTHRESH: Thresh specifies an absolute diference (not sigma)
;      ERRFLAG:   Set if a sutible fit is not found
;
;      Debugging help:
;      VERBOSE:   Set to 3 to halt the function after every iteration
;      PLOTTYPE:  Set to 3 to plot each fit itteration result
;      VECTNUM:   The vector number of the data. Sent to screen for debugging.
;
;      Fitting:
;      XVALS:     Indicies of the data to fit
;      FUNC:      The function to fit to the data.  The first five be xvals,
;                 datav, varv, multv, eval.  It should return an
;                 estimation, and in the sixth input the coefficents
;      PARM:      The additional paramaters of FUNC.
;      COEFFV:    Will contain the coeffiencents of the fit
;      OPVAR:     The optimally extracted spectrum's variance
;
;
; OUTPUTS:
;      Returns the output of func after cosmic rays have been eliminated
;
;
; SIDE EFFECTS:
;      MASKV:   Will set to 0 at the position of cosmic rays found
;      VARV:    Unless noupdate is set, will update the variance based
;               on the estimation each iteratation, and at the end for
;               all xvals
;
;
; RESTRICTIONS:
;      All vectors must be the same length
;
;
; PROCEDURE:
;      uses sigma rejection to determine positions of cosmic rays,
;      then return polynomial evaluated
;
;
; EXAMPLE:
;      ; Create a simple profile to fit a polynomial to
;      nx = 1000
;      xvals = findgen(nx)
;      profv = ((-(xvals/nx*2 - 1)^2 + 1)*.75)/nx*2
;      specl = 10000
;      gain = 5
;      basey = profv * specl
;      datav = basey + sqrt(basey*gain)*randomn(seed, nx)/gain
;      ; Initiailize variance
;      varv  = datav / gain
;      badloc = randomu(seed, 4)*nx
;      datav[badloc] = randomu(seed, 4)*max(datav)*2
;      oldvar = varv
;      maskv = datav/datav
;      crv =  datav - datav
;      est = PROCVECT (datav, maskv = maskv, varv = varv, Q = gain)
;      plot, datav & oplot, est
;      print, badloc[sort(badloc)], where(maskv eq 0)
;      plot, oldvar & oplot, varv
;
;
; MODIFICATION HISTORY:
;      Written by:    John Dermody, 25 Jun 2003
;      15 Jul 2003:   added optimally exctacting spectrum option
;      10 Aug 2003:   allows Gaussian fit
;      20 Aug 2003:   allows RSA fit
;      1 Sept 2003:   cleaning and commenting
;      19 Oct 2003:   cleaning
;      18 Dec 2003:   moved out functions into passing routines and cleaning
;
;      $Log: procvect.pro,v $
;      Revision 3.0  2004/07/08 20:07:51  jfd28
;      Moving up to 3.0 (reflects this version anyway)
;
;      Revision 2.0  2004/07/08 20:05:01  jfd28
;      Tagging all files to 2.0
;
;      Revision 1.5  2004/07/08 18:42:29  jfd28
;      added comments for stuff
;      plottype is either an array or 0
;
;      Revision 1.4  2004/07/02 18:33:39  jfd28
;      added an example
;      check inputs for proper formation
;      updated the error threshold to not reject when the input mask is mostly
;      bad pixels
;
;      Revision 1.3  2004/05/27 20:51:31  jfd28
;      updated rejection criteria
;
;-

function PROCVECT, datav, absthresh = absthresh, bgv = bgv, $
                   bpct = bpct, coeffv = coeffv, crv = crv, errflag = $
                   errflag, extract = extract, func = func, maskv = $
                   maskv, multv = multv, noupdate = noupdate, opvar = $
                   opvar, parm = parm, plottype = plottype, q = q, $
                   skyvarv = skyvarv, thresh = thresh, v0 = v0, varv = $
                   varv, vectnum = vectnum, verbose = verbose, xvals = $
                   xvals

; SET DEFAULTS AND CHECK INPUTS

nx = (size(datav))[1]           ; size of vector

if not keyword_defined(xvals   ) then xvals    = dindgen(nx)
if not keyword_defined(varv    ) then varv     = dblarr(nx) + 1.
if not keyword_defined(multv   ) then multv    = dblarr(nx) + 1.
if not keyword_defined(maskv   ) then maskv    = bytarr(nx) + 1
if not keyword_defined(bgv     ) then bgv      = dblarr(nx)
if not keyword_defined(skyvarv ) then skyvarv  = dblarr(nx)
if not keyword_defined(crv     ) then crv      = dblarr(nx)
if not keyword_defined(thresh  ) then thresh   = 5
if not keyword_defined(Q       ) then Q        = 1.
if not keyword_defined(v0      ) then v0       = 0.
if not keyword_defined(bpct    ) then bpct     = 0.5
if not keyword_defined(func    ) then func     = 'polyfunc'
if not keyword_defined(verbose ) then verbose  = 0
if not keyword_defined(plottype) then plottype = [0, 0, 0, 0]

if (n_elements(plottype) eq 1) then plottype = [0, 0, 0, 0]

str = ""
str += nx ne n_elements(varv)    ? "varv"    : ""
str += nx ne n_elements(maskv)   ? "maskv"   : ""
str += nx ne n_elements(multv) $
     && 1 ne n_elements(multv)   ? "multv"   : ""
str += nx ne n_elements(bgv)     ? "bgv"     : ""
str += nx ne n_elements(skyvarv) ? "skyvarv" : "" 
str += nx ne n_elements(crv)     ? "crv"     : ""
str += 4 ne n_elements(plottype) ? "plottype must have 4 elements" : ""
str += thresh le 0               ? "thresh less than 0" : ""
str += Q le 0                    ? "Q is less or equal to 0" : ""
str += v0 lt 0                   ? "v0 is less than 0" : ""
str += bpct gt 1 || bpct lt 0    ? "bpct is not betwee 0 and 1" : ""
if str ne "" then message, "Poorly formed inputs: " + str

; if a row or column number is given, use that for the title of the
; debugging plot
if keyword_defined(vectnum) then begin
    vectnum_s = " at Vector # " + strtrim(vectnum,1)
endif else begin
    vectnum_s = ""
endelse

; if an absolute threshold is used, use that threshold.  Else, square
; the sigma threshold so it can be used with variance calculations
vthresh = keyword_set(absthresh) ?  thresh : thresh * thresh

; set the error threshold to be the greatest of 6 pixel, 10% of the
; total pixels, or the given percentage good of pixels passed in
errorthresh = n_elements(where(maskv[xvals] eq 1)) * (1 - bpct) > $
              n_elements(xvals) * (0.10) > 6
          
errflag = 0
coeffv  = 0                     ; initialize coeffv
funcdone = 1                    ; initialize for while loop
funccount = 0                   ; number of loops

; MAIN LOOP 
; on each iteration first check to make sure there is enough good
; pixels left to fit the data.  Next pull out the good pixels and pass
; them to the function. Calculate the residuals of each pixel, using
; either the difference between actual and estimated, or if sigma
; rejection is used square the difference and divide by the variance
; of the pixel. If the user requests a summary plot, plot the
; estimated versus actual and the residual.  If needed, update the
; variance to reflect the new estimation. Reject the pixel with the
; largest residual larger then the threshold.  If no bad pixels are
; found, exit the loop.

while (funcdone eq 1) do begin
  funcdone = 0
  goodvals = where(maskv[xvals] eq 1)
  if (n_elements(goodvals) lt (errorthresh)) then begin
    fiteval = call_function(func, dindgen(nx), datav, varv, $
                            multv*maskv, 1, coeffv, parm)
    if (verbose gt 2) then print, "Too many pixels rejected" + vectnum_s
    errflag = 1
    return, fiteval
  endif
  fitx = xvals[goodvals]        ; xvals for good pixels
  fitdata = datav[fitx]         ; data for good pixels
  fitvar = varv[fitx]           ; variance of good pixels
  fitmult = multv[fitx]         ; multiplier for good pixels
  est = call_function(func, fitx, fitdata, fitvar, fitmult, 0, coeffv, parm)
  if keyword_set(absthresh) then begin
    crv[fitx] =  abs(fitdata / fitmult - est)
  endif else begin 
    crv[fitx] = (fitdata - fitmult * est)^2 / (fitvar > 1.e-6) ; no divide by 0
  endelse
  badpix = where (crv[fitx] gt vthresh) ; get bad locations
  if plottype[2] or (verbose eq 5) then begin
    device, window_state = ws
    if not ws[13] then window, 13 else wset, 13
    !p.multi = [0,1,2,1,1]
    plot, fitx, fitdata, $
          title='Actual vs. Fitted' + vectnum_s, $
          ytitle='Data Values (if applicable / Spec)', $
          xtitle='Pixel Locations', charsize=1.1
    oplot, fitx, est * fitmult
    if keyword_set(absthresh) then begin
      plot, fitx, crv[fitx], title='Residual', $
            ytitle='Abs of Data - Expected', $
            xtitle='Pixel Locations', charsize=1.1, $
            yrange = [0, max([thresh, max(crv[fitx])])]
      oplot, fitx, fitdata/fitdata*thresh
    endif else begin
      plot, fitx, sqrt(crv[fitx]), title='Residuals', $
            ytitle='Sigma Difference of Data vs. Expected', $
            xtitle='Pixel Locations', charsize=1.1
    endelse
    !p.multi = 0
    wait, 0.01
  endif
  if (verbose eq 5) then stop
  if not keyword_set(noupdate) then begin
    varv[fitx] = (abs(fitmult * est + bgv[fitx])) / Q + v0 + skyvarv[fitx]
  endif
  if (badpix ne [-1]) then begin
    badx = fitx[badpix]
    maxpos = where(crv[badx] eq  max(crv[badx])) ; only elimate max pixel
    maxx = badx[maxpos]
    funccount = funccount + n_elements(maxx) ; add count for bad pix
    maskv[maxx] = 0             ; mask bad pix
    funcdone = 1                ; set so sequence loops
  endif
endwhile

fiteval = call_function(func, dindgen(nx), datav, varv, $
                        multv*maskv, 1, coeffv, parm)

if not keyword_set(noupdate) then begin
  varv = (abs(multv * fiteval + bgv)) / Q + v0 + skyvarv 
                                ; get var for all pixels
end

return, fiteval

end
