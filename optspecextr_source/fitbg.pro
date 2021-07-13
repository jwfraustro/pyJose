;+ NAME:
;      FITBG
;
;
; PURPOSE:
;      creates a sky background image from a data image,
;      interpolated across the spectrum, to create an image of sky
;      lines only
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector Setup functions
;
;
; CALLING SEQUENCE:
;      Result = FITBG (dataim, x1, x2, bgdeg = bgdeg, bgmask = bgmask, bgres $
;                      = bgres, bpct = bpct, bthresh = bthresh, errvect = $
;                      errvect, gotovect = gotovect, inmask = inmask, nobgfit $
;                      = nobgfit, plottype = plottype, q = q, skyvar = $
;                      skyvar, varim = varim, verbose = verbose, v0 = v0)
;
;
; INPUTS:
;      DATAIM:  The data image, with spectrum
;      X1,X2:   boundaries in x which contain spectrum (spectrum is
;               inclusive of X1,X2)
;
;
; OPTIONAL KEYWORDS:
;      NOBGFIT:  set to not fit a background and assume the average
;
;      Cut into rows and passed to PROCVECT:
;      INMASK:   the main mask used by all functions
;      VARIM:    the array that contains the variance of the locations
;      SKYVAR:   the variance image for the sky subtraction
;
;      Passed Through to PROCVECT
;      BGDEG:    degree of polynomial interpolation (default is 1)
;      BTHRESH:  the threshold for sigma rejection
;      Q:        the electron per adu of the data
;      v0:       rn^2 of data
;      BPCT:     the percentage of bad pixels that kicks out of an iteration
;
;      Made from PROCVECT
;      BGRES:    the residuals for cosmic ray rejection
;      BGMASK:   the output mask of cosmic rays found
;
;      Debugging Help:
;      VERBOSE:  level of printed output
;      PLOTTYPE: a four element array with the type of plot to show
;      GOTOVECT: the row at which to stop the loop
;      ERRVECT:  will contain the rows that exited with bad pixels
;
;
; OUTPUT:
;      an array of size dataim in which, for each wavelength, the
;      spectrum from X1 to X2 has been removed, interpolated over
;      (with coefficients COEFFS) and polynomial has been evaluated at
;      all x values
;
;
; RESTRICTIONS:
;      All images must be the same length. 0 <= X1 <= X2 <= NX-1.
;      Degsky >= 0, bthresh > 0, Q > 0, v0 >= 0,
;
;
; PROCEDURE:
;      Horne Step 3
;
;
; EXAMPLE:
;      ; create a object frame using the synthetic spectrum package
;      obj = synthspec(objheader=objheader, /nosmile, /notilt)
;      rn = sxpar(objheader, 'RDNOISE')
;      gain = sxpar(objheader, 'EPADU')
;      dataim = obj
;      ; estimate X1 and X2 by finding the base of the tallest peak
;      t = total(dataim, 2)
;      mx = max(t, cnt)
;      x1 = cnt - 5 & while (t[x1-1] lt t[x1]) do x1--
;      x2 = cnt + 5 & while (t[x2+1] lt t[x2]) do x2++
;      varim = obj / gain
;      bgim = fitbg(dataim, x1, x2, varim = varim, Q = gain, v0 = rn^2, $
;                   bgmask = bgmask)
;      erase & tvscl, dataim < 1000
;      erase & tvscl, bgim
;      erase & tvscl, dataim - bgim < 1000
;
;
; FUNCTIONS CALLED:
;      procvect
;
;
; MODIFICATION HISTORY:
;       Written by:      Dara Zeehandelaar, 25 March 2003
;       8 April 2003:    finished interpolation routine
;                        added keyword COEFFS, made routine to return skyim
;       15 April 2003:   removed dependence on function concat
;       16 April 2003:   debug
;       6 May 2003:      takes mask as input, put skyproc into
;                        iteration within driver
;       20 May 2003:     changed default degree to 1 (gives best fit)
;       16 Jun 2003:     John Dermody: adding variance image keyword
;                        to allow weighted fits.  Also added basic
;                        sigma rejection scheme
;       8 Jul 2003:      JD: removed fit and rejection to its own
;                        function, and renamed file
;       19 Oct 2003:     JD: Matched header to actual calls
;       18 Dec 2003:     JD: Plots, verbosity, and a bunch of other stuff
;
;       $Log: fitbg.pro,v $
;       Revision 3.0  2004/07/08 20:07:51  jfd28
;       Moving up to 3.0 (reflects this version anyway)
;
;       Revision 2.0  2004/07/08 20:05:01  jfd28
;       Tagging all files to 2.0
;
;       Revision 1.5  2004/07/08 18:40:07  jfd28
;       plottype takes an array
;       fixed ranges on shade_surf plots
;       no recovery
;
;       Revision 1.4  2004/07/02 18:32:17  jfd28
;       removed all spec confusion.  now just assumes it is only running once
;       and doesn't try to fit the background using spectrum removal.
;       Also cleaned up header information and added a check for
;       properly formed inputs
;
;       Revision 1.3  2004/05/27 20:53:01  jfd28
;       If there is only one or no good pixels on either side then set polyfit
;       to a zero degree
;
;-

function FITBG, dataim, x1, x2, bgdeg = bgdeg, bgmask = bgmask, bgres $
                = bgres, bpct = bpct, bthresh = bthresh, errvect = $
                errvect, gotovect = gotovect, inmask = inmask, nobgfit $
                = nobgfit, plottype = plottype, q = q, skyvar = $
                skyvar, varim = varim, verbose = verbose, v0 = v0


; FILL INPUTS AND CHECK
dims = size(dataim)
nx = dims[1]
ny = dims[2]

if (verbose gt 1) then print, "Into background fitting"

if not keyword_defined(bgdeg   ) then bgdeg    = 1
if not keyword_defined(bthresh ) then bthresh  = 3.
if not keyword_defined(verbose ) then verbose  = 0
if not keyword_defined(plottype) then plottype = [0, 0, 0, 0]
if not keyword_defined(gotovect) then gotovect = -1
if not keyword_defined(inmask)   then inmask   = bytarr(nx, ny) + 1
if not keyword_defined(varim)    then varim    = dblarr(nx, ny) + 1
if not keyword_defined(skyvar)   then skyvar   = dblarr(nx, ny)

str = ""
str += x1 lt 0 || x1 gt x2                                ? "x1 "      : ""
str += x2 gt nx - 1 || x2 lt x1                           ? "x2 "      : ""
str += (size(varim ))[1] ne nx || (size(varim ))[2] ne ny ? "varim "   : ""
str += (size(inmask))[1] ne nx || (size(inmask))[2] ne ny ? "inmask "  : ""
str += (size(skyvar))[1] ne nx || (size(skyvar))[2] ne ny ? "skyvar "  : ""
str += 4 ne n_elements(plottype)                          ? "plottype" : ""
if str ne "" then message,  "Poorly formed inputs: " + str

xvals1 = dindgen(x1)            ; ignore x1:x2
xvals2 = dindgen(nx-1-x2)+x2+1  ; starts at xvalue x2+1, goes to nx-1
xvals = [xvals1,xvals2]

; Subtract off the bias so profile fitting doesn't have to worry about it
if keyword_set(nobgfit) then begin
  bgim = (dblarr(nx,ny) + 1) * median(dataim[xvals, *])
  return, bgim
endif


; PREPARE FOR ROW BY ROW FITTING
yvals  = indgen(ny)
bgim   = dataim
bgres  = dblarr(nx, ny)
bgmask = bytarr(nx, ny) + 1
func = "polyfunc"
errvect = bytarr(ny) + 1
bgim[x1:x2, *] = 0
allx = indgen(nx)


; FIT BY ROW
; cut up the data into rows and pass each to procvect.  tell procvect
; to use polyfunc.pro to estimate the background.  Procvect will
; handle bad pixel rejection, and return the polynomial over all x.
for i = 0, ny-1 do begin
  datav = dataim[*, i]          ; cut array into rows
  maskv = inmask[*, i]
  varv  = varim [*, i]
  bcrv  = bgres [*, i]
  skyvarv = skyvar[*, i]
  parm = (total(maskv[0:x1-1]) lt 2 or total(maskv[x2+1:nx-1]) lt 2) ? $
         0 : bgdeg             ; if too few good pixels then average
  if (i eq gotovect) then verbose = 5

  if plottype[1] or verbose eq 5 then begin ; plot debugging info
    device, window_state = ws   ; get window state
    if not ws[12] then window, 12 else wset, 12
    !p.multi = [0, 2, 2, 1, 1]
    plot, datav,   title = 'Data Vector for BG Fitting', /ystyle
    plot, maskv,   title = 'Input Mask', yrange = [0.0, 1.1]
    plot, varv,    title = 'Variance', /ystyle
    plot, skyvarv, title = 'Sky variance', /ystyle
    !p.multi = 0
    wait, 0.01
  endif
  if (verbose eq 5) then stop
  bgim[*, i] = procvect(datav, varv = varv, maskv = maskv, crv = bcrv, $
                        thresh = bthresh, xvals = xvals, errflag = errflag, $
                        Q = Q, v0 = v0, skyvarv = skyvarv, verbose = verbose, $
                        vectnum = i, func = func, parm = parm, bpct = bpct, $
                        plottype = plottype)
  if (errflag) then errvect[i] = 0 ; there was a problem fitting
  if plottype[3] then begin
    device, window_state = ws
    if not ws[14] then window, 14 else wset, 14
    sy1 = (i - 10) > 0
    sy2 = (i + 10) < ny-1
    shade_surf, bgim[*, sy1:sy2], allx, yvals[sy1:sy2], charsize = 3,  $
                title = "Background after # " + strtrim(i, 1), $
                xtitle = "X pixel location", $
                ytitle = "Y (Prefit and Postfit)"
    wait, 0.001
  endif
  varim[xvals, i] = varv[xvals] ; only update the variance for background
  bgmask[*, i] = maskv          ; save for user viewing
  bgres[*, i] = bcrv * maskv    ; also save
endfor

return, bgim

end
