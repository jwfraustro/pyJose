;+
; NAME:
;      FINDSHIFT
;
;
; PURPOSE:
;      Find the shift in the trace for each row
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Shifting Procedures
;
;
; CALLING SEQUENCE:
;      Result = FINDSHIFT (dataim, centerdeg = centerdeg, $
;                          centerfit = centerfit, centroid = centroid, $
;                          gaussth = gaussth, height = height, inmask = $
;                          inmask, plottype = plottype, shiftth = shiftth, $
;                          shmask = shmask, spec = spec, tracemask = $
;                          tracemask, varim = varim, verbose = verbose, $
;                          widthth = widthth, widthv = widthv, $
;                          x1 = x1, x2 = x2)
;
;
; INPUTS:
;      DATAIM: sky-subtracted processed image, horizontal is pixel, 
;              vertical is wavelength
;
;
; KEYWORD PARAMETERS:
;      General
;      VARIM:      the variance image for weighting of the poly_fit
;      INMASK:     the mask used for all functions
;      SPEC:       array holding extracted spectra, same dimension as
;                  vertical of datatim
;      X1,X2:      boundaries in x which contain spectrum (spectrum is
;                  inclusive of X1,X2)
;      VERBOSE:    level of output to screen (0)
;      PLOTTYPE:   The type of plot to output while running (0)
;      
;      Center and width finding
;      WIDTHV:     The width of the profile
;      HEIGHT:     The height of the profile
;      TRACEMASK:  The mask of bad values found durring center estimation
;      CENTROID:   Set to use centroid fitting instead
;      GAUSSTH:    The absolute value of allowable error when fitting a
;                  gaussian to the data (3 percent)
;
;      Center and width correction
;      CENTERFIT:  Set to fit the estimated centers and not accept the
;                  estimation to be correct
;      CENTERDEG:  Degree of fit for estimating shift (2)
;      SHIFTTH:    The absolute allowable error between a row's center
;                  and the fitted trace (0.5 pixels)
;      WIDTHTH:    The absolute allowable error between a row's width
;                  and the average width
;      SHMASK:     A mask with 0 for each row whose center appears to
;                  be miscalculated
;
;
; OUTPUTS:
;      A vector of each wavelength's shift in the trace at that wavelength
;
;
; RESTRICTIONS:
;      Dataim, varim and inmask must be the same size.  0 <= X1 <=
;      X2 <= NX-1.  Spec must be the same length as NY. Level >= 1.
;      Tracedeg >= 0. Gaussth, shiftth > 0.
;
;
; PROCEDURE:
;      Loop over all wavelengths finding the center of mass or fitting
;      a Gaussian and saving the center, marking bad fits. Polyfit the
;      vector of centers, and throw away those with more than shiftth
;      difference, return the poly fitted trace id tracefit, or 
;      return the origional estimated trace.
;
;
; EXAMPLE:
;      ; create a object and sky frame using the synthetic spectrum package
;      obj = synthspec(sky, objheader=objheader, skyheader=skyheader)
;      rn = sxpar(objheader, 'RDNOISE')
;      gain = sxpar(objheader, 'EPADU')
;      dataim = obj - sky
;      ; estimate X1 and X2 by finding the base of the tallest peak
;      t = total(dataim, 2)
;      mx = max(t, cnt)
;      x1 = cnt - 5 & while (t[x1-1] lt t[x1]) do x1--
;      x2 = cnt + 5 & while (t[x2+1] lt t[x2]) do x2++
;      varim = obj / gain + rn^2 + 2*sky / gain + rn^2
;      ; standard extraction of spectrum
;      spec = total(obj[x1:x2, *], 1)
;      shift = findshift(dataim, varim=varim, spec=spec, x1=x1, x2=x2)
;
;
;
; MODIFICATION HISTORY:
;      Written By John Dermody, Cornell Jun 05 2004
;
;      $Log: findshift.pro,v $
;      Revision 3.1  2004/07/13 15:39:20  jfd28
;      fixed bug so that it checks for centerfit not tracefit
;      added gotovect feature
;
;      Revision 3.0  2004/07/08 20:07:51  jfd28
;      Moving up to 3.0 (reflects this version anyway)
;
;      Revision 2.0  2004/07/08 20:05:01  jfd28
;      Tagging all files to 2.0
;
;      Revision 1.4  2004/07/08 18:42:14  jfd28
;      fixed centrioding
;      plottype is an array
;      added comments
;
;      Revision 1.3  2004/07/02 18:28:15  jfd28
;      use CENTERFIT instead of NOTRACEFIT for simplification
;
;      Revision 1.2  2004/06/24 15:00:36  jfd28
;      pulled varim, inmask, spec, x1, and x2 into keywords
;
;      Revision 1.1  2004/06/23 15:39:33  jfd28
;      Taken from fitprof.
;
;      Uses absolute threshholds for Gaussian fitting and trace fitting
;      output the width and height from the Gaussian as well
;
;      Allows center of mass calculation, then throws out rows with bad pixels
;
;-

function FINDSHIFT, dataim, centerdeg = centerdeg, $
                    centerfit = centerfit, centroid = centroid, $
                    gaussth = gaussth, height = height, inmask = $
                    inmask, plottype = plottype, shiftth = shiftth, $
                    shmask = shmask, spec = spec, tracemask = $
                    tracemask, varim = varim, verbose = verbose, $
                    widthth = widthth, widthv = widthv, x1 = x1, x2 = x2


; SET DEFAULTS AND CHECK INPUTS
nx = (size(dataim))[1]
ny = (size(dataim))[2]

; Defaults
if not keyword_defined(gaussth ) then gaussth  = 0.03
if not keyword_defined(shiftth ) then shiftth  = 0.5
if not keyword_defined(widthth ) then widthth  = 0.5
if not keyword_defined(tracedeg) then tracedeg = 4
if not keyword_defined(verbose ) then verbose  = 0
if not keyword_defined(plottype) then plottype = [0, 0, 0, 0]
if not keyword_defined(varim   ) then varim    = abs(dataim)
if not keyword_defined(inmask  ) then inmask   = bytarr(nx, ny)+1
if not keyword_defined(spec    ) then spec     = total(dataim, 1)
if not keyword_defined(x1      ) then x1       = 0
if not keyword_defined(x2      ) then x2       = nx-1
if not keyword_defined(gotovect) then gotovect = -1

; Checks
str = ""
str += n_elements(spec) ne ny                             ? "spec "   : ""
str += x1 lt 0 || x1 gt x2                                ? "x1 "     : ""
str += x2 gt nx - 1 || x2 lt x1                           ? "x2 "     : ""
str += (size(varim ))[1] ne nx || (size(varim ))[2] ne ny ? "varim "  : ""
str += (size(inmask))[1] ne nx || (size(inmask))[2] ne ny ? "inmask " : ""
str += 4 ne n_elements(plottype)                          ? "plottype" : ""
if str ne "" then message, "Poorly formed inputs: " + str

; Initialize
shiftv = dblarr(ny)             ; the center of each row's profile
widthv = dblarr(ny)             ; the width of each row's profile
height = dblarr(ny)             ; the heigth of each row's profile
xrange = indgen(x2-x1+1) + x1   ; the range of pixel values to examine
shmask = bytarr(ny) + 1         ; mask of bad estimates of profile
tracemask = inmask              ; the mask of rejected pixels

saveverbose = verbose
for i = 0,  ny-1 do begin
  mainmv = inmask[x1:x2, i]
  datav = dataim[x1:x2, i] > 0
  varv =  varim [x1:x2, i] * mainmv
  specv = spec[i] 
  if (i eq gotovect) then verbose = 5
  if keyword_set (centroid) then begin
    func = 'centermass'
    shiftv[i] = centermass(xrange, datav, varv, specv)
    widthv[i] = 1
    if (verbose eq 5) then stop
  endif else begin
    if plottype[1] or (verbose eq 5) then begin
      device, window_state = ws
      if not ws[12] then window, 12 else wset, 12
      !p.multi = [0, 1, 2, 1, 1]
      plot, datav,   title = 'Data Vector for Finding Center', /ystyle
      plot, mainmv,   title = 'Input Mask', yrange = [-0.1, 1.1]
      !p.multi = 0
      wait, 0.01
    endif

    if (verbose eq 5) then stop

    func = "gaussfunc"
    parm = 1
    foo = procvect(datav, varv = varv, /noupdate, thresh = gaussth, $
                   maskv = mainmv, func = func, parm = parm, $
                   verbose = verbose, plottype = plottype, coeffv = coeff, $
                   multv = specv, vectnum = i, /absthresh)
    tracemask[x1:x2, i] = mainmv
    ;foo = gaussfunc(xrange, datav, varv, specv, 0, coeff)
    if (coeff[0] eq 0) then begin 
      shmask[i] = 0
    endif else begin
      shmask[i] = 1
      shiftv[i] = coeff[1]
      widthv[i] = coeff[2]
      height[i] = coeff[0]
    endelse
  endelse
endfor

varv = dblarr(ny) + (1./10)^2
func = "polyfunc"
parm = tracedeg

if keyword_set(centroid) then begin 
  brows = where(total(inmask[xrange, *], 1) lt n_elements(xrange), rcount)
  if rcount ne 0 then shmask[brows] = 0
endif

if plottype[1] or verbose eq 5 then begin
  device, window_state = ws
  if not ws[12] then window, 12 else wset, 12
  !p.multi = [0, 1, 2, 1, 1]
  plot, shiftv,  title = 'Shift Vector', /ystyle
  plot, shmask,   title = 'Input Mask', yrange = [-0.1, 1.1]
  !p.multi = 0
  wait, 0.5
endif

verbose = saveverbose
if (verbose eq 5) then stop

shiftest = procvect(shiftv, varv = varv, /noupdate, thresh = shiftth, $
                    maskv = shmask, func = func, parm = parm, $
                    verbose = verbose, plottype = plottype, $
                    /absthresh)

if keyword_set(centerfit) then begin
  shiftfinal = shiftest
endif else begin
  shiftfinal = shiftv * shmask + shiftest * (1 - shmask)
endelse

if keyword_defined(widthth) then begin
  badw = where(abs(widthv - median(widthv)) gt widthth, wcount)
  if wcount gt 0 then widthv[badw] = median(widthv)
endif

if plottype[1] or verbose eq 5 then begin
  device, window_state = ws
  if not ws[12] then window, 12 else wset, 12
  plot, shiftv, title = 'Shift in PSF vs Estimate', xtitle = 'Wavelength', $
        ytitle = 'Shift amount in pixels'
  oplot, shiftest
  wait, 1.5
endif

return, shiftfinal
end
