;+ NAME:
;      EXTRSPEC
;
;
; PURPOSE:
;      Optimally extract spectra using weighted profiles
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector setup functions
;
;
; CALLING SEQUENCE:
;      function EXTRSPEC (dataim, profim, varim, v0, Q, x1, x2, $
;                         bgim = bgim, bpct = bpct, errvect = errvect, $
;                         ethresh = ethresh, exmask = exmask, exres = exres, $
;                         gotovect = gotovect, inmask = inmask, opvar = $
;                         opvar, plottype = plottype, skyvar = skyvar, $
;                         verbose = verbose)
;
;
; INPUTS:
;      DATAIM: sky-subtracted, processed image to extract spectrum
;              from.  Horizontal is pixel position, vertical is
;              wavelength, assumes no curvature.
;      PROFIM: image of spatial profiles
;      VARIM:  variance image from processed image
;      V0:     root(v0) is squared readout noise in DN
;      Q:      effective number of photons per DN
;      X1,X2:  boundaries in x which contain spectrum (spectrum is
;              inclusive of X1,X2)
;
;
; OPTIONAL KEYWORDS:
;      Sent to Procvect:
;      BGIM:   sky background image, interpolated across the spectrum
;      INMASK: input bad pixel mask for this array (default is all 1's)
;      SKYVAR: variance image of sky frame
;      ETHRESH: sigma threshold for rejection as cosmic ray, (5)
;      BPCT:   percentage of bad pixels that halt extraction
;
;      Returned from Procvect:
;      OPVAR:    The optimally extracted spectrum's variance
;      EXRES:    holds dataim-spec*profim*profim/varim
;      EXMASK:   the mask of cosmic rays found in this routine
;      ERRVECT:  the mask with 0's meaning the data was not able to be
;                extracted at that row
;
;      Debugging helpers:
;      GOTOVECT: Stop at that row (-1)
;      PLOTTYPE: The type of plot to output to screen ([0,0,0,0])
;      VERBOSE:  the level of output to the screen   (0)
;
;
; OUTPUTS:
;      returns fopt, an array of optimally extracted spectra, and their
;      variances
;
;
; RESTRICTIONS:
;      Dataim, profim, varim, bgim, inmask and skyvar must be the same
;      size.  0 <= X1 <= X2 <= NX-1.  V0 >= 0, Q > 0.  0 < BPCT < 1.
;      Ethresh > 0
;
;
; PROCEDURE:
;      Horne Step 8
;      Cut the images into rows and optimally extract the spectrum from each,
;      using procvect.  
;
;
; FUNCTIONS CALLED:
;      procvect
;
;
; EXAMPLE:
;      ; create a object frame using the synthetic spectrum package
;      obj = synthspec(sky, objheader=objheader, proframe = profim)
;      rn = sxpar(objheader, 'RDNOISE')
;      gain = sxpar(objheader, 'EPADU')
;      dataim = obj - sky
;      ; estimate X1 and X2 by finding the base of the tallest peak
;      t = total(dataim, 2)
;      mx = max(t, cnt)
;      x1 = cnt - 5 & while (t[x1-1] lt t[x1]) do x1--
;      x2 = cnt + 5 & while (t[x2+1] lt t[x2]) do x2++
;      varim = obj / gain + rn^2 + sky / gain + rn^2
;      optspec = extrspec(dataim, profim, varim, v0, Q, x1, x2)
;      plot, optspec & oplot, total(dataim[x1:x2, *], 1)
;
;
; MODIFICATION HISTORY:
;      Written by:      Dara Zeehandelaar, 4 March 2003
;      25 March 2003:   omitted loops in favor of matrix math
;      10 Dec 2003:     JD: added verbosity levels
;
;      $Log: extrspec.pro,v $
;      Revision 3.1  2004/07/16 17:45:08  jfd28
;      Fixed so exres and exmask return the correct image
;
;      Revision 3.0  2004/07/08 20:07:51  jfd28
;      Moving up to 3.0 (reflects this version anyway)
;
;      Revision 2.0  2004/07/08 20:05:01  jfd28
;      Tagging all files to 2.0
;
;      Revision 1.4  2004/07/08 18:45:22  jfd28
;      added comments
;      uses ethresh instead of mthresh
;      creates exmask instead of crmask
;      creates exres instead of crim
;      plottype takes an array
;
;-

function EXTRSPEC, dataim, profim, varim, v0, Q, x1, x2, $
                   bgim = bgim, bpct = bpct, errvect = errvect, $
                   ethresh = ethresh, exmask = exmask, exres = exres, $
                   gotovect = gotovect, inmask = inmask, opvar = $
                   opvar, plottype = plottype, skyvar = skyvar, $
                   verbose = verbose

; SET DEFAULTS AND CHECK INPUTS

dims = size(dataim)
nx = dims[1]
ny = dims[2]

if not keyword_defined(inmask  ) then inmask   = bytarr(nx,ny) + 1
if not keyword_defined(ethresh ) then ethresh  = 5
if not keyword_defined(skyvar  ) then skyvar   = dblarr(nx,ny)
if not keyword_defined(bgim    ) then bgim     = dblarr(nx,ny)
if not keyword_defined(gotovect) then gotovect = -1
if not keyword_defined(verbose ) then verbose  = 0
if not keyword_defined(plottype) then plottype = [0, 0, 0, 0]

if (verbose gt 1) then print, "Into optimal extraction"

str = ""
str += x1 lt 0 || x1 gt x2                                ? "x1 "     : ""
str += x2 gt nx - 1 || x2 lt x1                           ? "x2 "     : ""
str += (size(varim ))[1] ne nx || (size(varim ))[2] ne ny ? "varim "  : ""
str += (size(profim))[1] ne nx || (size(profim))[2] ne ny ? "profim " : ""
str += (size(inmask))[1] ne nx || (size(inmask))[2] ne ny ? "inmask " : ""
str += (size(bgim  ))[1] ne nx || (size(bgim  ))[2] ne ny ? "bgim "   : ""
str += (size(skyvar))[1] ne nx || (size(skyvar))[2] ne ny ? "skyvar " : ""
str += 4 ne n_elements(plottype)                          ? "plottype" : ""
if str ne "" then message, "Poorly formed inputs: " + str

; INITIALIZE
exres = dblarr(nx, ny)           ; an image with the residuals for each pixel
exmask = bytarr(nx, ny)+1       ; bad pixel mask created by this routine
crv = dblarr(x2-x1+1)           ; the residuals for a row
optspec = dblarr(ny)            ; the extracted spectrum
opvar = dblarr(ny)              ; the variance of the optimal extraction
errvect = bytarr(ny) + 1        ; the location of very bad rows
func = "extractfunc"


; LOOP THRROUGH ROWS
; Cut up the images into rows and pass each to procvect.  Tells
; procvect to use extractfunc.pro which returns the optimal
; extraction.  It also return in parm the variance of the extraction.
; Procvect will handle the sigma rejection.
for i = 0, ny-1 do begin
  varv    = varim [x1:x2, i]    ; cut up array into rows: variance row
  maskv   = inmask[x1:x2, i]    ; input mask row
  datav   = dataim[x1:x2, i]    ; data row
  multv   = profim[x1:x2, i]    ; profile row
  bgv     = bgim  [x1:x2, i]    ; background image row
  skyvarv = skyvar[x1:x2, i]    ; sky variance row
  if (i eq gotovect) then verbose = 5
  if plottype[1] or verbose eq 5 then begin
    device, window_state = ws   ; get the present window state
    if not ws[12] then window, 12 else wset, 12 ; open window 12
    !p.multi = [0, 2, 3, 1, 1]  ; fit all six inside
    !p.charsize = 2.0
    !x.margin = [8, 2]
    !y.margin = [2, 2]
    plot, datav,   title = 'Data Vector for Extraction', /ystyle
    plot, multv,   title = 'Profile', /ystyle
    plot, maskv,   title = 'Input Mask', yrange = [-0.1, 1.1]
    plot, varv,    title = 'Variance', /ystyle
    plot, skyvarv, title = 'Sky variance', /ystyle
    plot, bgv,     title = 'Background', /ystyle
    !p.multi = 0
    !p.charsize = 1.0
    !x.margin = [10, 3]
    !y.margin = [4, 2]
    wait, 0.01                  ; pause to allow user to view
  endif
  if verbose eq 5 then stop     ; stop for debugging purposes
  optspec[i] = procvect(datav, varv = varv, multv = multv, verbose = verbose, $
                        thresh = ethresh, maskv = maskv, crv = crv, /extract, $
                        bgv = bgv, Q = Q, v0 = v0, skyvarv = skyvarv, $
                        errflag = errflag, opvar = opvarv, vectnum = i, $
                        func = func, parm = parm, bpct = bpct, $
                        plottype = plottype)
  if (errflag) then errvect[i] = 0
  opvar[i] = parm               ; the optimal spectrum's variance
  exmask[x1:x2, i] = maskv      ; the mask created during extraction
  varim[x1:x2, i]  = varv       ; the final variance from processed values
  exres[x1:x2, i]   = crv * maskv ; the residual image
endfor

return, optspec

end
