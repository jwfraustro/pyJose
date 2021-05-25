;+ NAME:
;      FITPROF
;
;
; PURPOSE:
;      Creates an image of normalized spatial profiles, enforces positivitiy,
;      normalization.  If needed, can also shift and expand the data image
;      prior to fitting.  By default fits a polynomial in the spacial
;      direction
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector Setup Function
;
; CALLING SEQUENCE:
;      Result = FITPROF (Dataim, Spec, x1, x2, adjfunc = adjfunc, $
;                        adjparms = adjparms, adjoptions = adjoptions, $
;                        bgim = bgim, boxcarhw = boxcarhw, bpct = bpct, $
;                        difpmask = difpmask, errvect=errvect, fitboxcar = $
;                        fitboxcar, fitgauss = fitgauss, gotovect = gotovect, $
;                        inmask = inmask, noproffit = noproffit, plottype = $
;                        plottype, profdeg = profdeg, profmask = profmask, $
;                        profres = profres, pthresh = pthresh, q = q, skyvar $
;                        = skyvar, varim = varim, verbose = verbose, v0 = v0)
;
;
; INPUTS:
;      DATAIM:  sky-subtracted processed image, horizontal is pixel,
;               vertical is wavelength
;      SPEC:    array holding extracted spectra, same dimension as
;               vertical of datatim
;      X1,X2:   boundaries in x which contain spectrum (spectrum is
;               inclusive of X1,X2)
;
;
; OPTIONAL KEYWORDS:
;      General
;      VARIM:      the variance image for weighting of the poly_fit (1 array)
;      SKYVAR:     the cooresponding sky variance image (0 array)
;      BGIM:       the background image to subtracted off the data (0 array)
;      INMASK:     the mask used for all functions (1 array)
;      Q:          the gain of the image (1)
;      v0:         the readnoise ^2 of the data (0)
;      BPCT:       The percentage of allowable bad pixels before halting (.50)
;
;      Geometry Modification
;      ADJFUNC:    function to use to conduct any geometry manipulation
;      ADJPARMS:   the parameters of the geometry for this data set
;      ADJOPTIONS: any options that can be set for the adjfunc
;      
;      Smoothing the Profile Options
;      NOPROFFIT:  if set, will not smooth profile image
;      FITGAUSS:   set to fit profile in spacial direction with a
;                  gaussian
;      FITBOXCAR:  set to fit the profile with a boxcar filter
;      BOXCARHW:   the width of the median in boxcar average (3)
;      PTHRESH:    the threshold for sigma rejection (3)
;      PROFDEG:    degree of smoothing for profile using a polyfit, (2)
;
;      Output from Smoooth Profile Fitting
;      PROFMASK:   the mask modified by the fit vector routine
;      PROFRES:    the final residuals for fitting
;      ERRVECT:    Will contain the vectors that are recovered
;      DIFPMASK:   The pixels rejected during the profile fitting
;                  routine.
;
;      Debugging Help:
;      VERBOSE:    level of output to screen (0)
;      PLOTTYPE:   The type of plot to output while running ([0,0,0,0])
;      GOTOVECT:   stop when this vector is reached (none)
;
;
; OUTPUTS:
;      returns an image, same dimensions as dataim,containing the
;      spatial profile, (smoothed?), positivity enforced, normalized
;
;
; EXAMPLE:
;      ; create an object frame using the
;      ; synthetic spectrum package
;      obj = synthspec(sky, objheader = objheader, /notrace, tilt =  1)
;      rn =  sxpar(objheader, 'RDNOISE')
;      gain = sxpar(objheader, 'EPADU')
;      dataim = obj - sky
;      ; estimate X1 and X2 by finding the
;      ; base of the tallest peak
;      t = total(dataim,  2)
;      mx = max(t,  cnt)
;      x1 = cnt - 5 & while (t[x1-1] lt t[x1]) do x1--
;      x2 = cnt + 5 & while (t[x2+1] lt t[x2]) do x2++
;      varim = obj / gain + rn^2 + sky / gain + rn^2
;      skyvar = sky*2 / gain + rn^2
;      spec = total(dataim[x1:x2,  *],  1)
;      profim = fitprof(dataim, spec, x1, x2, varim = varim,  Q = gain,  $
;                       v0 = rn^2, profdeg = 2,  /straighten)
;
;
; PROCEDURE:
;      Horne Step 5
;      Divide the sky-subtracted data image by the spectrum to get the initial
;      profile image. If noproffit is set return this.  Else, if straighten is
;      set, find the center of the profile for each wavelength (using
;      findshift), and expand and straighten the array so these centers line
;      up (using sampleshift).  Now, fit a polynomial, Gaussian, or running
;      average to the image, and use sigma rejection to reject bad pixels
;      (using the procvect procedure).  Returned the fitted data, resizing if
;      necessary (using sampleshift again).
;
;
; FUNCTIONS CALLED:
;      procvect, findshift, centriod, sampleshift
;
;
; MODIFICATION HISTORY:
;       Written by:      Dara Zeehandelaar, 2 March 2003
;       25 March 2003:   omitted loops in favor of matrix math
;                        removed dimension check (will be done in driver)
;       1 April 2003:    refined matrix math
;       15 April 2003:   added degprof keyword, finished smoothing
;       19 April 2003:   debug normalizing
;       16 June 2003:    John Dermody: inserted weighting and a sigma
;                        clip algorthm inton the smoothing function
;       8 Jul 2003:      JD: removed fit and rejection to its own
;                        function
;       1 Sept 2003:     JD: added spacial fit and cleaning
;       19 Oct 2003:     JD: fixed header
;       03 Dec 2003:     JD: Added recovery
;       18 Dec 2003:     JD: cleaning
;
;       $Log: fitprof.pro,v $
;       Revision 3.0  2004/07/08 20:07:51  jfd28
;       Moving up to 3.0 (reflects this version anyway)
;
;       Revision 2.0  2004/07/08 20:05:01  jfd28
;       Tagging all files to 2.0
;
;       Revision 1.11  2004/07/08 18:40:51  jfd28
;       no recovery
;       fixed ranges on shade_surf
;       plottype takes an array
;       profdeg instead of degprof
;
;       Revision 1.10  2004/07/02 18:30:21  jfd28
;       removed all shifting options and statements to adjgauss.pro.  Sets up
;       a structure to pass the necessary info across and assumes the data
;       coming in from adjfunc is good.  Also removed recovery information
;
;       Revision 1.9  2004/06/22 21:26:38  jfd28
;       cleaned up the header added catagory and example.  organized keywords
;
;       made nosmooth return with correctly initialized outputs
;
;       pulled out the shift finding procedure and placed into its
;       own routine (findshift.pro)
;
;       pulled out the expansion and shift calculations so now just send
;       in the data an array of x values, and an array of locations
;       where values are wanted
;
;       Revision 1.8  2004/06/02 21:27:24  jfd28
;       pushed degcontr up andd removed the extraction of expanded data
;       removed remaking the spectrum from expanded data
;       expanded origional data is saved, so varout can be placed inside
;
;       Revision 1.7  2004/05/29 20:59:46  jfd28
;       Shift the whole array, then pull out the profile
;       gotovect sets verbose to 5
;       added feature to extract from the expanded array
;
;       Revision 1.6  2004/05/27 21:00:32  jfd28
;       sending gaussth and shiftth back to top
;
;       Revision 1.5  2004/05/27 20:56:46  jfd28
;       if gotovect is hit then set verbose to 5
;
;+

function FITPROF, Dataim, Spec, x1, x2, adjfunc = adjfunc, $
                  adjparms = adjparms, adjoptions = adjoptions, $
                  bgim = bgim, boxcarhw = boxcarhw, bpct = bpct, $
                  difpmask = difpmask, errvect = errvect, fitboxcar = $
                  fitboxcar, fitgauss = fitgauss, gotovect = gotovect, $
                  inmask = inmask, noproffit = noproffit, plottype = $
                  plottype, profdeg = profdeg, profmask = profmask, $
                  profres = profres, pthresh = pthresh, q = q, skyvar $
                  = skyvar, varim = varim, verbose = verbose, v0 = v0

; SET DEFUALTS AND CHECK INPUTS
dims = size(dataim)
nx = dims[1]
ny = dims[2]

if not keyword_defined(varim   ) then varim    = dblarr(nx, ny) + 1
if not keyword_defined(inmask  ) then inmask   = bytarr(nx, ny) + 1
if not keyword_defined(pthresh ) then pthresh  = 3
if not keyword_defined(bgim    ) then bgim     = dblarr(nx, ny)
if not keyword_defined(skyvar  ) then skyvar   = dblarr(nx, ny)
if not keyword_defined(profdeg ) then profdeg  = 3
if not keyword_defined(boxcarhw) then boxcarhw = 3
if not keyword_defined(verbose ) then verbose  = 0
if not keyword_defined(plottype) then plottype = [0, 0, 0, 0]
if not keyword_defined(gotovect) then gotovect = -1

str = ""
str += n_elements(spec) ne ny                             ? "spec "   : ""
str += x1 lt 0 || x1 gt x2                                ? "x1 "     : ""
str += x2 gt nx - 1 || x2 lt x1                           ? "x2 "     : ""
str += (size(varim ))[1] ne nx || (size(varim ))[2] ne ny ? "varim "  : ""
str += (size(inmask))[1] ne nx || (size(inmask))[2] ne ny ? "inmask " : ""
str += (size(bgim  ))[1] ne nx || (size(bgim  ))[2] ne ny ? "bgim "   : ""
str += (size(skyvar))[1] ne nx || (size(skyvar))[2] ne ny ? "skyvar " : ""
str += 4 ne n_elements(plottype)                          ? "plottype" : ""
if str ne "" then message, "Poorly formed inputs: " + str

if (verbose gt 1) then print, "Into profile fitting"


; NO PROFILE SMOOTHING
; if the user doens't want to fit the profile, he can just divide the data
; image by the spectrum at that row and hope for the best.  This is not
; recomended because it is sensitive to the standard extraction which is
; sensitive to bad pixels.  The image is normalized and made everywhere
; greater than 0.
if keyword_set(noproffit) then begin
  specim = spec ## (dblarr(nx) + 1)
  profim = dataim / specim
  profim = profim * inmask > 0
  t = total(profim[x1:x2, *], 1)
  t = t ## (dblarr(nx) + 1)
  profim = profim / t
  profmask = inmask             ; initialize unused outputs
  difpmask = profmask - profmask
  errvect = bytarr(ny) + 1
  return, profim
endif


; ADJUST FRAMES
; If the user wishes to modify the geometry of the images before profile
; fitting (say straightening the trace) pass all the info known so far, and
; the images needed profile fitting to adjfunc.
if keyword_set(AdjFunc) then begin
  ; Pass information as a structure for simplicity
  OptInfo = {Dataim:dataim, Varim:varim, Bgim:Bgim, Skyvar:skyvar, Spec:Spec, $
             Inmask:Inmask, X1:X1, X2:X2, Verbose:Verbose, Plottype:Plottype}
  if keyword_defined(bpct) then OptInfo = create_struct('BPCT', bpct, OptInfo)
  if keyword_defined(Q   ) then OptInfo = create_struct('Q'   , Q   , OptInfo)
  if keyword_defined(V0  ) then OptInfo = create_struct('V0'  , v0  , OptInfo)
  inarray  = dblarr(x2-x1+1, ny, 5)
  inarray[*, *, 0] = dataim[x1:x2, *]
  inarray[*, *, 1] = varim [x1:x2, *]
  inarray[*, *, 2] = inmask[x1:x2, *]
  inarray[*, *, 3] = skyvar[x1:x2, *]
  inarray[*, *, 4] = bgim  [x1:x2, *]
  saveverbose  = verbose
  saveplottype = plottype
  outarray = call_function(AdjFunc, InArray, OptInfo, AdjParms = AdjParms, $
                           AdjOptions = AdjOptions)
  verbose  = saveverbose
  plottype = saveplottype
  pdataim = outarray[*, *, 0]
  pvarim  = outarray[*, *, 1]
  pmask   = outarray[*, *, 2]
  pskyvar = outarray[*, *, 3]
  pbgim   = outarray[*, *, 4]
  ; estimate the influence of a bad pixel on its neighbors
  pmask = (pmask gt 0.99) and (pmask lt 1.01)
endif else begin
  pdataim = dataim[x1:x2, *]
  pvarim  = varim [x1:x2, *]
  pmask   = inmask[x1:x2, *]
  pskyvar = skyvar[x1:x2, *]
  pbgim   = bgim [x1:x2, *]
endelse


; INITIALIZE
; Because the images may be a different size now, new initilization is
; needed.  If a Gaussian profile is assumed, transverse by rows.  If
; polynomial fitting or running average fitting is used loop by columns.
pnx = (size(pdataim))[1]
pny = (size(pdataim))[2]

if keyword_set(fitgauss) then begin
  f1 = 0                        ; fit over spatial dimension
  f2 = pny - 1
  nvect = pnx
  ndata = pny
  func = "gaussfunc"
  parm = 1
endif else begin
  f1 = 0                        ; fit over spectral dimension
  f2 = pnx - 1
  nvect = pny
  ndata = pnx
  if keyword_set(fitboxcar) then begin
    func = "boxcarfunc"
    parm = boxcarhw
  endif else begin
    func = "polyfunc"
    parm = profdeg
  endelse
endelse

profmask = pmask
profres  = dblarr(pnx, pny)
pprofim = dblarr(pnx, pny)
pspecim = spec ## (dblarr(pnx) + 1)
datav   = dblarr(nvect)         ; initialize so dimensions are correct
maskv   = dblarr(nvect)
varv    = dblarr(nvect)
multv   = dblarr(nvect)
skyvarv = dblarr(nvect)
bgv     = dblarr(nvect)
errvect = dblarr(ndata) + 1
index = dindgen(pnx, pny)
pprofim = pdataim / pspecim
yvals  = indgen(pny)
xvals  = indgen(pnx)


; LOOP THROUGH ROWS OR COLUMNS
; Cut up the data according to the type of fitting and pass the info into
; procvect.  Provect will take care of bad pixel rejection and pass the
; smoothed profile back.
for i = f1, f2 do begin
  if keyword_set(fitgauss) then begin
    is = index[*, i]
  endif else begin
    is = index[i, *]
  endelse
  datav[*]   = pdataim[is]
  maskv[*]   = pmask[is]
  varv[*]    = pvarim[is]
  crv        = profres[is]
  multv[*]   = pspecim[is]
  skyvarv[*] = pskyvar[is]
  bgv[*]     = pbgim[is]
  if (i eq gotovect) then verbose = 5
  if plottype[1] or (verbose eq 5) then begin
    device, window_state = ws   ; Debugging Plot #1
    if not ws[12] then window, 12 else wset, 12
    !p.multi = [0, 2, 3, 1, 1]
    !p.charsize = 2.0
    !x.margin = [8, 2]
    !y.margin = [2, 2]
    plot, datav,   title = 'Data Vector for Profile Fitting', /ystyle
    plot, multv,   title = 'Std Spectrum', /ystyle
    plot, maskv,   title = 'Input Mask', yrange = [-0.1, 1.1]
    plot, varv,    title = 'Variance', /ystyle
    plot, skyvarv, title = 'Sky variance', /ystyle
    plot, bgv,     title = 'Background', /ystyle
    !p.multi = 0
    !p.charsize = 1.0
    !x.margin = [10, 3]
    !y.margin = [4, 2]
    wait, 0.01
  endif
  if (verbose eq 5) then stop
  errflag = 0
  pprofim[is] = procvect(datav, varv = varv, thresh = pthresh, multv = multv, $
                         maskv = maskv, crv = crv, bgv = bgv, Q = Q, v0 = v0, $
                         skyvarv = skyvarv, verbose = verbose, $
                         errflag = errflag, vectnum = i, func = func, $
                         parm = parm, bpct = bpct, plottype = plottype)
  if (errflag) then errvect[i] = 0
  if plottype[3] then begin     ; Debugging Plot #3
    device, window_state = ws
    if not ws[14] then window, 14 else wset, 14
    shx1 = ((i - 10) > f1)
    shx2 = ((i + 10) < f2)
    if keyword_set(fitgauss) then begin
      sis = index[*, shx1:shx2]
      shade_surf, pprofim[sis], xvals, yvals[shx1:shx2], charsize = 3, $
                  title = "Profile after # " + strtrim(i, 1), $
                  xtitle = "X pixel location", $
                  ytitle = "Y (Prefit and Postfit)"
    endif else begin
      sis = index[shx1:shx2, *]
      shade_surf, pprofim[sis], xvals[shx1:shx2], yvals, charsize = 3, $
                  title = "Profile after # " + strtrim(i, 1), $
                  xtitle = "X (prefit and postfit)", $
                  ytitle = "Y pixel location"
    endelse
    wait, 0.01
  endif
  profmask[is] = maskv
  pvarim[is] =   varv
  profres[is] =  crv * maskv
endfor


; REVERT TO ORIGINAL ARRAY
profim = dblarr(nx, ny)
if keyword_set(AdjFunc) then begin
  inarray = dblarr(pnx, pny, 2)
  inarray[*, *, 0] = pprofim
  inarray[*, *, 1] = pvarim
  outarr = call_function(AdjFunc, InArray, OptInfo, AdjParms = AdjParms, $
                         AdjOptions = AdjOptions, /revert)
  profim = dblarr(nx, ny)
  profim[x1:x2, *] = outarr[*, *, 0]
  varim [x1:x2, *] = outarr[*, *, 1]
endif else begin
  varim [x1:x2, *] = pvarim
  profim[x1:x2, *] = pprofim
endelse


; NORMALIZE
profim = profim > 0
t = total(profim[x1:x2, *], 1)
t = t ## (dblarr(nx) + 1)
profim = profim / t
difpmask = (float(pmask) - profmask)

return, profim

end
