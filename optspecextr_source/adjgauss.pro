;+
; NAME:
;     ADJGAUSS
;
;
; PURPOSE:
;     Stadardize the pixel locations when the profile image is a
;     Gaussian with quickly changing center and/or height.
;
;
; CATEGORY:
;     Optimal Spectrum Extraction Package
;         - Image Adjustment
;
;
; CALLING SEQUENCE:
;      Result = ADJGAUS (InArray, OptInfo, AdjParms=AdjParms,
;                        AdjOptions=AdjOptions, Revert=Revert)
;
;
; INPUTS:
;     INARRAY: An [nx, ny, nimages] array of images to adjust.
;     OPTINFO: A structure of images and info from the optimal extraction:
;       .DATAIM:   the image to use when finding te profile
;       .VARIM:    the variance image for weighting of the poly_fit
;       .SKYVAR    the cooresponding sky variance image (0 array)
;       .BGIM:     the background image to subtracted off the data (0 array)
;       .INMASK:   the mask used for all functions (1 array)
;       .Q:        the gain of the image (1)
;       .V0:       the readnoise ^2 of the data (0)
;       .X1:       the start of the object's profile (0)
;       .X2:       the end of the object's profile (nx-1)
;       .SPEC:     the standard extracted spectrum (total(dataim[x1:x2,*],1)
;       .BPCT:     The percentage of allowable bad pixels before halting (.5)
;       .VERBOSE:  level of output to screen (0)
;       .PLOTTYPE: the type of plot to output while running ([0, 0, 0, 0])
;
;
; OPTIONAL INPUTS:
;      Profile Fitting - Straightening:
;      STRAIGHTEN: Set to straighten the trace
;      EXPAND:     level of expansion in straightening routine (11)
;      CENTROID:   Set to straighten with centroid instead of Gaussian
;      TRACEDEG:   degree of fitting for estimating shift in trace (2)
;      NOTRACEFIT: Set to not fit the estimated trace
;      GAUSSTH:    error allowed during row's center calculation (0.05)
;      SHIFTTH:    error allowed between the center and the fitted trace (0.5)
;      DEGCONTR:   degree of fit to use when contracting the array (2)
;
;
; KEYWORD PARAMETERS:
;     ADJPARMS:   A structure for paramaters of the adjustment with 
;       .ORIGX:      the x values for each pixel in the original images
;       .ADJX:       the x values for each pixel in the adjusted images
;       .TRACEEST:   the estimated center for each row
;       .WIDTHEST:   the estimated width for each row;
;     ADJOPTIONS: A structure with options for the gaussian adjustment
;       .LEVEL:      The amount to expand the array
;       .CENTER:     Set to 1 to adjust each column so the centers line up (1)
;       .WIDTH:      Adjust to reflect each pixel's relative distance 
;                    from center according to that row's Gaussian width (0)
;       .CENTROID:   Set to use centroid fitting instead
;       .GAUSSTH:    The absolute value of allowable error when fitting a
;                    gaussian to the data (3 percent)
;       .SHIFTTH:    The absolute allowable error between a row's center
;                    and the fitted trace (0.5 pixels)
;       .DEGCONTR:   degree of fit to use when contracting the array (2)
;       .CENTERDEG:  Degree of fit for estimating shift (2)
;       .CENTERFIT:  Set to not fit the estimated centers and not accept the
;                    fitting to be correct
;     REVERT = Set to use the data in AdjParms to convert the images
;              in the input array back into the origional shape
;
;
; OUTPUTS:
;     An array ((x2 -x1 + 1)*level by ny by nimages) where [*, *, i]
;     standardized image from the corresponding image in the input array
;
;
; PROCEDURE:
;     Unpack the structures first.  Then find the width and centers
;     using findshift.  Finally adjust each image in InArray so that
;     the profiles line up.
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;     Written by John Dermody 24 Jun 2004
;
;     $Log: adjgauss.pro,v $
;     Revision 3.1  2004/07/13 15:38:34  jfd28
;     fixed header to reflect centerest not traceest as the keyword
;
;     Revision 3.0  2004/07/08 20:07:51  jfd28
;     Moving up to 3.0 (reflects this version anyway)
;
;     Revision 2.0  2004/07/08 20:05:01  jfd28
;     Tagging all files to 2.0
;
;     Revision 1.2  2004/07/08 18:41:44  jfd28
;     check to see if tags in optinfo are defined before using
;     made a plot for in coordinate array and out coordinate array
;     plottype is an array
;
;     Revision 1.1  2004/07/02 18:26:37  jfd28
;     initial implementation
;
;     pulled out all of fitprof's straightening statements and gave it its
;     own file so that optimal spectrum extraction can continue with no
;     knowledge of what is occuring here.  To do that I needed to use a
;     structure to pass the parmaters of the fit around inside fitprof
;     and use a structure to pass options down to this function.
;
;     uses a data cube to pass down the frames fitprof would like to modify
;     and a structure for all the info optspectextr has up to that point in
;     running the program.
;
;-

function ADJGAUSS, InArray, OptInfo, AdjParms = AdjParms, $
  AdjOptions = AdjOptions, Revert = Revert


; INITIALIZE AND CHECK INPUTS

; Sizes
nx     = (size(inarray))[1]
ny     = (size(inarray))[2]
nimage = (size(inarray))[3]

; Unpack Optimal Spectrum Extraction Info
otagn    = Tag_Names(OptInfo)
dataim   = OptInfo.Dataim
varim    = OptInfo.Varim 
verbose  = in("VERBOSE",  otagn) ? OptInfo.Verbose  : [0, 0, 0, 0]
plottype = in("PLOTTYPE", otagn) ? OptInfo.Plottype : 0
if in("INMASK",   otagn) then inmask = OptInfo.Inmask
if in("SKYVAR",   otagn) then skyvar = OptInfo.Skyvar
if in("BGIM",     otagn) then bgim   = OptInfo.Bgim
if in("X1",       otagn) then X1     = OptInfo.X1
if in("X2",       otagn) then X2     = OptInfo.X2
if in("SPEC",     otagn) then Spec   = OptInfo.Spec

; Unpack Adjustment Procedure Options
atagn = keyword_defined(AdjOptions) ? Tag_Names(AdjOptions) : [""]
if in("LEVEL",     atagn) then level     = AdjOptions.Level  else level  = 7
if in("CENTER",    atagn) then center    = AdjOptions.Center else center = 1
if in("WIDTH",     atagn) then width     = AdjOptions.Width  else width  = 0
if in("CENTROID",  atagn) then centroid  = AdjOptions.Centroid
if in("GAUSSTH",   atagn) then gausth    = AdjOptions.Gaussth
if in("SHIFTTH",   atagn) then shiftth   = AdjOptions.Shiftth
if in("CENTERDEG", atagn) then centerdeg = AdjOptions.Centerdeg
if in("CENTERFIT", atagn) then centerfit = AdjOptions.Centerfit
if in("WIDTHTH",   atagn) then widthth   = AdjOptions.Widthth
if in("DEGCONTR",  atagn) then degcontr  = AdjOptions.Degcontr

; Check inputs
if level lt 1 then message, "Level must be greater or equal to 1"


; REVERT TO ORIGINAL IMAGE SPECIFICATIONS
; Grab the old geometry's coordinate array and the new coordinate array and
; use sampleshift to estimate the pixel value at the old geometry using
; closest pixels in the new geometry
if keyword_set(revert) then begin
  origx = AdjParms.origx
  adjx  = AdjParms.adjx
  nxo = (size(origx))[1]
  outarr = dblarr(nxo, ny, nimage)
  if level le 2 then begin      
    fitpoly = 0                 ; If there is less than 3 in pixels per output
    fitspline = 1               ; pixel then the polynomials formed for
  endif else begin              ; estimation will be werid.
    fitspline = 0               ; Use splines insted
    fitpoly = 1 
  endelse
  for i = 0, nimage - 1 do begin
    outarr[*, *, i] = sampleshift(inarray[*, *, i], adjx, origx, $
                                  fitpoly = fitpoly, degcontr = degcontr, $
                                  fitspline = fitspline)
  endfor
  return, outarr
endif


; FIND THE PROFILE'S CENTER AND WIDTH
; Findshift will return the center and width for each row. Center value is the
; x value in the data array.
saveverbose = verbose
traceest = findshift(dataim, varim = varim, inmask = inmask, spec = spec, $
                     verbose = verbose, plottype = plottype, $
                     x1 = x1, x2 = x2, $
                     centroid = centroid, centerdeg = centerdeg, $
                     gaussth = gaussth, shiftth = shiftth, $
                     centerfit = centerfit, height = height, $
                     widthv = widthv, tracemask = tracemask)
verbose = saveverbose
if (verbose eq 5) then stop

; Check and see if the user wants to adjust row by row
if not center then traceest = dblarr(ny) + nx/2
if not width  then widthv   = dblarr(ny) + mean(widthv)

; ADJUST THE GEOMETRY OF THE FRAMES

; Create the coordinate grids
nxo = nx * level
origx = (dindgen(nx) # (dblarr(ny) + 1) - traceest ## (dblarr(nx) + 1)) / $
        (widthv ## (dblarr(nx) + 1))
maxx = max(origx)
minx = min(origx)
adjx  = dindgen(nxo, ny) mod (nxo) / nxo * (maxx-minx) + minx

if plottype[1] then begin
  device, window_state = ws
  if not ws[12] then window, 12 else wset, 12  
  !p.multi = [0, 1, 2, 1, 1]
  origy = findgen(ny) ## (fltarr(nx) + 1)
  adjy = findgen(ny) ## (fltarr(nxo) + 1)
  plot, origy, origx, psym = 3, title = 'Old X coordinates', /ystyle
  plot, adjy, adjx, psym = 3, title = 'New X coordinates', /ystyle
  !p.multi = 0
  wait, 1.5
endif

; Adjust the images in InArray
outarray = dblarr(nxo, ny, nimage)
for i = 0, nimage-1 do begin
  outarray[*, *, i] = sampleshift(inarray[*, *, i], origx, adjx, /fitspline)
endfor

; Pack up parameters of adjustment
AdjParms = {origx:origx, adjx:adjx, traceest:traceest, widthest:widthv, $
           mask:tracemask}


return, outarray
end
