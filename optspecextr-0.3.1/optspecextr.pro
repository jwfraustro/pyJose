;+
; NAME:
;      OPTSPECEXTR
;
; PURPOSE:
;      data pipeline for image processing and optimal extraction
;
; CALLING SEQUENCE:
;      Result = OPTSPECEXTR (data, var, rn, q, x1, x2, $
;        adjfunc = adjfunc, adjoptions = adjoptions, adjparms = adjparms,   $
;        berrvect = berrvect, bgdeg = bgdeg, bgim = bgim, bgmask = bgmask, $
;        bgotovect = bgotovect, bgres = bgres, boxcarhw = boxcarhw, bpct = $
;        bpct, bthresh = bthresh, debughead = debughead, difpmask = $
;        difpmask, eerrvect =  eerrvect, egotovect = egotovect, ethresh = $
;        ethresh, exmask = exmask, exres = exres, fitboxcar = fitboxcar, $
;        fitgauss = fitgauss, inmask = inmask, integrate = integrate, $
;        nobgfit = nobgfit, noproffit = noproffit, opvar = opvar, perrvect $
;        = perrvect, pgotovect = pgotovect, plottype = plottype, profdeg = $
;        profdeg, profim = profim, profmask = profmask, profres = profres, $
;        pthresh = pthresh, skyvar = skyvar, stdspec = stdspec, stdvar = $
;        stdvar, varout = varout, verbose = verbose)
;
; INPUTS:
;      DATA:  processed image (raw-sky)/flat for IR, (raw-bias)/flat
;             for optical (by Horne step 1)
;      VAR:   variance image of data in DN^2 (i.e., Data/Q+Rn^2, Horne step 2)
;      RN:    RN is read noise of data frame in DN
;      Q:     effective number of photons per DN for data frame
;      X1,X2: boundaries in x of spectrum, spectrum is inclusive of X1,X2
;
; OPTIONAL KEYWORDS: (Input)
;      General:
;      SKYVAR:     Variance image of sky (by Horne step 2) (default 0s)
;      INMASK:     User input mask with good pixels=1, bad=0 (default 1s)
;      BPCT:       Set the percentage of bad pixels before kicking out (0.5)
;      INTEGRATE:  Set to integrate overbad pixels when calculating
;                  spec for use in the profile fitting
;
;      Background Fitting:
;      NOBGFIT:    set to not fit the background
;      BTHRESH:    sigma threshold for cosmic ray rejection for bg, (3)
;      BGDEG:      degree of polynomial fit for interpolation in fitbg (1)
;
;      Profile Fitting - Geometry Adjustment
;      ADJFUNC:    the name of the function to call with images to modify
;      ADJOPTIONS: any options needed for the function
;      ADJPARMS:   any parmaters of the fit needing to be saved
;
;      Profile Fitting - Smoothing:
;      NOPROFFIT:  will not smooth spatial profile if set
;      PTHRESH:    sigma threshold for profile image, (3)
;      PROFDEG:    degree to smooth spatial profile, (3)
;      FITGAUSS:   Set to smooth the profile image using a Gaussian fit
;                  in the spatial direction, instead a polynomial fit in
;                  wavelength direction
;      FITBOXCAR:  Set to smooth the profile using a median filter in
;                  spectral dimension
;      BOXCARHW:   The half-width of the boxcar fit (5)
;
;      Optimal Extraction:
;      ETHRESH:    sigma threshold for main mask, (5)
;
;      Debuging Assitance:
;      VERBOSE:    Set to level of output: 0 - fatal, 1 - collect
;                  stats, 2 - where you are in the program, 3 - output
;                  every warning, 4 - plot how the data are being used,
;                  5 - stop after each iteration.
;      PLOTTYPE:   Set to the type of plot outputed to the screen:
;                  0 - no plots
;                  1 - show a summary of all processed frames
;                  2 - show the data used as input to procvect
;                  3 - show how procvect is fitting the data
;                  4 - use shade_surf to show the fitted and raw data
;                  5 - show all plots
;                  Also can take a 4 element array where
;                  PLOTTYPE[0] - show a summary of all processed frames
;                  PLOTTYPE[1] - show the data used as input to procvect
;                  PLOTTYPE[2] - show how procvect is fitting the data
;                  PLOTTYPE[3] - use shade_surf to show the fitted and raw
;      BGOTOVECT:  A tool to halt execution at that row in fitbg.
;      PGOTOVECT:  A tool to halt execution at that vector in fitprof.
;      EGOTOVECT:  A tool to halt execution at that row in specextr.
;      DEBUGHEAD:  Text to place at the top of the proccessed array window
;
; OPTIONAL KEYWORDS: (Output)
;      Background Fitting:
;      BGIM:     interpolated sky background image
;      BGMASK:   background's mask
;      BGRES:    the sigma difference between the actual and expected
;      BERRVECT: mask of rows where iteration stopped during fitbg
;
;      Standard Extraction:
;      STDSPEC:  standard spectrum
;      STDVAR:   returns variance of standard spectrum
;
;      Profile Fitting:
;      PROFIM:   calculated profile image
;      PROFMASK: the calculated profile cosmic ray mask
;      DIFPMASK: pixels rejected during the profile fitting routine.
;      PROFRES:  the profile residual image
;      PERRVECT: mask of rows/columns where iteration stopped during fitprof
;
;      Optimal Extraction:
;      OPVAR:    returns optimally extracted variance
;      VAROUT:   returns the final caculated variance
;      EXMASK:   calculated cosmic ray mask
;      EXRES:    cosmic ray image: dataim-spec*profim*profim/varim
;      EERRVECT: mask of rows where iteration stopped during extrspec
;
;
; OUTPUTS:
;      Returns the optimally extracted spectrum
;
; PROCEDURE:
;      implements Horne optimal extraction algorithm
;
; FUNCTIONS CALLED:
;      Horne step 3: fitbg
;      Horne step 4: stdexr
;      Horne step 5, 6: fitprof
;      Horne step 6, 7, 8: extrspec
;      summaryopt
;
; MODIFICATION HISTORY:
;
;       Written by:      Dara Zeehandelaar, 13 March 2003
;       23 April 2003    took out Horne steps 1 and 2, user does
;                        outside
;       30 April 2003    added optional user mask input
;       6 May 2003:      put skyproc within final iteration in order
;                        to use cosmic ray mask
;       7 May 2003:      added optional keyword edge
;       16 May 2003:     added keywords stdspec, crmask, bgim, crim,
;                        profim, itnum to return intermediate steps if desired
;       10 Aug 2003:     John Dermody: reorganized to have only 3
;                        functions called
;       20 Aug 2003:     JD: added main loop for extreme cases
;       1 Sept 2003:     JD: removed edge and cleaning up.
;       14 Oct 2003:     JD: turned everything into doubles
;       10 Dec 2003:     JD: set up verbosity levels
;       15 May 2004:     JD: set up cvs
;
;       $Log: optspecextr.pro,v $
;       Revision 3.0  2004/07/08 20:07:51  jfd28
;       Moving up to 3.0 (reflects this version anyway)
;
;       Revision 2.0  2004/07/08 20:05:01  jfd28
;       Tagging all files to 2.0
;
;       Revision 1.11  2004/07/08 19:42:58  jfd28
;       default values for plottype and verbose
;
;       Revision 1.10  2004/07/08 18:39:28  jfd28
;       final for 0.3.0
;
;       changed variable names to be consistent (bg, prof, ex).
;       Removed recovery stuff for simplicity.  Plottype uses an array
;
;       Revision 1.9  2004/07/02 18:35:09  jfd28
;       removed all straightening to an outside set of functions.
;
;       Revision 1.8  2004/07/02 16:21:40  jfd28
;       Simplified calling sequences
;       removed recovery and repitition of optimal extraction and
;       profile fitting
;
;       improved documentation
;
;       pulled all summary output to its own function
;
;       Revision 1.7 2004/06/02 21:25:22 jfd28
;       pushed degree of poly fit in contraction to the top so it can
;       be user defined removed bigdif and big spec stuff added
;       optional header to print to top of processed arrays window
;       scaled data and varout outputs to 96% of median instead of
;       full range
;
;       Revision 1.6  2004/05/29 21:36:20  jfd28
;       made gaussth and shiftth can be user-defined
;       pulled in code to:
;       if plottype eq 4 then plot optspec, stdspec, trace and
;       tv bgim, prof, var, data, and masks, and made it window 10.
;
;       Revision 1.5  2004/05/27 21:03:08  jfd28
;       Passing Gaussth and shiftth through
;
;       Revision 1.4  2004/05/27 20:55:40  jfd28
;       output cosmic ray mask also show pixels rejected during
;        background fitting
;
;-

function OPTSPECEXTR, data, var, rn, q, x1, x2, $
  adjfunc = adjfunc, adjoptions = adjoptions, adjparms = adjparms,   $
  berrvect = berrvect, bgdeg = bgdeg, bgim = bgim, bgmask = bgmask, $
  bgotovect = bgotovect, bgres = bgres, boxcarhw = boxcarhw, bpct = $
  bpct, bthresh = bthresh, debughead = debughead, difpmask = $
  difpmask, eerrvect =  eerrvect, egotovect = egotovect, ethresh = $
  ethresh, exmask = exmask, exres = exres, fitboxcar = fitboxcar, $
  fitgauss = fitgauss, inmask = inmask, integrate = integrate, $
  nobgfit = nobgfit, noproffit = noproffit, opvar = opvar, perrvect $
  = perrvect, pgotovect = pgotovect, plottype = plottype, profdeg = $
  profdeg, profim = profim, profmask = profmask, profres = profres, $
  pthresh = pthresh, skyvar = skyvar, stdspec = stdspec, stdvar = $
  stdvar, varout = varout, verbose = verbose

; REFORM AND CHECK INPUTS
if not keyword_defined(verbose ) then verbose  = 0
if not keyword_defined(plottype) then plottype = [0, 0, 0, 0]
if n_elements(plottype) eq 1 then begin
  case plottype of
    0:    plottype = [0, 0, 0, 0]
    1:    plottype = [1, 0, 0, 0]
    2:    plottype = [0, 1, 0, 0]
    3:    plottype = [0, 0, 1, 0]
    4:    plottype = [0, 0, 0, 1]
    else: plottype = [1, 1, 1, 1]
  endcase
endif

v0 = RN * RN                    ; read noise --> V0 (variances)
x1 = round(x1)                  ; x1 and x2 should be pixel values
x2 = round(x2)
varim = var                     ; don't modify user's variance

; Horne Step 3 - fit sky background
saveverbose = verbose
saveplottype = plottype
if (verbose eq 5) then stop
bgim = fitbg(data, x1, x2, bgdeg = bgdeg, Q = Q, v0 = v0, $
             inmask = inmask, varim = varim, bthresh = bthresh, $
             bgres = bgres, errvect = Berrvect, bgmask = bgmask, $
             skyvar = skyvar, verbose = verbose, nobgfit = nobgfit, $
             gotovect = Bgotovect, bpct = bpct, plottype = plottype)
verbose = saveverbose
plottype = saveplottype

; sky-subtract
dataim = data - bgim

; Horne Step 4 - extract standard spectrum, compute variance
stdspec = stdextr(dataim, varim, x1, x2, inmask = inmask, stdvar = stdvar, $
               adjspec = adjspec)
spec = keyword_set(integrate) ? adjspec : stdspec


; Horne Step 5 - construct spatial profile
; WITH STEP 6 - revise variance estimates
saveverbose = verbose
saveplottype = plottype
if (verbose eq 5) then stop
profim = fitprof(dataim, spec, x1, x2, noproffit = noproffit, $
                 profdeg = profdeg, varim = varim, inmask = inmask, $
                 profmask = profmask, pthresh = pthresh, bgim = bgim, $
                 Q = Q, v0 = v0, skyvar = skyvar, profres = profres, $
                 fitgauss = fitgauss, fitboxcar = fitboxcar, $
                 boxcarhw = boxcarhw, gotoVect = Pgotovect, $
                 verbose = verbose, errvect = Perrvect, bpct = bpct, $
                 plottype = plottype, difpmask = difpmask, $
                 adjfunc = adjfunc, adjparms = adjparms, $
                 adjoptions = adjoptions)
verbose = saveverbose
plottype = saveplottype

; Horne Steps 6, 7, 8 - Optimally Extract Spectrum
; Also revise variance estimates and mask cosmic rays
if (verbose eq 5) then stop
optspec = extrspec(dataim, profim, varim, v0, Q, x1, x2, bgim = bgim, $
                   inmask = inmask, opvar = opvar, ethresh = ethresh, $
                   exres = exres, skyvar = skyvar, verbose = verbose, $
                   gotovect = Egotovect, errvect = Eerrvect, bpct = bpct, $
                   plottype = plottype, exmask = exmask)

varout = varim                  ; Send to user

; Send summary information to screen based on user request
summaryopt, data, varout, bgim, profim, inmask, bgmask, exmask, $
            difpmask, optspec, stdspec, traceest, $
            berrvect, perrvect, eerrvect, $
            verbose, plottype, adjparms = adjparms, debughead = debughead

return, optspec

end
