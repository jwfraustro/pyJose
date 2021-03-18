;+
; NAME:
;     SUMMARYOPT
;
;
; PURPOSE:
;     Plot summary information from the optimal extraction to the
;     screen so the user can make sure the process is running cleanly
;
;
; CATEGORY:
;     Optimal Spectrum Extraction Package
;         - Miscellaneous
;
;
; CALLING SEQUENCE:
;     optsummary, data, varout, bgim, profim, mask, bgmask, exmask, $
;                 difpmask, optspec, stdspec, traceest, $
;                 debughead = debughead
;
;
; INPUTS:
;      Data:     The original data image
;      Varout:   The calculated variance image
;      Bgim:     The fitted backgroun
;      Profim:   The fitted profile
;      Mask:     The input mask
;      Bgmask:   The mask created during background fitting
;      Exmask:   The mask created during optimal extraction
;      Difpmask: The mask created during profile fitting
;      Optspec:  The optimally extracted spectrum
;      Stdspec:  The standard extracted spectum
;      Traceest: The estimated trace
;      Berrvect: Background fitting error columns
;      Perrvect: Profile fitting error vector
;      Eerrvect: Optimal Extraction error rows
;
;
; KEYWORD PARAMETERS:
;      Debughead: General information about the array being processed
;
;
; OUTPUTS:
;     None
;
;
; PROCEDURE:
;     If plottype[0] is on, send all the data frames to window 11.
;     Then send the standard spectrum and the optimal spectrum to the
;     default IDL window.
;
;
; MODIFICATION HISTORY:
;       Written by John Dermody Jun 14 2004
;
;            $Log: summaryopt.pro,v $
;            Revision 3.0  2004/07/08 20:07:51  jfd28
;            Moving up to 3.0 (reflects this version anyway)
;
;            Revision 2.0  2004/07/08 20:05:01  jfd28
;            Tagging all files to 2.0
;
;            Revision 1.3  2004/07/08 18:49:02  jfd28
;            plottype is an array
;            added check for adjparms first
;
;            Revision 1.2  2004/07/02 18:36:49  jfd28
;            since straightening is done outside, check to see if it
;            is  possible to plot the info from it.
;
;            Revision 1.1  2004/06/25 19:36:56  jfd28
;            initial implimentation.  Plotting data to te screen
;
;-

pro SUMMARYOPT, data, varout, bgim, profim, mask, bgmask, exmask, $
                difpmask, optspec, stdspec, traceest, $
                berrvect, perrvect, eerrvect, $
                verbose, plottype, adjparms = adjparms, debughead = debughead

nx = (size(data))[1]
ny = (size(data))[2]
pnx = (size(difpmask))[1]

rdata = optspec ## (fltarr(nx) + 1) * profim + bgim

brv = strjoin(strtrim(where(Berrvect eq 0), 1), " ")
prv = strjoin(strtrim(where(Perrvect eq 0), 1), " ")
Erv = strjoin(strtrim(where(Eerrvect eq 0), 1), " ")

atagn = keyword_defined(AdjParms) ? Tag_Names(AdjParms) : [""]
admask = in("MASK",    atagn) ? AdjParms.mask : bytarr(nx, ny) + 1

if verbose eq 5 then stop
if plottype[0] then begin
  device, get_screen_size = sz
  device, window_state = ws
  xbuf = 1
  xsze = xbuf+nx
  xwin = (xsze)*9 + pnx
  ybuf = 15

  if (not ws[10]) then begin &$
    window, 10, xsize = xwin, ysize = ny+5*ybuf, title = 'Processed Arrays', $
            retain = 2, xpos = sz[0]-xwin-50, ypos = 100 &$
  endif else begin
    wset, 10
  endelse

  erase
  ts = !D.TABLE_SIZE
  dif = (data - rdata)*exmask*bgmask
  sep = fltarr(1, ny)+ts/2
  nel = n_elements(data)
  st = sort(data) & mnd = data[st[nel*.02]] & mxd = data[st[nel*.98]]
  tvdata = (ts-1) * (data    < mxd > mnd - mnd) / (mxd - mnd)
  st = sort(data) & mnv = varout[st[nel*.02]] & mxv = varout[st[nel*.98]]
  tvvar  = (ts-1) * (varout  < mxv > mnv - mnv) / (mxv - mnv)
  tv, tvdata,           0*xsze, ybuf  &  tv, sep, 1*xsze-1, ybuf, 3
  tvscl, bgim,          1*xsze, ybuf  &  tv, sep, 2*xsze-1, ybuf, 3
  tvscl, profim,        2*xsze, ybuf  &  tv, sep, 3*xsze-1, ybuf, 3
  tv, tvvar,            3*xsze, ybuf  &  tv, sep, 4*xsze-1, ybuf, 3
  tvscl, dif,           4*xsze, ybuf  &  tv, sep, 5*xsze-1, ybuf, 3
  tv, -mask,            5*xsze, ybuf  &  tv, sep, 6*xsze-1, ybuf, 3
  tv, mask - bgmask -1, 6*xsze, ybuf  &  tv, sep, 7*xsze-1, ybuf, 3
  tv, mask - exmask -1, 7*xsze, ybuf  &  tv, sep, 8*xsze-1, ybuf, 3
  tv, mask - admask -1, 8*xsze, ybuf  &  tv, sep, 9*xsze-1, ybuf, 3 
  tv, difpmask - 1,     9*xsze, ybuf  

  xpos = [indgen(9)*xsze + nx/2, nx*9+pnx/2]
  ypos = intarr(10) + 5
  names = ['data', 'bgim', 'profile', 'var', 'dif', 'inmask', 'bgmask', $
           'exmask', '(adjmask)', 'profmask']
  xyouts, xpos, ypos, names, /device, alignment = 0.5

  highs = [string(mxd,         format = '(F7.0)'), $
           string(max(bgim),   format = '(F7.0)'), $
           string(max(profim), format = '(F7.5)'), $
           string(mxv,         format = '(F7.0)'), $
           string(max(dif),    format = '(F7.0)')]
  lows =  [string(mnd,         format = '(F7.0)'), $
           string(min(bgim),   format = '(F7.0)'), $
           string(min(profim), format = '(F7.5)'), $
           string(mnv,         format = '(F7.0)'), $
           string(min(dif),    format = '(F7.0)')]
  foo = where(mask eq 0, incount)
  foo = where(mask - bgmask gt 0, bgcount)
  foo = where(mask - exmask gt 0, excount)
  foo = where(mask - admask gt 0, adcount)
  foo = where(difpmask, dfcount)
  nbads = [string(incount, format = '(I7)'), $
           string(bgcount, format = '(I7)'), $
           string(excount, format = '(I7)'), $
           string(adcount, format = '(I7)'), $
           string(dfcount, format = '(I7)')]
  bnumbers = [lows, nbads]
  ypos = intarr(10) + 5 + ybuf + ny
  xyouts, xpos, ypos, bnumbers, /device, alignment = 0.5
  tnumbers = [highs]
  ypos = intarr(5) + 5 + 2*ybuf + ny
  xyouts, xpos[0:4], ypos, tnumbers, /device, alignment = 0.5

  xyouts, nx * 5/2,  5 + 3*ybuf + ny, "Image Ranges", /device, alignment = 0.5
  xyouts, nx * 19/2, 5 + 3*ybuf + ny, "Number of Bad Pixels", /device, $
          alignment = 0.5
  if keyword_defined(debughead) then begin
    xyouts, xwin/2, 5 + 4*ybuf + ny, debughead, /device, alignment = 0.5
  endif 
  wset, 0
  yvals = indgen(ny)
  if keyword_defined(adjparms) then begin
    !p.multi = [0, 2, 2, 1, 1]
    if in("TRACEEST", Tag_Names(adjparms)) then begin
      plot, adjparms.traceest, yvals, title = 'Estimated Trace', $
            xtitle = 'Pixels right of X1', /xstyle, /ystyle
    endif
    if in("WIDTHEST", Tag_Names(adjparms)) then begin
      plot, adjparms.widthest, yvals, title = 'Estimated Width', $
            /xstyle, /ystyle
    endif
  endif else begin
    !p.multi = [0, 1, 2, 1, 1]
  endelse
  plot, optspec, title = 'Optimally Extracted Spectrum'
  plot, stdspec, title = 'Standard Extraction Spectrum'
  !p.multi = 0
endif

initf = "Eliminated half the pixels at "
fsummary = initf
esummary = ""
if (brv ne "-1") then fsummary = fsummary + "rows " + brv + " in BG fit, "
if (prv ne "-1") then fsummary = fsummary + "cols " + prv + " in Prof fit. " 
if (erv ne "-1") then esummary = "Unable to extract spectrum at rows " + erv

if (fsummary eq initf) then begin
  summary = esummary
endif else begin
  summary = fsummary + esummary
endelse
if ((summary ne "") and (verbose gt 0)) then print, summary


end
