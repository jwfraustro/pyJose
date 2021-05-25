;+
; NAME:
;      SAMPLESHIFT
;
;
; PURPOSE:
;      Expands the array, shifts the data according to an input
;      vector, and shrinks the data
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Shifting Procedures
;
;
; CALLING SEQUENCE:
;      Result = SAMPLESHIFT (DATA, INX, OUTX)
;
;
; INPUTS:
;      DATA:  The array to be modified
;      INX:   Array of x postions in the reference frame for each
;             pixel in data.
;      OUTX:  The x positon values for the output array. If
;             defined, shift is ignored.  May be modified to reflect the
;             true positions.
;
;
; KEYWORD PARAMETERS:
;      FITSAMPLE:  Set to sample the data when changing shape
;      FITINTERP:  Set to use linear interpolation
;      FITSPLINE:  Set to use a spline when upsizing
;      FITAVERAGE: Set to use nearest neighbor sampling when downsizing
;      FITPOLY:    Set to use polynomail fitting when downsizing
;      FITCUBIC:   Set to use a cubic convolution 
;      DEGCONTR:   Degree of fit to use when downsizing
;
;
; OUTPUTS:
;      Returns the data array that has been expanded or contracted
;      and/or shifted 
;
;
; PROCEDURE:
;      Expand, shift, then contract
;
;
; EXAMPLE:
;      ; create a object frame and trace using the synthetic spectrum package
;      nx = 100
;      ny = 400
;      data = synthspec(traceout = traceout, nx = nx, ny = ny)
;      shift = traceout[0:ny-1] - nx/2  ; adjust output for sampleshift
;      level = 5
;      bigdata = sampleshift(data, level, -shift)
;      erase & tvscl,  bigdata   ; large data
;      redata = sampleshift(bigdata, 1./level, shift*level)
;      erase & tvscl, redata     ; smaller remade data
;
;
; MODIFICATION HISTORY:
;      Written by:    John Dermody, 22 Sep 2003
;
;      $Log: sampleshift.pro,v $
;      Revision 3.0  2004/07/08 20:07:51  jfd28
;      Moving up to 3.0 (reflects this version anyway)
;
;      Revision 2.0  2004/07/08 20:05:01  jfd28
;      Tagging all files to 2.0
;
;      Revision 1.7  2004/07/08 18:47:19  jfd28
;      simplified so it checks for the first fitting option
;      added comments
;
;      Revision 1.5  2004/06/17 14:49:41  jfd28
;      In process between shifting calculation uses
;      interpoltation, and complete rewrite using
;      outx and inx of any sizes / grid
;
;      Revision 1.4  2004/06/02 21:28:45  jfd28
;      allowed degree of poly fit during contraction to be user defined
;      set a user defined variable for nearest neighbor averaging to be used
;
;      Revision 1.3  2004/05/29 20:47:21  jfd28
;      Turned expand and contract into a single variable.
;      Made spline the default interpolation technique
;      new array pixels 0-level interpolate old array pixel 0
;      made polyfit the default contraction technique
;
;+

function SAMPLESHIFT, data, inx, outx, fitsample = fitsample, $
  fitinterp = fitinterp, fitaverage = fitaverage, degcontr = degcontr, $
  fitspline = fitspline, fitpoly = fitpoly, fitcubic = fitcubic


; SET DEFUALTS AND CHECK INPUTS

nx = (size(data))[1]
ny = (size(data))[2]
nxo = (size(outx))[1]
nxi = (size(inx))[1]

if not keyword_defined(degcontr) then degcontr = 2

if degcontr lt 0 then message, "Degree of fit must be ge 0"
if nx ne nxi or ny ne (size(inx))[2] then begin
  message, "inx must be equal in size to data"
endif 
if ny ne (size(outx))[2] then begin
  message, "outx must have equal size y dimension to data"
endif 

; INITIALIZE
; make a large array of differences between the input coordinate and
; output coordinate for each row.
if keyword_set(fitsample) or keyword_set(fitinterp) $
  or keyword_set(fitaverage) then begin
  inxa = rebin(inx, nxi, ny, nxo)
  outxa = rebin(outx, nxo, ny, nxi)
  inxa = transpose(inxa, [2, 1, 0])
  dif = (outxa - inxa)
endif

; SAMPLE
; Just return the pixel with the closest coordiante to whats needed
if keyword_set(fitsample) then begin
  foo = min(abs(dif), rmap, dimension=3)
  rmap = rmap / nxo
  mdata = data[rmap]
endif

; INTERPOLATE
; Interpolate using the closest coordiantes above and below.  Find the
; pixels above and below, and add them together weighted by relative
; distance from the needed coordinate.
if keyword_set(fitinterp) then begin
  foo = min(1. / dif, lmap, dimension=3, subscript_max=hmap)
  lweight = dif[hmap] / (dif[hmap] - dif[lmap]) ; the weight
  lmap = lmap / nxo             ; the map of old pixels to new
  hmap = hmap / nxo 
  mdata = data[hmap] * (1 - lweight)
  mdata += data[lmap] * lweight
  mdata[where(1 - lweight gt 1)] = data[hmap[where(1 - lweight gt 1)]]
  mdata[where(    lweight gt 1)] = data[lmap[where(    lweight gt 1)]]
endif

; SPLINE
; Use the spline function on each row, then repeat the first or last
; pixel so that extrapolation does not occur.
if keyword_set(fitspline) then begin
  mdata = dblarr(nxo, ny)
  for i = 0, ny - 1 do begin &$
    mdata[*,i] = spline(inx[*,i], data[*,i], outx[*,i], 0.1) &$
  endfor
  lez = where(outx le inx[0,     *] ## (fltarr(nxo) + 1), lc)
  gez = where(outx ge inx[nxi-1, *] ## (fltarr(nxo) + 1), gc)
  if lc ne 0 then mdata[lez] = data[0,     lez/nxo]
  if gc ne 0 then mdata[gez] = data[nxi-1, gez/nxo]
endif

; CUBIC
; Use the interpolate function with a cubic convolution of -0.5.
; First calculate the output coordinates when the input coordiante is
; just the pixel value.
if keyword_set(fitcubic) then begin
  mdata = dblarr(nxo, ny)
  for i = 0, ny-1 do begin
    minx = min(inx[*, i])
    maxx = max(inx[*, i])
    xvals = (outx[*, i] - minnx) / (maxx - minx) * (nxi - 1)
    mdata[*, i] = interpolate(data, xvals, cubic = -0.5)
  endfor
endif

; AVERAGE
; Average the closest pixels to the needed coordiante.  Add together
; those pixels with a difference less than range, and divide by the
; number of images
if keyword_set(fitaverage) then begin
  dataa = rebin(data, nxi, ny, nxo)
  dataa = transpose(dataa, [2, 1, 0])
  range = mean(outx[1:nxo-1, *] - outx[0:nxo-2, *]) / 2
  dmask = (dif ge -range) * (dif le range)
  mdata = total(dataa*dmask, 3) / total(dmask, 3)
endif

; POLYNOMIAL FITTING
; Loop through every pixel in the output array.  Grab the closest
; pixels and if there is enough around the needed coordinate fit a
; polynomial to them. Evaluate the polynomial and the output coordinate
if keyword_set(fitpoly) then begin
  mdata = outx - outx
  for i = 0, ny-1 do begin
    range = mean(outx[1:nxo-1, i] - outx[0:nxo-2, i]) / 2
    minx = min(inx[*, i])
    maxx = max(inx[*, i])
    for j = 0, nxo-1 do begin
      x1 = min(where(inx[*, i] ge outx[j, i]-range))
      x2 = max(where(inx[*, i] le outx[j, i]+range))
      if minx gt outx[j, i] then begin
        mdata[j, i] = data[x1, i]
      endif else if maxx lt outx[j, i] then begin
        mdata[j, i] = data[x2, i]
      endif else if x1 eq x2 then begin
        mdata[j, i] = data[x1, i]
      endif else if x1 gt x2 then begin
        dif = inx[x2, i] - inx[x1, i]
        y1 = data[x2, i] * (outx[j, i] - inx[x1, i]) / dif
        y2 = data[x1, i] * (inx[x2, i] - outx[j, i]) / dif 
        mdata[j, i] = y1 + y2
      endif else if (x2-x1+1) le degcontr then begin
        co = linfit(inx[x1:x2, i], data[x1:x2, i], /double)
        mdata[j, i] = co[0] + co[1]*outx[j, i]
      endif else begin
        co = poly_fit(inx[x1:x2, i], data[x1:x2, i], degcontr, $
                      /double, yfit = est, status = status)
        mdata[j, i] = poly(outx[j, i], co)
      endelse
    endfor
  endfor
endif

return, mdata
end
