;+ NAME:
;      STDEXTR
;
;
; PURPOSE: 
;      standard box extraction of spectrum, returns spectrum and
;      variance image, and estitmates over bad pixels
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Horne Step Functions
;
;
; CALLING SEQUENCE:
;      Result = STDEXTR(Dataim,Varim,X1,X2,Inmask,Var=var)
;
;
; INPUTS:
;      DATAIM:  sky-subtracted, processed image to extract spectrum
;               from.  Horizontal is pixel position, vertical is 
;               wavelength, assumes no curvature.
;      VARIM:   variance image from processed image
;      X1,X2:   boundaries in x over which spectrum will be summed
;
;
; OPTIONAL KEYWORDS:
;      INMASK:  the cosmic ray mask for the data image
;      STDVAR   variance of extracted spectra, length n where n is the
;               number of wavelengths taken of the image (vertical direction)
;      ADJSPEC: set to lineraly interpolate the data, then extract
;               the standard spectrum.  Only useful if input mask has 
;               bad pixels
;
;
; OUTPUT:
;      returns stdspec, an array of length n containing extracted spectra,
;      where n is the number of wavelengths of the image (vertical
;      direction)
;
;
; RESTRICTIONS: 
;      Dataim, varim and inmask must be the same size.  O <= X1 <= X2
;      <= NX-1
;
;
; PROCEDURE:
;      Horne Step 4 (standard box extraction)
;
;
; MODIFICATION HISTORY:
;       Written by:      Dara Zeehandelaar, 2 March 2003
;
;       $Log&

;       18 Dec 2003, John Dermody:
;                  Cleaned up, added the integration.
;
;       25 March 2003:   removed loops in favor of matrix math
;
; 
;- 

function stdextr, dataim, varim, x1, x2, inmask = inmask, stdvar = stdvar, $
                  adjspec = adjspec

; CHECK INPUTS
ny = (size(dataim))[2]
nx = (size(dataim))[1]
if not keyword_defined(inmask) then inmask = bytarr(nx, ny) + 1

str = ""
str += x1 lt 0 || x1 gt x2                                ? "x1 "     : ""
str += x2 gt nx - 1 || x2 lt x1                           ? "x2 "     : ""
str += (size(varim ))[1] ne nx || (size(varim ))[2] ne ny ? "varim "  : ""
str += (size(inmask))[1] ne nx || (size(inmask))[2] ne ny ? "inmask " : ""
if str ne "" then message, "Poorly formed inputs: " + str

; INTERPOLATE OVER BAD PIXELS
if arg_present(adjspec) then begin
  adjspec = dblarr(ny) 
  for i = 0, ny-1 do begin
    bv = where(inmask[x1:x2, i] eq 0)
    gv = where(inmask[x1:x2, i] eq 1)
    datav = dataim[x1:x2, i]
    if (bv ne [-1]) then datav[bv] = interpol(datav[gv], gv, bv)
    adjspec[i] = total(datav)
  endfor
endif

; STANDARD EXTRACTION
stdspec = total((dataim * inmask)[x1:x2, *], 1)
stdvar  = total((varim  * inmask)[x1:x2, *], 1)


return, stdspec

end
