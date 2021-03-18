; NAME:
;      TESTOPTEXTR
;
; PURPOSE:
;      Test the installation of the optimal extraction package
;
; CALLING SEQUENCE:
;      TESTOPTEXTR, Precision
;
; INPUTS:
;      PRECISION: Allowable difference between test and control (1.e-6)
;
; PROCEDURE:
;      Run the two test files to make sure the procedue is operating correctly
;
; FUNCTIONS CALLED:     
;      Optspecextr
;
; MODIFICATION HISTORY:
;       Written by:      John Dermody, Sept 16, 2003
;
;       $Log: testoptextr.pro,v $
;       Revision 3.0  2004/07/08 20:07:51  jfd28
;       Moving up to 3.0 (reflects this version anyway)
;
;       Revision 2.0  2004/07/08 20:05:01  jfd28
;       Tagging all files to 2.0
;
;       Revision 1.3  2004/07/02 21:56:34  jfd28
;       Grab locations, complete extaction, and compare to good versions
;
;-

pro testoptextr, precision

on_error, 1

if not keyword_defined(precision) then precision = 1.e-6

locframe1 = file_which('optextr1frame.fits')
locframe2 = file_which('optextr2frame.fits')
locspec1 = file_which('optextr1spec.fits')
locspec2 = file_which('optextr2spec.fits')

; simple frame

frame1 = readfits(locframe1, header1, /silent)
Q =      sxpar(header1, 'EPADU') 
rn =     sxpar(header1, 'RDNOISE') / Q 
objpos = sxpar(header1, 'OBJPOS') 
nx =     sxpar(header1, 'NAXIS1') 
ny =     sxpar(header1, 'NAXIS2') 
seeing = sxpar(header1, 'SEEING') 
trace =  sxpar(header1, 'TRACEAMP') 
tilt =   sxpar(header1, 'TILT') 
tiltadj = (ny / 2) * abs(tan(tilt * !pi / 180)) 
x1 = objpos - trace - 1.5 * seeing - tiltadj 
x2 = objpos + trace + 1.5 * seeing + tiltadj 

dataim = frame1
varim = abs(frame1)/Q + rn^2
opspec1 = optspecextr(dataim, varim, rn, Q, x1, x2, stdspec=stdspec1)

print, "Creating simple test spectrum"
writefits, 'test1spec.fits', opspec1, header1

; default frame

frame2 = readfits(locframe2, header2, /silent)
Q =      sxpar(header2, 'EPADU') 
rn =     sxpar(header2, 'RDNOISE') / Q 
objpos = sxpar(header2, 'OBJPOS') 
nx =     sxpar(header2, 'NAXIS1') 
ny =     sxpar(header2, 'NAXIS2') 
seeing = sxpar(header2, 'SEEING') 
trace =  sxpar(header2, 'TRACEAMP') 
tilt =   sxpar(header2, 'TILT') 
tiltadj = (ny / 2) * abs(tan(tilt * !pi / 180)) 
x1 = objpos - trace - 1.5 * seeing - tiltadj 
x2 = objpos + trace + 1.5 * seeing + tiltadj 

dataim = frame2
varim = abs(frame2)/Q + rn^2
opspec2 = optspecextr(dataim, varim, rn, Q, x1, x2, $
                      /fitgauss, stdspec=stdspec2)

print, "Creating default test data"
writefits, 'test2spec.fits', opspec2, header2

; test them

if (locspec1 ne "")  then begin
  control1 = readfits(locspec1, controlheader1, /silent)  ; restore control1
  if (max(abs(control1 - opspec1)) gt precision) then begin
    print, "FAIL - Simple Frame"
  endif else begin
    print, "PASS - Simple Frame"
  endelse
endif else begin
  message, "Good simple test files not found"
endelse


;defaults

if (locspec2 ne "") then begin
  control2 = readfits(locspec2, controlheader2, /silent)  ; restore control2
  if (max(abs(control2 - opspec2)) gt precision) then begin
    print, "FAIL - Default Frame"
  endif else begin
    print, "PASS - Default Frame"
  endelse
endif else begin
  message, "Good default test files not found"
endelse

end
