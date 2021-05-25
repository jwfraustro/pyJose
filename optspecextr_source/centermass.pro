;+ NAME:
;      CENTERMASS
;
;
; PURPOSE:
;      A function which returns the center of mass of datav
;
;
; CATEGORY:
;      Optimal Spectrum Extraction Package
;          - Vector fitting functions
;
;
; CALLING SEQUENCE:
;      Result = CENTERMASS (XVALS, DATAV, VARV, SPECV)
;
;
; INPUTS:
;      XVALS: The x values for the data
;      DATAV: The sky subtracted data
;      VARV:  The variance
;      SPECV: The spectrum
;
;
; RESTRICTIONS:
;      All vectors must be the same length
;
;
; OUTPUTS:
;      Returns the center of mass of the data vector / spec vector.
;
;
; MODIFICATION HISTORY:		
;		 Written by: John Dermody, 22 Sep 2003
;
;       Thu Dec 18 22:53:00 2003, John Dermody
;              cleaned up a bit
; 
;       $Log: centermass.pro,v $
;       Revision 3.0  2004/07/08 20:07:51  jfd28
;       Moving up to 3.0 (reflects this version anyway)
;
;       Revision 2.0  2004/07/08 20:05:01  jfd28
;       Tagging all files to 2.0
;
;       Revision 1.7  2004/07/08 18:44:39  jfd28
;       added comments
;
;       Revision 1.6  2004/06/03 16:33:09  jfd28
;       chenged print in error message to a message command
;
;       Revision 1.5  2004/06/03 15:48:46  jfd28
;       fixed typo in error check
;
;       Revision 1.4  2004/06/03 15:43:27  jfd28
;       Cleaned up the header
;       Added check for proper inputs
;
;       Revision 1.3  2004/05/27 18:55:09  jfd28
;       testing log
;
;-

function centermass, xvals, datav, varv, specv

nx = n_elements(xvals)

if nx ne n_elements(datav) then begin
  message, "length of xvals not equal to length of datav"
endif

multv = total(datav / specv)

return, total(xvals * datav / specv / multv)

end
