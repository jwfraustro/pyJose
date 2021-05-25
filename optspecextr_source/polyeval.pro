function polyeval,coeffs,xvals

;
; NAME:
;       polyeval
;
; PURPOSE:
;       Takes an array of a set of polynomial coefficients and
;       evaluates a given set of xvalues using them
;
; CALLING SEQUENCE:
;       Result = polyeval(coeffs,xvals)
;
; INPUTS:
;       COEFFS: an array containing a set of polynomial coeffecients
;       XVALS:  locations to evaluate. It's type is the return type
;
; OUTPUTS:
;       Returns an array, dimention n_elements(xvals) containing xvals
;       evalated with the polynomial coefficients
; 
; MODIFICATION HISTORY:
;      Written by:    Dara Zeehandelaar
;      1 Sept 2003:   John Dermody, Setting return type to xvals
 

degree=n_elements(coeffs)

polys = xvals
polys[*] = 0.

;evaluationg polynomial
polys[*]=coeffs[0]
if degree gt 0 then for d=1,degree-1 do polys[*]=polys[*]+coeffs[d]*xvals^d

return,polys

end
