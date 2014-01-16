; $Id: //depot/idl/IDL_71/idldir/lib/amoeba.pro#1 $
; Copyright (c) 1996-2009, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
; IDL(r) is a registered trademark of Exelis Inc. and IDL(r) 
; software is copyrighted by Exelis Visual Information Solutions, Inc.

; Modified by R.D. Hanna, 2009
; 
; Copied from above location in IDL 7.1 distribution, changing scale 
;   keyword so it will accept 2 bounds for each dimension and inserted
;   a slope check
;

Function BHC_slope_check, ptry
  
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
        LUN, show, y0, yN, bhcParams, matMeans, reconutct
  
  ptry_orig = ptry
  ndims = n_elements(ptry)
    
  ; add endpoints to ptry and adjust ndims since ends fixed
  ptry = [y0, ptry, yN]
  ndims += 2
  
  ptest = BHC_Construct_Transform(ptry)
  x = *(bhcParams.x)
  
  ; check slopes and re-adjust if necessary
  for i=ndims-2,1,-1 do begin
    slopes = (ptest[i+1] - ptest[0:i])/(x[i+1]-x[0:i])    
    maxSlope = Max(slopes,maxPos)
    if (maxPos NE i) then begin
      ptest_prev = ptest
      ptest[i] = ptest[i+1] - maxSlope*(x[i+1]-x[i]) 
    endif
  endfor

  ptry = BHC_Construct_Transform(ptest,/INVERSE)
  
  return, ptry[1:ndims-2]
END


Function BHC_amotry, p, y, psum, func, ihi, fac, bounds
; Extrapolates by a factor fac through the face of the simplex, across
; from the high point, tries it and replaces the high point if the new
; point is better.
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
        LUN, show, y0, yN, bhcParams, matMeans, reconutct
           
  compile_opt hidden
  fac1 = (1.0 - fac) / n_elements(psum)
  fac2 = fac1  - fac
  ptry = psum * fac1 - p[*,ihi] * fac2
  
  ; check boundaries of simplex
  ; lower bound changes increases i*deltaY, but both upper and lower bounds stay within
  ;     the global bounds defined.
  ndims = n_elements(ptry)
  if n_elements(bounds) ne 0 then $
    for i=0,ndims-1,1 do begin
      if ptry[i] lt bounds[0,i] then ptry[i] = bounds[0,i] $
        else if ptry[i] gt bounds[1,i] then ptry[i] = bounds[1,i]
    endfor 
 
  ; check slopes
   if (bhcParams.pos_2nd_der) then ptry = BHC_slope_check(ptry)
  
  ytry = call_function(func, ptry)  ;Eval fcn at trial point
  if ytry lt y[ihi] then begin    ;If its better than highest, replace highest
    y[ihi] = ytry
    psum = psum + ptry - p[*,ihi]
    p[0,ihi] = ptry
    endif
  return, ytry
end



Function BHC_Amoeba, ftol, FUNCTION_NAME=func, FUNCTION_VALUE=y, $
	NCALLS = ncalls, NMAX = nmax, P0 = p0, SCALE=scale, SIMPLEX=p
; The Numerical Recipes procedure Amoeba, with some embellishments.
;
; Reference: Numerical Recipes, 2nd Edition, Page 411.
;  P = fltarr(ndim, ndim+1)
;+
; NAME:
;	AMOEBA
;
; PURPOSE:
;	Multidimensional minimization of a function FUNC(X), where
;	X is an N-dimensional vector, using the downhill simplex
;	method of Nelder and Mead, 1965, Computer Journal, Vol 7, pp 308-313.
;
;	This routine is based on the AMOEBA routine, Numerical
;	Recipes in C: The Art of Scientific Computing (Second Edition), Page
;	411, and is used by permission.
;
; CATEGORY:
;	Function minimization/maximization. Simplex method.
;
; CALLING SEQUENCE:
;	Result = AMOEBA(Ftol, ....)
; INPUTS:
;    FTOL:  the fractional tolerance to be achieved in the function
;	value.  e.g. the fractional decrease in the function value in the
;	terminating step.  This should never be less than the
;	machine's single or double precision.
; KEYWORD PARAMETERS:
;    FUNCTION_NAME: a string containing the name of the function to
;	be minimized.  If omitted, the function FUNC is minimized.
;	This function must accept an Ndim vector as its only parameter and
;	return a scalar single or double precision floating point value as its
;	result. 
;    FUNCTION_VALUE: (output) on exit, an Ndim+1 element vector
;	containing the function values at the simplex points.  The first
;	element contains the function minimum. 
;    NCALLS: (output) the of times the function was evaluated. 
;    NMAX: the maximum number of function evaluations allowed
;	before terminating.  Default = 5000.
;    P0: Initial starting point, an Ndim element vector.  The starting
;	point must be specified using either the keyword SIMPLEX, or P0 and
;	SCALE.  P0 may be either single or double precision floating.
;	For example, in a 3-dimensional problem, if the initial guess
;	is the point [0,0,0], and it is known that the function's
;	minimum value occurs in the interval: -10 <
;	X(0) < 10, -100 < X(1) < 100, -200 < X(2) < 200, specify: P0=[0,0,0],
;	SCALE=[10, 100, 200]. 
;    SCALE: a scalar or Ndim element vector contaiing the problem's
;	characteristic length scale for each dimension.
;	SCALE is used with P0 to form an initial (Ndim+1) point simplex.
;	If all dimensions have the same	scale, specify a scalar.
;    SIMPLEX: (output and/or optional input) On input, if P0 and SCALE
;	are not set, SIMPLEX contains the Ndim+1 vertices, each of
;	Ndim elements, of starting simplex, in either single or double
;	precision floating point, in an (Ndim, Ndim+1) array. On output,
;	SIMPLEX contains the simplex, of dimensions (Ndim, Ndim+1), enclosing
;	the function minimum.  The first point, Simplex(*,0), corresponds to
;	the function's minimum.
;
; OUTPUTS:
;   Result: If the minimum is found, an Ndim vector, corresponding to
;	the Function's minimum value is returned.  If a function minimum
;	within the given tolerance, is NOT found in the given number of
;	evaluations, a scalar value of -1 is returned.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; PROCEDURE:
;	This procedure implements the Simplex method, described in
;	Numerical Recipes, Section 10.4.  See also the POWELL procedure.
;
;	Advantages:  requires only function evaluations, not
;	derivatives, may be more reliable than the POWELL method.
;	Disadvantages: not as efficient as Powell's method, and usually
;	requires more function evaluations.
;
;	Results are performed in the mode (single or double precision)
;	returned by the user-supplied function.  The mode of the inputs P0,
;	SCALE, or SIMPLEX, should match that returned by the function. The
;	mode of the input vector supplied to the user-written function, is
;	determined by P0, SCALE, or SIMPLEX.
;
; EXAMPLE:
;	Use Amoeba to find the slope and intercept of a straight line fitting
;	a given set of points minimizing the maximum error:
;
;	The function to be minimized returns the maximum error,
;	given p(0) = intercept, and p(1) = slope:
; FUNCTION FUNC, p
; COMMON FUNC_XY, x, y
; RETURN, MAX(ABS(y - (p(0) + p(1) * x)))
; END
;
;	Put the data points into a common block so they are accessible to the
;	function: 
; COMMON FUNC_XY, x, y
;	Define the data points:
;   x = findgen(17)*5
;   y = [ 12.0,  24.3,  39.6,  51.0,  66.5,  78.4,  92.7, 107.8, 120.0, $
;        135.5, 147.5, 161.0, 175.4, 187.4, 202.5, 215.4, 229.9]
;
;	Call the function.  Fractional tolerance = 1 part in 10^5, 
;	Initial guess = [0,0], and the minimum should be found within
;	a distance of 100 of that point: 
;   r = AMOEBA(1.0e-5, SCALE=1.0e2, P0 = [0, 0], FUNCTION_VALUE=fval)
;
;	Check for convergence:
;   if n_elements(r) eq 1 then message,'AMOEBA failed to converge'
;	Print results.
;   print, 'Intercept, Slope:', r, 'Function value (max error): ', fval(0)
;Intercept, Slope:      11.4100      2.72800
;Function value:       1.33000
;
; MODIFICATION HISTORY:
;	DMS, May, 1996.	Written.
;	
;	Romy D. Hanna - 2009 -  changed scale 
;   keyword so it will accept 2 bounds for each dimension and inserted a 
;   slope check 
;-
COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
        LUN, show, y0, yN, bhcParams, matMeans, reconutct
        
if keyword_set(scale) then begin  ;If set, then p0 is initial starting pnt
   if (n_elements(scale) eq 1) then begin     ; scalar scale
      ndim = n_elements(p0)*2
      bounds = replicate(1.0, 2, ndim)
      for i=0, ndim-1 do begin
        bounds[0,i] = -1.0*scale[0]
        bounds[1,i] = scale[0]
      endfor
   endif else begin
      temp = scale
      ndim = n_elements(p0)
      bounds = replicate(1.0, 2, ndim)
      scale = replicate(1.0, ndim)
      for i=0, ndim-1 do begin
        scale[i] = p0[i] + ((p0[i]-temp[2*i]) GT (temp[2*i+1]-p0[i]) ? -1 : 1) * 0.5*(temp[2*i+1]-temp[2*i])
        bounds[0,i] = temp[2*i]
        bounds[1,i] = temp[2*i+1]
      endfor
   endelse
   ndim = n_elements(p0)
   p = p0 # replicate(1.0, ndim+1)
   ; Actually seems to result in a better answer if allowed to start
   ;    out of bounds a bit 
   for i=0, ndim-1 do p[i,i+1] = scale[i < (n_elements(scale)-1)]
endif

s = size(p)
if s[0] ne 2 then message, 'Either (SCALE,P0) or SIMPLEX must be initialized'
ndim = s[1]			;Dimensionality of simplex
mpts = ndim+1			;# of points in simplex
if n_elements(func) eq 0 then func = 'FUNC'
if n_elements(nmax) eq 0 then nmax = 5000L

; Make sure all initial functions are legal
if(bhcParams.pos_2nd_der) then begin 
  for i=0,ndim do begin
    ptry = BHC_slope_check(p[*,i])
    p[*,i] = ptry
  endfor
endif
  
y = replicate(call_function(func, p[*,0]), mpts)  ;Init Y to proper type
for i=1, ndim do y[i] = call_function(func, p[*,i]) ;Fill in rest of the vals
ncalls = 0L
psum = total(p,2)

while ncalls le nmax do begin		;Each iteration
  s = sort(y)
  ilo = s[0]		;Lowest point
  ihi = s[ndim]		;Highest point
  inhi = s[ndim-1]	;Next highest point
  d = abs(y[ihi]) + abs(y[ilo]) ;Denominator = interval
  if d ne 0.0 then rtol = 2.0 * abs(y[ihi]-y[ilo])/d $
  else rtol = ftol / 2.         ;Terminate if interval is 0

  if rtol lt ftol then begin ;Done?
      t = y[0] & y[0] = y[ilo] & y[ilo] = t ;Sort so fcn min is 0th elem
      t = p[*,ilo] & p[*,ilo] = p[*,0] & p[*,0] = t
      return, t                 ;params for fcn min
  endif
    
  ncalls = ncalls + 2
  ytry = BHC_amotry(p, y, psum, func, ihi, -1.0, bounds)
  if ytry le y[ilo] then ytry = BHC_amotry(p,y,psum, func, ihi, 2.0, bounds) $
  else if ytry ge y[inhi] then begin
    ysave = y[ihi]
    ytry = BHC_amotry(p,y,psum,func, ihi, 0.5, bounds)
    if ytry ge ysave then begin
      for i=0, ndim do if i ne ilo then begin
	       psum = 0.5 * (p[*,i] + p[*,ilo])
         p[*,i] = psum
         y[i] = call_function(func, psum)
         endif
         ncalls = ncalls + ndim
         psum = total(p,2)
         endif		;ytry ge ysave
    endif else ncalls = ncalls  - 1
  endwhile
  return, -1		;Here, the function failed to converge.
end
