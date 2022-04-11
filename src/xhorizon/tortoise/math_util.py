
import numpy as np




def pos_real_roots(p,imtol=1e-6):
	"""
	Return the sorted positive real roots of a polynomial.
	p is an array of coefficients of the polynomial
		p[n] * x^0 + p[n-1] * x^1 + ... + p[0] * x^n
	(aka highest order first).
	"""
	## get roots of polynomial
	r = np.roots(p)
	## find real values to within tolerance
	mask = ( np.abs(r.imag) < imtol )
	r = r[mask]
	r = np.real(r)
	## get positive values and sort
	r = np.sort( r[r>0] )
	## return
	return r








