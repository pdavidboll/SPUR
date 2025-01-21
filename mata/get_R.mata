*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix get_R(real matrix sig, real scalar qmax)
{

	real matrix V, R
	real colvector d, ds
	real colvector ii
	
	symeigensystem(sig, V, d)
	d = d'
	
	ii = order(d,-1)
	d = d[ii]
	V = V[.,ii]
	
	ds = d[1..qmax]
    R = V[., 1..qmax]
	
	return(R)
}

end