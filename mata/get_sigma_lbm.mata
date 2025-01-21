*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

real matrix get_sigma_lbm(real matrix distmat)
{	
// Compute the LBM covariance matrix from distmat using the first location as the origin

	real scalar n
	real matrix sigma_lbm
	
	n = rows(distmat)
	
	sigma_lbm = 0.5 * ( J(1,n,distmat[.,1]) + J(n,1,distmat[1,.]) - distmat )
	
	return(sigma_lbm)
}	
end