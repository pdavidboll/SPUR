*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix get_sigma_residual(real matrix distmat, real scalar c, real matrix M)
{
	real matrix sigma, sigma_dm
	
	sigma = exp(-c*distmat)
	sigma_dm = M*sigma*M'
	
	return(sigma_dm)
	
}

end