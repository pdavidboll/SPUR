*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

real matrix get_sigma_lbm_dm(real matrix distmat)
{	
// Demean

	real matrix sigma_lbm
	real matrix sigma_lbm_dm
	
	sigma_lbm = get_sigma_lbm(distmat)
	sigma_lbm_dm = demeanmat(sigma_lbm)
	
	return(sigma_lbm_dm)
}	
end