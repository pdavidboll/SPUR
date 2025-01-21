*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

real matrix getdistmat_normalized(real matrix s)
{	
// Get matrix of distances from location matrix (using scpc function)
// Normalize the matrix so the largest distance is 1

	external real scalar latlongflag
	real matrix distmat
	
	distmat = getdistmat(s)
	distmat = distmat / max(distmat)
	
	return(distmat)
}	
end