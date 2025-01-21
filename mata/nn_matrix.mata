*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

real matrix nn_matrix(real matrix s)
{
// Compute the NN matrix 
// Input:
//   s: (n,2) matrix of locations
//   latlongflag: 1 if the input is latitute and longitude, 0 if euclidean distance is to be used
	
	external real scalar latlongflag
	real matrix distmat
	real matrix distmat_nozeros
	real colvector rowmins
	real matrix rowmins_mat
	real matrix NN_mat

		
	distmat = getdistmat_normalized(s)
	
	distmat_nozeros = distmat
	distmat_nozeros = distmat_nozeros + (distmat_nozeros :== 0) :* 1e10

	rowmins = rowmin(distmat_nozeros)

	rowmins_mat = distmat_nozeros :== rowmins

	rowmins_mat = rowmins_mat :/ rowsum(rowmins_mat)

	NN_mat = I(rows(distmat)) - rowmins_mat
	
	return(NN_mat)
	
}
end