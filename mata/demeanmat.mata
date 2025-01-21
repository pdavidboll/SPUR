*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix demeanmat(real matrix mat)
// demeans matrix
{
	mat=mat:-(rowsum(mat)/cols(mat))
	mat=mat:-(colsum(mat)/rows(mat))
	return(mat)
}

end