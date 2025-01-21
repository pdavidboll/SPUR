*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix make_transform(real matrix s, string scalar transformation, | real scalar radius, real vector cluster)
{
	real matrix H
	external real scalar latlongflag
	
	if (transformation == "nn") {
		H = nn_matrix(s)
	} else if (transformation == "iso") {
		H = iso_matrix(s, radius)
	} else if (transformation == "cluster") {
		H = cluster_matrix(cluster)
	} else if (transformation == "lbmgls") {
		H = lbm_gls_matrix(s)
	} else {
		_error("Invalid transformation.")
	}
	
	return(H)
}

real matrix transform(string scalar varname, real matrix H, string scalar touse, string scalar transformation)
{
 
    real vector y, hy
	real scalar standardflag
	
	y = st_data(., varname, touse)
	
	hy = H * y
	
	if (transformation == "lbmgls") {
		hy = hy:-mean(hy)
	}
	
	return(hy)
 
}

end