*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix get_cluster_matrix(string scalar clustervar, string scalar touse)
{
	real matrix cluster
		
	cluster = st_data(., clustervar, touse)
	
	return(cluster)
}

end