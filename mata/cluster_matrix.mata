*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

real matrix cluster_matrix(real colvector cluster)
{
// Compute the cluster matrix 
// Input:
//   cluster: vector of cluster ids
	
	real scalar n
	real matrix clust_mat
	
	n = rows(cluster)
	
	stata("disp as text"+char(34)+"Number of observations: "+strofreal(n, "%9.0g")+", number of clusters "+strofreal(rows(uniqrows(cluster)), "%9.0g"))
		
	clust_mat = J(1,n,cluster) :== J(n,1,cluster')
	
	clust_mat = clust_mat :/ rowsum(clust_mat)
				
	clust_mat = I(n) - clust_mat
	
	return(clust_mat)
	
}
end