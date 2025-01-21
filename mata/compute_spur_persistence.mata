*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

// Define the struct for Spatial_I1_Test
struct struct_spatial_persistence {
    real scalar hl_ci_lower
	real scalar hl_ci_upper
}

struct struct_spatial_persistence scalar spatial_persistence(real matrix Z, real matrix distmat, real matrix emat, real scalar level) {
	
	real scalar n_hl, n_hl_ha, n_hl_grid_ho, ci_l, ci_u
	real matrix hl_grid_ho, hl_grid_ha, c_grid_ho, c_grid_ha, pv_mat, ii, hl
	
	n_hl = 100
	hl_grid_ho = rangen(0.001,1,n_hl)
	hl_grid_ho = hl_grid_ho \ rangen(1.01,3,30)
	hl_grid_ho = hl_grid_ho \ 100
	n_hl_grid_ho = length(hl_grid_ho)
	n_hl_ha = 50
	hl_grid_ha = rangen(0.001,1.0,n_hl_ha)
	
	c_grid_ho = -log(0.5) :/ hl_grid_ho
	c_grid_ha = -log(0.5) :/ hl_grid_ha
	
	pv_mat=c_ci(Z,distmat,emat,c_grid_ho,c_grid_ha)
    ii = pv_mat :> 1 - level
    hl = select(hl_grid_ho,ii)
    ci_l = min(hl)
    ci_u = max(hl)
	
	
	// Struct output
     struct struct_spatial_persistence scalar SP
     SP.hl_ci_lower = ci_l
     SP.hl_ci_upper = ci_u

     return(SP)
	
	
	
}

real matrix c_ci(real matrix Y, real matrix distmat, real matrix emat, 
                 real vector c_grid_ho, real vector c_grid_ha) {
				 	
    // Variable declarations
    real scalar q, n, n_c_ho, n_c_ha, i, j
    real scalar rho_bm, c, const_den, den_ho_X, den_ha_avg_X, lr_X
    real matrix sigdm_bm, R, lam, om, sigdm, omi
    real matrix ch_om_mat, ch_omi_mat, const_den_mat
    real matrix ch_omi_mat_ha, const_den_mat_ha
    real matrix sz_vec, cv_mat, pv_mat, X, den_ho_X_mat, den_ha_avg_X_mat
    real matrix e, xc, Xc, den_ho, den_ha_mat, den_ha_mat_X, lr, cvalue

    // Initialization
    q = rows(emat)
    n = rows(distmat)
    n_c_ho = length(c_grid_ho)
    n_c_ha = length(c_grid_ha)
    

    // BM covariance matrix
	rho_bm = 0.999
    c_bm = getcbar(rho_bm, distmat)
    sigdm_bm = get_sigma_dm(distmat, c_bm)

    // Construct R and eigenvectors
    R = get_R(sigdm_bm, q)

    // Compute Omega matrices for c_grid_ho
    ch_om_mat = asarray_create("real",1)
    ch_omi_mat = asarray_create("real",1)
    const_den_mat = J(n_c_ho, 1, .)

    for (i = 1; i <= n_c_ho; i++) {
        c = c_grid_ho[i]
        //om = I(q)
        sigdm = get_sigma_dm(distmat, c)
        om = R' * sigdm * R
		asarray(ch_om_mat, i, cholesky(om)')
        omi = invsym(om)
		asarray(ch_omi_mat, i, cholesky(omi)')
        const_den_mat[i,1] = sqrt(det(omi)) * 0.5 * gamma(q / 2) / (pi()^(q / 2))
    }

    // Compute Omega matrices for c_grid_ha
	ch_omi_mat_ha = asarray_create("real",1)
    const_den_mat_ha = J(n_c_ha, 1, .)

    for (i = 1; i <= n_c_ha; i++) {
        c = c_grid_ha[i]
        //om = I(q)
        sigdm = get_sigma_dm(distmat, c)
        om = R' * sigdm * R
        omi = invsym(om)
		asarray(ch_omi_mat_ha, i, cholesky(omi)')
        const_den_mat_ha[i,1] = sqrt(det(omi)) * 0.5 * gamma(q / 2) / (pi()^(q / 2))
    }

    // Construct critical values and p-value matrices
	sz_vec = (0.10, 0.05, 0.01)
    cv_mat = J(length(sz_vec), n_c_ho, .)
    pv_mat = J(n_c_ho, 1, .)
    X = R' * Y
    den_ho_X_mat = J(n_c_ho, 1, .)
    den_ha_avg_X_mat = J(n_c_ho, 1, .)

    for (i = 1; i <= n_c_ho; i++) {
		ch_null = asarray(ch_om_mat,i)
        e = ch_null' * emat
        ch_omi = asarray(ch_omi_mat,i)
        const_den = const_den_mat[i,1]
        xc = ch_omi * e
        den_ho = const_den * ((colsum(xc:^2)) :^ (-q / 2))'
        den_ha_mat = J(cols(emat), n_c_ha, .)
		
        // Compute for data
        Xc = ch_omi * X
        den_ho_X = const_den * ((colsum(Xc:^2)) :^ (-q / 2))'
        den_ha_mat_X = J(n_c_ha, 1, .)

        for (j = 1; j <= n_c_ha; j++) {
			ch_omi = asarray(ch_omi_mat_ha,j)
            const_den = const_den_mat_ha[j,1]
            xc = ch_omi * e
            den_ha_mat[.,j] = const_den * ((colsum(xc:^2)) :^ (-q / 2))'
            Xc = ch_omi * X
            den_ha_mat_X[j,1] = const_den * ((colsum(Xc:^2)) :^ (-q / 2))
        }
	
        den_ha_avg = mean(den_ha_mat')'
        lr = den_ha_avg :/ den_ho
		cvalue = mm_quantile(lr,1, (1 :- sz_vec'))
        cv_mat[.,i] = cvalue
        den_ha_avg_X = mean(den_ha_mat_X)'
        lr_X = den_ha_avg_X / den_ho_X
        pv_mat[i,1] = mean(lr :> lr_X)
        den_ha_avg_X_mat[i,1] = den_ha_avg_X
        den_ho_X_mat[i,1] = den_ho_X
		
		
    }

    return(pv_mat)
}




end