*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

// Define the struct for Spatial_I1_Test
struct struct_spatial_i0_test_residual {
    real vector LR
    real vector pvalue
    real vector cvalue
    real scalar ha_parm
	real matrix cvalue_mat
	real matrix pvalue_mat
	real vector rho_grid
}

struct struct_spatial_i0_test_residual scalar spatial_i0_test_residual(real matrix Y, real matrix X_in, real matrix distmat, real matrix emat) {
    // Declare all variables
    real scalar q, n, rho, n_y, rho_min, rho_max, n_rho
    real vector lam, rho_grid, sz_vec, cv_vec, LR, pvalue, cvalue
    real matrix sigdm_bm, R, sigdm_rho, sigma, sigma_dm, sigdm_wn, om_rho, om_wn, om_bm
    real matrix om_i0, om_ho, om_ha, ch_omi_ho, ch_omi_ha, ch_om_ho_mat, pvalue_mat, cvalue_mat, ch_om_ho
    real matrix X, P, y_P_ho, y_P_ha, q_P_ho, q_P_ha, y_ho, y_ho_ho, y_ho_ha, q_ho_ho, q_ho_ha, lr_ho, sigdm_ho
    real scalar ha_parm, c, i, ir, j
	real matrix M
	real scalar rho_bm, c_bm
	

    q = rows(emat)
    n = rows(distmat)
	
	// Form M
	M = I(n)-X_in*(invsym(X_in'*X_in))*X_in'
	
	// BM covariance matrix (approximation for demeanded value)
	rho_bm = 0.999
	c_bm = getcbar(rho_bm,distmat)
	sigdm_bm = get_sigma_residual(distmat,c_bm,M)

    // Construct R .. eigenvectors for low-frequency weights
    R = get_R(sigdm_bm, q)

    // Value of om_ho for rhobar = rhobar_max
    rho = 0.001
    c = getcbar(rho, distmat)
    sigdm_rho = get_sigma_residual(distmat, c, M)

    // Matrices used in analysis
    om_rho = R' * sigdm_rho * R
    om_bm = R' * sigdm_bm * R

    // Find value of Ha parameter that yields (roughly 50% power)
    om_i0 = om_rho
    om_ho = om_rho
    ha_parm = get_ha_parm_I0(om_ho, om_i0, om_bm, emat)
    om_ha = om_i0 + ha_parm * om_bm

    // Matrices for carrying out tests
    ch_omi_ho = cholesky(invsym(om_ho))'
    ch_omi_ha = cholesky(invsym(om_ha))'

    // Get LR for data
    n_y = cols(Y)
    LR = J(n_y, 1, .)
    for (i = 1; i <= n_y; i++) {
        X = Y[,i] :- mean(Y[,i])
        P = R' * X
        y_P_ho = ch_omi_ho * P
        y_P_ha = ch_omi_ha * P
        q_P_ho = colsum(y_P_ho:^2)'
        q_P_ha = colsum(y_P_ha:^2)'
        LR[i] = q_P_ho / q_P_ha
    }

    // Construct om_ho for a grid of values of rho
    rho_min = 0.0001
    rho_max = 0.03
    n_rho = 30
    rho_grid = rangen(rho_min, rho_max, n_rho)
    ch_om_ho_mat = asarray_create("real",1)
    for (i = 1; i <= n_rho; i++) {
        om_ho = I(q)
        rho = rho_grid[i]
        if (rho > 0) {
            c = getcbar(rho, distmat)
            sigdm_ho = get_sigma_residual(distmat, c, M)
            om_ho = R' * sigdm_ho * R
        }
		asarray(ch_om_ho_mat, i, cholesky(om_ho)')
    }

    pvalue_mat = J(n_rho, n_y, .)
    sz_vec = (0.01, 0.05, 0.10)
    cvalue_mat = J(n_rho, 3, .)

    for (ir = 1; ir <= n_rho; ir++) {
        ch_om_ho = asarray(ch_om_ho_mat,ir)
        y_ho = ch_om_ho' * emat
        y_ho_ho = ch_omi_ho * y_ho
        y_ho_ha = ch_omi_ha * y_ho
        q_ho_ho = colsum(y_ho_ho:^2)
        q_ho_ha = colsum(y_ho_ha:^2)
        lr_ho = q_ho_ho :/ q_ho_ha
        cv_vec = mm_quantile(lr_ho',1, (1 :- sz_vec'))
        cvalue_mat[ir,] = cv_vec'
        for (j = 1; j <= n_y; j++) {
            pvalue_mat[ir,j] = mean(lr_ho' :> LR[j])
        }
    }
	

    cvalue = colmax(cvalue_mat)'
    pvalue = colmax(pvalue_mat)'

	
	
	// Struct output
     struct struct_spatial_i0_test scalar SP
     SP.LR = LR
     SP.pvalue = pvalue
     SP.cvalue = cvalue
     SP.ha_parm = ha_parm
	 SP.cvalue_mat = cvalue_mat
	 SP.pvalue_mat = pvalue_mat
	 SP.rho_grid = rho_grid

     return(SP)
}


end