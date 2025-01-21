*! version 1.0.0  21jan2025 pdavidboll
program define _spurtest_i0, rclass
	version 14
	syntax varname(numeric) [if] [in] ,  [q(int 15) nrep(int 100000) latlong ]
	marksample touse
	
	tempname p teststat g_a cv
		
	mata: q = `q'
	mata: nrep = `nrep'
	mata: emat = rnormal(q, nrep, 0, 1)
	
	mata: y = st_data(., "`varlist'", "`touse'")
	mata: s = get_s_matrix("`touse'", "`latlong'")
	mata: distmat = getdistmat_normalized(s)
	mata: test = spatial_i0_test(y, distmat, emat)
	
	mata: st_numscalar("`p'", test.pvalue)
    mata: st_numscalar("`teststat'", test.LR)
	mata: st_numscalar("`g_a'", test.ha_parm)
	mata: st_matrix("`cv'", test.cvalue)
	
	return scalar p=`p'
	return scalar teststat=`teststat'
	return scalar ha_param=`g_a'
	return matrix cv = `cv'
	
	// Display results
    di as text "Spatial I(0) Test Results"
    di as result "---------------------------------------"
    di as result "Test Statistic (LFST) : " %9.4f `teststat'
    di as result "P-value               : " %9.4f `p'
    di as result "---------------------------------------"
	
end