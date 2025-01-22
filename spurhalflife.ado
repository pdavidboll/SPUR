*! version 1.0.0  21jan2025 pdavidboll
mata mata clear
mata mata set matastrict on

program define spurhalflife, rclass
	version 14
	syntax varname(numeric) [if] [in] , [q(int 15) nrep(int 100000) Level(real 95) latlong NORMdist]
	marksample touse
	
	// check moremata installed
	cap which moremata.hlp
	if _rc {
		display as error in smcl `"Please install package {it:moremata} from SSC in order to use this command;"' _newline ///
        `"you can do so by clicking this link: {stata "ssc install moremata, replace":auto-install moremata}"'
		exit 199
	}
	
	tempname ci_l ci_u max_dist
	
	if `level' >= 100 | `level' <= 0 {
		di as error "Invalid level."
		exit
	}
		
	mata: q = `q'
	mata: nrep = `nrep'
	mata: emat = rnormal(q, nrep, 0, 1)
	
	mata: y = st_data(., "`varlist'", "`touse'")
	mata: s = get_s_matrix("`touse'", "`latlong'")
	mata: distmat = getdistmat_normalized(s)
	mata: ci = spatial_persistence(y, distmat, emat, `level' / 100)
	
	
	mata: st_numscalar("`ci_l'", ci.hl_ci_lower)
    mata: st_numscalar("`ci_u'", ci.hl_ci_upper)
	
	if `ci_u' >= 100 {
		local ci_u .
	}
	
	mata: max_dist = max(getdistmat(s)) * 3.14159265359 * 6371000.009 * 2 // great circle distance
	mata: st_numscalar("`max_dist'", max_dist)
	if "`normdist'" == "" {
		local ci_l = `ci_l' * `max_dist' 
		local ci_u = `ci_u' * `max_dist'
	}
	
	return scalar ci_l=`ci_l'
	return scalar ci_u=`ci_u'
	return scalar max_dist = `max_dist'
	
	// Display results
	if "`normdist'"!=""{
		local units : di "fractions of maximum distance " %9.4f `max_dist'
	}
	else {
		if "`latlong'"!="" {
			local units "metres"
		}
		else {
			local units "unit of original coordinates"
		}
	}
    di as text "Spatial half life " "`level'" "% confidence interval (in " "`units'" ")"
    di as result "---------------------------------------"
    di as result "Lower bound: " %9.4f `ci_l'
	if `ci_u' != . {
		di as result "Upper bound: " %9.4f `ci_u'
	}
    else {
		di as result "Upper bound:    inf" 
	}
    di as result "---------------------------------------"
	
end