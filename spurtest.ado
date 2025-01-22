*! version 1.0.0  21jan2025 pdavidboll
mata mata clear
mata mata set matastrict on

program spurtest, prop(twopart)
	version 14
	
	// check moremata installed
	cap which moremata.hlp
	if _rc {
		display as error in smcl `"Please install package {it:moremata} from SSC in order to use this command;"' _newline ///
        `"you can do so by clicking this link: {stata "ssc install moremata, replace":auto-install moremata}"'
		exit 199
	}
	
	local alltests i1 i0 i1resid i0resid

	local test : word 1 of `0'
	
	local 0 : subinstr local 0 "`test'" ""
	
	// if user specifies test but no variable, the comma stays with
	// test, issuing wrong error message
	local test : subinstr local test "," "" , count(local comma)
	if `comma' {
		local 0 , `0'
	}
	
	local valid : list posof "`test'" in alltests
	if !`valid' {
		// remove comma before options
		local test : subinstr local test "," ""
		di as error "test `test' not allowed"
		exit 198
	}
	
	syntax varlist(numeric) [if] [in], [*]
	
	// run test
	_spurtest_`test' `varlist' `if' `in', `options'

end