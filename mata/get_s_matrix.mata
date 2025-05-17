*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix get_s_matrix(string scalar touse, string scalar latlong)
{
	real matrix s
    real scalar ns
	external real scalar latlongflag
		
	// FIX BUG WITH COORDINATE VARIABLE ORDERING --------------------------
	
	stata("unab slist: s_*")  // create local macro slist containing a list of all s_* variables

	slist = st_local("slist") // import to mata
	slist = tokens(slist)' // tokenize list into vector
	slist_nums = strtoreal(substr(slist,3,.)) // extract numbers for sorting
	slist_order = order(slist_nums,1) // extract numerical order (alphabetical order would be incorrect for d>9)
	slist = slist[slist_order] // re-order elements of slist
	slist_len = length(slist) // number of elements

	for (i=1; i<=slist_len; i++) {
		if (substr(slist[i],3,.) != strofreal(i)){ 	// check that i'th element of slist is equal to "s_i", prohibiting e.g. "s_2 s_3" or "s_1 s_3"
			stata("disp as text"+char(34)+"s_* variables not continuously numbered starting from 1 (s_1, s_2, etc.)"+char(34))
			exit(999)
		}
	}
	
	s=st_data(., slist', touse)		// Import slist variables, now correctly ordered and continuously numbered 
								// expects locations in s_1, s_2, s_3... etc. Dimension d is equal to the number of s_ variables 
	
	// ----------------------------------------------------------------------
	
	ns=rows(s)
	if(sum(rowmissing(s) :== 0) < ns){
		stata("disp as text"+char(34)+"missing value(s) in s_* variable; aborting"+char(34))
		exit(999)
	}
	if(ns<5){
		stata("disp as text"+char(34)+"too few locations found in variables s_1, s_2 etc; aborting"+char(34))
		exit(999)
	}		
	stata("disp as text"+char(34)+"found "+strofreal(rows(s), "%6.0f")+" observations and "+strofreal(cols(s), "%3.0f")+"-dimensional locations in s_*"+char(34))
	latlongflag=latlong!=""
	if(latlongflag==0) {
		stata("disp as text"+char(34)+"using Euclidan norm to compute distance between locations stored in s_*")
	}
	else
	{
		if(cols(s)!=2) {
			stata("disp as text"+char(34)+"with latlong option, there must only be s_1 and s_2 present")
			exit(999)
		}
		else{
			stata("disp as text"+char(34)+"Computing distances on surface of sphere treating s_1 as latitude and s_2 as longitude")
		}
	}
	
	return(s)
}

end