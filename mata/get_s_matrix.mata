*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix get_s_matrix(string scalar touse, string scalar latlong)
{
	real matrix s
    real scalar ns
	external real scalar latlongflag
		
	s 	 = st_data(., "s_*", touse)
	
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