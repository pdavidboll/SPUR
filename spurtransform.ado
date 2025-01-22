*! version 1.0.0  21jan2025 pdavidboll
mata mata clear
mata mata set matastrict on

program spurtransform, sortpreserve 
    version 14
 	syntax varlist(numeric) [if] [in] , PREfix(string) [ transformation(string) radius(real -1) clustvar(varname) latlong Replace separately] 
	
	// check moremata installed
	cap which moremata.hlp
	if _rc {
		display as error in smcl `"Please install package {it:moremata} from SSC in order to use this command;"' _newline ///
        `"you can do so by clicking this link: {stata "ssc install moremata, replace":auto-install moremata}"'
		exit 199
	}
	
	// check correct specification of options
	if "`transformation'"=="" {
		local transformation "lbmgls"
	}
	
	if "`transformation'"!="iso" & `radius'!=-1 {
		di as error "Option radius only allowed with transformation(iso)."
		exit
	}
	
	if "`transformation'"=="iso" & `radius'==-1 {
		di as error "Radius missing."
		exit
	}
	
	if "`transformation'"=="iso" & `radius'<=0 {
		di as error "Radius must be positive."
		exit
	}

	if "`transformation'"!="cluster" & "`clustvar'"!="" {
		di as error "Option clustvar only allowed with transformation(cluster)."
		exit
	}
	
	if "`transformation'"=="cluster" & "`clustvar'"=="" {
		di as error "Clustvar missing."
		exit
	}
	
	// make numeric clustervar
	tempvar clustervar
	if "`transformation'"=="cluster" & "`clustvar'"!="" {
		qui egen `clustervar' = group(`clustvar')
	}
	
	// select sample
	if "`separately'"!="" {
		marksample touse, novarlist
	}
	else {
		marksample touse
	}
	
	tempvar touse2
	
 	confirm name `prefix'
 	foreach name of local varlist {
		qui gen `touse2' = `touse'
		if "`separately'"!="" {
			qui replace `touse2' = 0 if missing(`name')
		}
		if "`transformation'"=="cluster" {
			qui replace `touse2' = 0 if missing(`clustervar')
		}
		
		
		capture confirm new variable `prefix'`name'
        if _rc {
            if "`replace'"!="" {
                drop `prefix'`name'
            }
            else {
                di as error "`prefix'`name' already defined or invalid name"
                exit 110
            }
        }
		
		
	   mata: s = get_s_matrix("`touse2'", "`latlong'")
	   
	   if "`transformation'"=="cluster" {
			mata: cluster = get_cluster_matrix("`clustervar'", "`touse2'")
	   }
	   else {
			mata: cluster = 0
	   }
 	   mata: H = make_transform(s, "`transformation'", `radius', cluster)
	   mata: hy = transform("`name'", H, "`touse2'", "`transformation'")
	   quietly mata: st_addvar("double", "`prefix'`name'")
	   mata: st_store(., "`prefix'`name'", "`touse2'", hy)
	   qui drop `touse2'
     }
	
	
end	