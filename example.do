clear all
mata mata clear

// check if scpc installed
cap which scpc
if _rc {
	display as error in smcl `"Please install package {it:scpc} in order to run this do-file"' _newline ///
        `"which you can find under: {browse "https://github.com/ukmueller/SCPC/"}"'
    exit 199
}

// import variable labels
import excel "./chetty_data_labels.xlsx", sheet("Sheet1") firstrow case(lower) clear
local i 0
foreach v of varlist * {
    local lab`++i' = `v'[1]
}

// import data
import excel "./chetty_data_1.xlsx", sheet("Sheet1") firstrow case(lower) clear

// assign variable labels
local i 0
foreach v of varlist * {
    label variable `v' `"`lab`++i''"'
}

// drop non-contiguous states
drop if state == "HI"
drop if state == "AK"

// rename lat and lon
rename lat s_1
rename lon s_2

// make list of covariates
local myvars "fracblack racseg segpov25 fraccom15 hipc gini incsh1 tsr tsperc hsdrop scind fracrel crimer fracsm fracdiv fracmar loctr colpc coltui colgrad manshare chimp tlfpr migirate migorate fracfor"


// loop over variables
foreach var of varlist am `myvars' {

		local label_`var': variable label `var'

	// i1 test
	spurtest i1 `var', latlong
		local tab_1_`var' = `r(p)'
	
	// i0 test
	spurtest i0 `var', latlong
		local tab_2_`var' = `r(p)'
	
	// half-life
	spurhalflife `var', latlong normdist nrep(10000)
		local tab_3_`var' = `r(ci_l)'
		local tab_4_`var' = `r(ci_u)'
	
	if "`var'"!="am" {
		
		preserve
		// Standardize variables
		qui sum am if !missing(am) & !missing(`var')
		qui replace am = (am - `r(mean)')/`r(sd)' if !missing(am) & !missing(`var')
		qui sum `var' if !missing(am) & !missing(`var')
		qui replace `var' = (`var' - `r(mean)')/`r(sd)' if !missing(am) & !missing(`var')
		
		// Naive OLS
		reg am `var', noconstant vce(cluster state)
			local tab_5_`var' = `e(r2)'
			matrix res = r(table)
			local tab_6_`var' = res[1,1]
			local tab_7_`var' = res[5,1]
			local tab_8_`var' = res[6,1]
		
		// Residual I(1) test
		spurtest i1resid am `var', latlong 
			local tab_9_`var' = `r(p)'
		
		// Residual I()) test (not in table)
		spurtest i0resid am `var', latlong 
		
		// LBMGLS transformation
		qui spurtransform am `var', prefix("h_") latlong replace
		
		// OLS on transformed
		qui reg h_am h_`var', noconstant robust
			local tab_10_`var' = `e(r2)'
		scpc, latlong
			matrix res = e(scpcstats)
			local tab_11_`var' = res[1,1]
			local tab_12_`var' = res[1,5]
			local tab_13_`var' = res[1,6]
		
		restore
		
	}

}


// write to Latex table
texdoc init "chetty_replication.tex", replace force
tex \begin{threeparttable}
tex \begin{tabularx}{1\textwidth}{X @{\extracolsep{3pt}}llllllll@{}}
tex & \multicolumn{3}{c}{Spatial Persistence Statistics} & \multicolumn{5}{c}{Regression of AMI onto Variable} \\
tex \cline{2-4} \cline{5-9}
tex & \multicolumn{2}{c}{\$p\$-Value of Test} & \multicolumn{1}{c}{Half-life} & \multicolumn{3}{c}{Levels} & \multicolumn{2}{c}{LBM-GLS} \\
tex \cline{2-3} \cline{4-4} \cline{5-7} \cline{8-9} 
tex & & & & & \multicolumn{1}{c}{\$\hat{\beta}\$[95\% CI]} & \multicolumn{1}{c}{\$p\$-Value} & & \multicolumn{1}{c}{\$\hat{\beta}\$[95\% CI]} \\
tex Variable & \multicolumn{1}{c}{\$I(1)\$} & \multicolumn{1}{c}{\$I(0)\$} & \multicolumn{1}{c}{95\% CI} & \multicolumn{1}{c}{\$R^2\$} & \multicolumn{1}{c}{Cluster} & \multicolumn{1}{c}{Resid. \$I(1)\$} & \multicolumn{1}{c}{\$R^2\$} & \multicolumn{1}{c}{C-SCPC} \\
tex \hline
foreach var of varlist am `myvars' {
	forval j = 1/13 {
		if "`tab_`j'_`var''"==""{
 			continue
 		}
		local tab_`j'_`var' : di %9.2f `tab_`j'_`var''
		if `tab_`j'_`var'' == . {
 			local tab_`j'_`var' "\phantom{--}\infty"
 		}
		di "`tab_`j'_`var''"
		di substr(strtrim("`tab_`j'_`var''"),1,1)
 		else if substr(strtrim("`tab_`j'_`var''"),1,1) != "-" {
 			local tab_`j'_`var' "\phantom{-} `tab_`j'_`var''"
 		}
	}
	if "`var'"=="am"{
		tex `label_`var'' & \$ `tab_1_`var'' \$ & \$ `tab_2_`var'' \$ & \$ [`tab_3_`var'', `tab_4_`var''] \$ & \multicolumn{1}{c}{NA} & \multicolumn{1}{c}{NA} & \multicolumn{1}{c}{NA} & \multicolumn{1}{c}{NA} & \multicolumn{1}{c}{NA} \\
	}
	else {
	tex `label_`var'' & \$ `tab_1_`var'' \$ & \$ `tab_2_`var'' \$ & \$ [`tab_3_`var'', `tab_4_`var''] \$ & \$ `tab_5_`var'' \$ & \$ `tab_6_`var'' [`tab_7_`var'', `tab_8_`var''] \$ & \$ `tab_9_`var'' \$ & \$ `tab_10_`var'' \$ & \$ `tab_11_`var'' [`tab_12_`var'', `tab_13_`var''] \$ \\
	}
}
tex \hline
tex \end{tabularx}
tex \end{threeparttable}
texdoc close

