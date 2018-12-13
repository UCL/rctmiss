// explore RCTmiss algorithm with different x in P and S models
// IW 7oct2016
adopath ++ C:\Users\ian\Dropbox\ado\rctmiss
pda
set tracedepth 1
set trace off

use H:\missing\sensanal\RCTmiss\QUATRO\QUATRO, clear
keep id sf_mcs* alloc centreid
xi i.centreid
drop if mi(sf_mcsba)
gen miss = mi(sf_mcs)
gen sf_mcs_bin =  sf_mcs>40 if !miss

foreach reg in logit regress {
	if "`reg'"=="logit" {
		local deltavals -8(.2)8
		local yvar sf_mcs_bin
	}
	else if "`reg'"=="regress" {
		local deltavals -20(.5)20
		local yvar sf_mcs
	}
	else exit 496
	qui gen `reg'_delta = .
	qui gen `reg'_neff = .
	qui gen `reg'_beta = .
	qui gen `reg'_se = .
	local i 0
	forvalues delta = `deltavals' {
		local ++i
		qui replace `reg'_delta=`delta' in `i'	
		di as input "delta=`delta'"
		pmm_glm3 `reg' `yvar', ///
			xslist(alloc sf_mcsba _I*) ///
			xplist(alloc sf_mcsba _I*) ///
			b(b) v(V) neff(neff) pstar(pstar) delta(`delta') 
		qui replace `reg'_neff=scalar(neff) in `i'
		qui replace `reg'_beta=b[1,1] in `i'
		qui replace `reg'_se=V[1,1] in `i'
		
	}
	foreach thing in neff beta se {
		line `reg'_`thing' `reg'_delta if !mi(`reg'_delta), name(`reg'_`thing', replace)
	}
}
graph combine logit_neff logit_beta logit_se regress_neff regress_beta regress_se
