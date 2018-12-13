// explore RCTmiss algorithm with different x in P and S models
// IW 7oct2016
adopath ++ C:\Users\ian\Dropbox\ado\rctmiss
pda
set tracedepth 1
set trace off

/*
use H:\missing\sensanal\RCTmiss\QUATRO\QUATRO, clear
keep id sf_mcs* alloc centreid
xi i.centreid
drop if mi(sf_mcsba)
gen sf_mcs_bin =  sf_mcs>40 if !mi(sf_mcs)
*/
use quatro2, clear

rctmiss, pmmdelta(0): logit sf_mcs_bin alloc sf_mcsba _I*

pmm_glm2 logit sf_mcs_bin, ///
	xslist(alloc sf_mcsba _I*) ///
	xplist(alloc sf_mcsba _I*) ///
	b(b) v(V) neff(neff) pstar(pstar) delta(-6) influ(inf)
ereturn post b V
ereturn display
graph box inf*, by(miss)

pmm_glm2 logit sf_mcs_bin, ///
	xslist(alloc) ///
	xplist(alloc sf_mcsba _I*) ///
	b(b) v(V) neff(neff) pstar(pstar) delta(-6) 
ereturn post b V
ereturn display

pmm_glm2 logit sf_mcs_bin, ///
	xslist(alloc) ///
	xplist(alloc) ///
	b(b) v(V) neff(neff) pstar(pstar) delta(-6) 
ereturn post b V
ereturn display

rctmiss, pmmdelta(-6): logit sf_mcs_bin alloc 

