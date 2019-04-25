/*
cscript for rctmiss
13dec2018: updated for v0.12.2 syntax with sens options as suboptions of sens()
	added all-options test
	added basemiss test
28jan-3feb2017: updated (log->exp, sandwich->fullsandwich, listopt->list)
IRW 16dec2016
*/

* confirm location
global root c:\ado\ian\rctmiss
confirm file $root/test/rctmiss_cscript.do
confirm file $root/package/rctmiss.ado
adopath ++ $root/package/
cap log close
log using $root/testlogs/rctmiss_cscript.log, replace

pda
set tracedepth 1
set trace off
set more off
local i 0

which rctmiss, all

use $root/test/QUATRO, clear
keep id sf_mcs* alloc centreid
xi i.centreid
gen miss = mi(sf_mcs)
gen sf_mcs_bin =  sf_mcs>40 if !miss

summ sf_mcsba, meanonly
gen sf_mcsba_imp = cond(mi(sf_mcsba),r(mean),sf_mcsba)

prog def store
	mat b0=e(b)
	mat V0=e(V)
	if e(cmd)=="regress" scalar df0=e(df_r)
	else scalar df0=.
end

prog def comp
	assert e(cmd)=="rctmiss"
	mat b1=e(b)
	mat V1=e(V)
	scalar df1=e(df_r)
	cap assert mreldif(b0,b1)<1E-6
	if _rc {
		di as error "Point estimates differ"
		mat l b0, title(Original b)
		mat l b1, title(New b)
		exit 499
	}
	cap assert mreldif(V0,V1)<1E-6
	if _rc {
		di as error "Variance matrices differ"
		mat l V0, title(Original V)
		mat l V1, title(New V)
		exit 499
	}
	if !mi(scalar(df0)) {
		cap assert reldif(df0,df1)<1E-6
		if _rc {
			di as error "Degrees of freedom differ"
			scalar list df0
			scalar list df1
			exit 499
		}
	}
	di as text "Comparison succeeded at 1E-6 level"
end

prog def dicmd
noi di as input `"`0'"'
`0'
end

// check perfect agreement 

* MAR, reg
dicmd reg sf_mcs alloc sf_mcsba_imp _I*
store
dicmd rctmiss, pmmdelta(0): reg sf_mcs alloc sf_mcsba _I*
comp	

* MAR, reg robust, no x
dicmd reg sf_mcs, robust
store
dicmd rctmiss, pmmdelta(0) fullsandwich: reg sf_mcs
comp
dicmd rctmiss, pmmdelta(0): reg sf_mcs
comp	
dicmd rctmiss, smdelta(0): reg sf_mcs
comp	

* MAR, reg robust
dicmd reg sf_mcs alloc sf_mcsba_imp _I*, robust
store
dicmd rctmiss, pmmdelta(0) fullsandwich: reg sf_mcs alloc sf_mcsba _I*
comp
dicmd rctmiss, pmmdelta(0): reg sf_mcs alloc sf_mcsba _I*, robust
comp	
dicmd rctmiss, smdelta(0): reg sf_mcs alloc sf_mcsba _I*
comp	

* ditto with constant as variable
gen one=1
dicmd reg sf_mcs alloc sf_mcsba_imp _I* one, robust noconst
store
dicmd rctmiss, pmmdelta(0) fullsandwich: reg sf_mcs alloc sf_mcsba _I* one, noconst
comp
dicmd rctmiss, pmmdelta(0): reg sf_mcs alloc sf_mcsba _I* one, robust noconst
comp	
dicmd rctmiss, smdelta(0): reg sf_mcs alloc sf_mcsba _I* one, noconst
comp	

* MAR, reg clustered
dicmd reg sf_mcs alloc sf_mcsba_imp, cluster(centreid)
store
dicmd rctmiss, pmmdelta(0) fullsandwich: reg sf_mcs alloc sf_mcsba, cluster(centreid)
comp	
dicmd rctmiss, pmmdelta(0): reg sf_mcs alloc sf_mcsba, cluster(centreid)
comp	

* MAR, no constant
dicmd reg sf_mcs alloc sf_mcsba_imp, robust noconst
store
dicmd rctmiss, pmmdelta(0) fullsandwich: reg sf_mcs alloc sf_mcsba, noconst
comp	
dicmd rctmiss, pmmdelta(0): reg sf_mcs alloc sf_mcsba, noconst robust
comp	
dicmd rctmiss, smdelta(0): reg sf_mcs alloc sf_mcsba, noconst
comp	

* MAR, logit
dicmd logit sf_mcs_bin alloc sf_mcsba_imp _I*, robust
store
dicmd rctmiss, pmmdelta(0): logit sf_mcs_bin alloc sf_mcsba _I*
comp	

* MAR, logit clustered
dicmd logit sf_mcs_bin alloc sf_mcsba_imp, cluster(centre)
store
dicmd rctmiss, pmmdelta(0): logit sf_mcs_bin alloc sf_mcsba, cluster(centre)
comp	

* M=F, logit
gen sf_mcs_bin_MF = sf_mcs_bin & !miss
dicmd logit sf_mcs_bin_MF alloc sf_mcsba_imp _I*, robust
store
dicmd rctmiss, pmmdelta(-99): logit sf_mcs_bin alloc sf_mcsba _I*
comp

* M=F, logit clustered
dicmd logit sf_mcs_bin_MF alloc sf_mcsba_imp, cluster(centre)
store
dicmd rctmiss, pmmdelta(-99): logit sf_mcs_bin alloc sf_mcsba, cluster(centre)
comp

* M=F with auxiliaries, logit
dicmd logit sf_mcs_bin_MF alloc, robust
store
dicmd rctmiss, pmmdelta(-99) auxiliary(sf_mcsba _I*): logit sf_mcs_bin alloc 
comp

* M=F with auxiliaries, logit clustered
dicmd logit sf_mcs_bin_MF alloc, cluster(centre)
store
dicmd rctmiss, pmmdelta(-99) auxiliary(sf_mcsba): logit sf_mcs_bin alloc, cluster(centre)
comp

// check handling of missing baselines
dicmd rctmiss, pmmdelta(-5): reg sf_mcs alloc sf_mcsba _I*
store
dicmd rctmiss, pmmdelta(-5): reg sf_mcs alloc sf_mcsba_imp _I*
comp


// check syntax more broadly	
foreach cmd in "reg sf_mcs" "logit sf_mcs_bin" {
	foreach delta in -5 0 {
		dicmd rctmiss, pmmdelta(`delta') fullsandwich: `cmd' alloc 
		dicmd rctmiss, pmmdelta(`delta'): `cmd' alloc 
		dicmd rctmiss, aux(sf_mcsba _I*) pmmdelta(`delta'): `cmd' alloc 
		dicmd rctmiss, aux(sf_mcsba _I*) smdelta(`delta'): `cmd' alloc 
	}
	dicmd rctmiss, pmmdelta(-10(2)0) ///
		sens(alloc, note(`cmd') name(PMMsens,replace) list savedta(PMMsens,replace)): ///
		`cmd' alloc 
	dicmd rctmiss, smdelta(-2(.1)0)  ///
		sens(alloc, note(`cmd') name(SMsens,replace) list savedta(SMsens,replace)): ///
		`cmd' alloc 
}

// similar analyses
dicmd logit sf_mcs_bin alloc sf_mcsba _I*
dicmd rctmiss, pmmdelta(1): logit sf_mcs_bin alloc sf_mcsba _I*
dicmd rctmiss, smdelta(1): logit sf_mcs_bin alloc sf_mcsba _I*

dicmd logit sf_mcs_bin alloc 
dicmd rctmiss, pmmdelta(1): logit sf_mcs_bin alloc 
dicmd rctmiss, smdelta(1): logit sf_mcs_bin alloc 
dicmd rctmiss, pmmdelta(1) aux(sf_mcsba _I*): logit sf_mcs_bin alloc 
dicmd rctmiss, smdelta(1) aux(sf_mcsba _I*): logit sf_mcs_bin alloc 

// basemiss
foreach basemiss in "mean" "mim" "mim, min(5)" "mim, min(30)" {
dicmd rctmiss, pmmdelta(-10(2)0) auxil(_Icentreid*) ///
		sens(alloc, name(basemiss`++i')) ///
		basemiss(`basemiss'): ///
		reg sf_mcs alloc sf_mcsba 
}

// many options
dicmd rctmiss, pmmdelta(-10(2)0) auxil(sf_mcsba) fulls ///
		sens(alloc, senstype(unequal) list(sep(3)) clear ///
		stagger(.1) colors(red orange yellow) lwidth(1) lpat(dash) msym(Oh) ciband ///
		note(This is my note) name(manyopts,replace) savedta(manyopts,replace)) ///
		basemiss(mean) eform(My_eform): ///
		reg sf_mcs alloc _Icentreid*

// tidy up
erase SMsens.dta
erase PMMsens.dta
erase manyopts.dta

log close
