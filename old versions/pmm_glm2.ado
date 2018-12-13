prog def pmm_glm2
version 10
* 7oct2016 - allow XP > XS
* 9nov2016 - add weights

/*******************************************
NOTES
    Avoid -predict, residual- which uses unexpected formulae.
TO DO	
	Make X's optional
23nov2016
	Found successful method for effective sample size
	Still to do: 
		make it faster and not require nxp matrices in Stata
			- best to do this by coding in Mata I think
		tidy up the code - it currently leaves matrices around
		integrate with rctmiss
*******************************************/

// PARSE
syntax anything [if] [in] [iweight/], ///
	XPlist(varlist) XSlist(varlist) delta(string) /// model specification
	[b(string) V(string) neff(string) pstar(string) /// returned values
	noCONStant /// optional model specification
	NEFFValue(real 0) /// optional analysis specification
	keepmat(string) /// optional returned values
	debug display INFluence(string)] // output settings
gettoken cmd y : anything
unabcmd `cmd'
local cmd = r(cmd)
if "`cmd'"=="logistic" local cmd logit
if !inlist("`cmd'","regress","logit","poisson") {
    di as error "Sorry, command `cmd' is not yet supported"
    exit 498
}
if "`debug'"=="" local ifdebug qui
di as text "Method:" _col(26) as result "mean score + joint sandwich variance"
if !missing("`weight'") {
	local wtexp [`weight'=`exp']
	local timesweight *sqrt(`exp')
}
foreach thing in b v neff pstar {
	if "``thing''"=="" local `thing' `thing'
}
if `neffvalue'<0 {
	di as error "neffvalue must be >0"
	exit 198
}

// SET UP
marksample touse, novarlist
tempvar rowmiss id residP predP residS predS ystar offsetvar residPS
tempname bP bCC Vdrop Vmiss Vfull Binv B BPS BPP BSP BSP0 BSS C CSP CSP0 CSS CPP
gen `id'=_n
unab xPlist : `xplist'
local nxPlist : word count `xPlist'
unab xSlist : `xslist'
local nxSlist : word count `xSlist'
`ifdebug' di as text "PM: " as result "`nxPlist'" as text " variables: " as result "`xPlist'"
`ifdebug' di as text "SM: " as result "`nxSlist'" as text " variables: " as result "`xSlist'"
local hascons = ("`constant'"!="noconstant")

// COUNT OBS
qui count if `touse'
local ntot = r(N)
qui count if `touse' & !mi(`y')
local nobs = r(N)

// FIT PATTERN-MIXTURE MODEL (P)
`ifdebug' di as text _new "*** Fitting imputation (pattern-mixture) model ***"
if "`cmd'" != "regress" {
    qui gen `offsetvar' = missing(`y')*`delta' if `touse'
    local offsetopt offset(`offsetvar')
}
`ifdebug' `cmd' `y' `xPlist' if `touse' `wtexp', `offsetopt' `constant'
qui predict `predP' if `touse'
if "`cmd'" == "regress" {
    qui replace `predP' = `predP' + missing(`y')*`delta' if `touse'
	local varP = e(rmse)^2
}
qui gen `residP' = cond(mi(`y'), 0, `y'-`predP') `timesweight' if `touse'
mat `bP' = e(b)
local pP = colsof(`bP')

// SUBSTANTIVE MODEL (S)
if "`cmd'"=="regress" {
    local cmd2 regress
}
else if "`cmd'"=="logit" {
    local cmd2 glm
    local opts family(binomial)
}
else if "`cmd'"=="poisson" {
    local cmd2 glm
    local opts family(poisson)
}
* CC analysis (only for calculating pS and effective sample size)
`ifdebug' di as text _new "*** Fitting CC analysis ***"
`ifdebug' `cmd2' `y' `xSlist' if `touse' `wtexp', `opts' `constant' robust
mat `bCC' = e(b)
local pS = colsof(`bCC')
scalar `pstar' = cond("`cmd'"=="regress",`pS',1)
mat `Vdrop' = e(V)*(e(N)-`pstar')/e(N)

* main analysis
qui gen `ystar' = cond(missing(`y'),`predP',`y') if `touse'
`ifdebug' di as text _new "*** Fitting substantive model ***"
`ifdebug' `cmd2' `ystar' `xSlist' if `touse' `wtexp', `opts' `constant' robust
qui predict `predS' if `touse'
qui gen `residS' = (`ystar' - `predS') `timesweight' if `touse'
mat `b' = e(b)
if "`cmd2'"=="regress" local scale = e(rmse)^2
else if "`cmd2'"=="glm" local scale = e(dispers_p)

`ifdebug' mat list `b', title(b)

// CONSTRUCT C MATRIX
`ifdebug' di as text _new "*** Constructing C matrix ***"
sort `id'
qui gen `residPS'=sqrt(`residP'*`residS') if `touse'
mat opaccum `CSS' = `xSlist' if `touse', group(`id') opvar(`residS') `constant'
mat opaccum `CPP' = `xPlist' if `touse', group(`id') opvar(`residP') `constant'
mat opaccum `CSP0' = `xSlist' `xPlist' if `touse', group(`id') opvar(`residPS') `constant'
local top    = 1
local bottom = `nxSlist'
local left   = `nxSlist' + 1
local right  = `nxSlist' + `nxPlist' + `hascons'
mat `CSP' = `CSP0'[`top'..`bottom',`left'..`right']
if `hascons' mat `CSP' = `CSP' \ `CSP0'[`right',`left'..`right']

mat `C' = (`CSS', `CSP' \ `CSP'', `CPP')
`ifdebug' mat list `C', title(C)

// CONSTRUCT B MATRIX 
`ifdebug' di as text _new "*** Constructing B matrix ***"
tempvar hprimeS hprimeP opSS opSP opPP 
if "`cmd'"=="logit" {
    qui gen `hprimeS' = `predS'*(1-`predS')
    qui gen `hprimeP' = `predP'*(1-`predP')
}
else if "`cmd'"=="regress" {  
    qui gen `hprimeS' = 1
    qui gen `hprimeP' = 1
}
else if "`cmd'"=="poisson" {
    qui gen `hprimeS' = `predS'
    qui gen `hprimeP' = `predP'
}
gen `opSS' = sqrt(`hprimeS') `timesweight'
gen `opSP' = sqrt(`hprimeP') * mi(`y') `timesweight'
gen `opPP' = sqrt(`hprimeP') * !mi(`y') `timesweight'
mat opaccum `BSS' = `xSlist' if `touse', group(`id') opvar(`opSS') `constant'
mat opaccum `BSP0' = `xSlist' `xPlist' if `touse', group(`id') opvar(`opSP') `constant'
mat `BSP' = `BSP0'[`top'..`bottom',`left'..`right']
if `hascons' mat `BSP' = `BSP' \ `BSP0'[`right',`left'..`right']
mat `BSP' = -`BSP'
mat opaccum `BPP' = `xPlist' if `touse', group(`id') opvar(`opPP') `constant'
mat `BPS' = J(`pP',`pS',0)
mat `B' = (`BSS', `BSP' \ `BPS', `BPP') 
`ifdebug' mat list `B', title(B)


// CALCULATE V MATRIX
`ifdebug' di as text _new "*** Constructing V matrix ***"
mat `Binv' = inv(`B')
mat `Vfull' = `Binv' * `C' * `Binv'' 
mat `v' = `Vfull'[1..`pS',1..`pS']
if "`debug'"=="debug" {
    mat l `Binv', title("Binv")
    mat l `Vfull', title("Vfull")
    mat l `v', title("v")
}
mat `Vmiss' = `v'

// EFFECTIVE SAMPLE SIZE
if `neffvalue'==0 {
	`ifdebug' di as text _new "*** Estimating effective sample size ***"
	if `hascons' {
		tempvar one
		gen `one'=1
	}
	if mi("`influence'") tempvar influence
	mkmat `residS', mat(residS)
	mkmat `xSlist' `one', mat(xS)
	mkmat `residP', mat(residP)
	mkmat `xPlist' `one', mat(xP)
	mat U = (diag(residS)*xS, diag(residP)*xP)
	mat A = (I(`pS'), J(`pS',`pP',0))
	mat ABinvU = A*`Binv'*U'
	mat vinv = syminv(`v')
	mat inf = vecdiag(ABinvU'*vinv*ABinvU)'
	svmat inf, name(`influence')
	rename `influence'1 `influence'obs

	mat BSSinv=inv(`BSS')
	mat inf2 = vecdiag(xS*BSSinv'*vinv*BSSinv*xS')'
	svmat inf2, name(`influence')
	rename `influence'1 `influence'full
	qui replace `influence'full =  `influence'full * (`residS'^2 + `scale'*`hprimeP')

	summ `influence'obs if missing(`y'), meanonly
	local wtobs = r(sum)
	summ `influence'full if missing(`y'), meanonly
	local wtfull = r(sum)
	`ifdebug' di as text "Weight for missing value = " as result (`wtobs'/`wtfull')
	local neffvalue = `nobs' + (`wtobs'/`wtfull')*(`ntot'-`nobs')
}
scalar `neff' = `neffvalue'
`ifdebug' di as text "Effective sample size = " as result `neff'

// DF CORRECTION
local factor = `neff' / (`neff' - `pstar')
`ifdebug' di as text "Small-sample correction factor = " as result `factor'
mat `v' = `v' * `factor'

// ROW AND COL NAMES
mat coleq `b' = ""
mat rownames `b' = `y'
mat roweq `v' = ""
mat coleq `v' = ""

// OPTIONALLY SAVE B AND C MATRICES
if "`keepmat'"!="" {
    tokenize "`keepmat'"
    mat `1' = `B'
    di as text "B matrix saved as `1'"
    mat `2' = `C'
    di as text "C matrix saved as `2'"
}

if !missing("`display'") {
	`ifdebug' di as text _new "*** Displaying final results ***"
	ereturn post `b' `v'
	di as text "Effective sample size = " as result `neff'
	ereturn display
}

end
