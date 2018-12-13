prog def pmm_glm2
version 10

/*******************************************
NOTES
    Avoid -predict, residual- which uses unexpected formulae.
*******************************************/

// PARSE
syntax anything [if] [in], XPlist(varlist) XSlist(varlist) delta(string) ///
	[debug noCONStant b(string) V(string) neff(string) pstar(string) keepmat(string)]
gettoken cmd y : anything
unabcmd `cmd'
local cmd = r(cmd)
if "`debug'"=="" local ifdebug qui
di as text "Method:" _col(26) as result "mean score + joint sandwich variance"

// SET UP
marksample touse, novarlist
tempvar rowmiss id residP predP residS predS ystar offsetvar residPS
tempname bP bCC Vpred Vdrop Vmiss Vfull Binv B BPS BPP BSP BSP0 BSS C CSP CSP0 CSS CPP
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

// PATTERN-MIXTURE MODEL (P)
`ifdebug' di as text _new "*** Fitting imputation (pattern-mixture) model ***"
if "`cmd'" != "regress" {
    qui gen `offsetvar' = missing(`y')*`delta' if `touse'
    local offsetopt offset(`offsetvar')
}
`ifdebug' `cmd' `y' `xPlist' if `touse', `offsetopt' `constant'
qui predict `predP' if `touse'
if "`cmd'" == "regress" {
    qui replace `predP' = `predP' + missing(`y')*`delta' if `touse'
}
qui gen `residP' = cond(mi(`y'), 0, `y'-`predP') if `touse'
mat `bP' = e(b)
local pP = colsof(`bP')

// SUBSTANTIVE MODEL (S)
if "`cmd'"=="regress" {
    local cmd2 regress
}
else if "`cmd'"=="logit" | "`cmd'"=="logistic" {
    local cmd2 glm
    local opts family(binomial)
}
else if "`cmd'"=="poisson" {
    local cmd2 glm
    local opts family(poisson)
}
else {
    di as error "Sorry, command `cmd' is not yet supported"
    exit 498
}
* CC analysis (only for calculating pS and effective sample size)
`ifdebug' di as text _new "*** Fitting CC analysis ***"
`ifdebug' `cmd2' `y' `xSlist' if `touse', `opts' `constant' robust
mat `bCC' = e(b)
local pS = colsof(`bCC')
scalar `pstar' = cond("`cmd'"=="regress",`pS',1)
mat `Vdrop' = e(V)*(e(N)-`pstar')/e(N)

* main analysis
qui gen `ystar' = cond(missing(`y'),`predP',`y') if `touse'
`ifdebug' di as text _new "*** Fitting substantive model ***"
`ifdebug' `cmd2' `ystar' `xSlist' if `touse', `opts' `constant' robust
mat `Vpred' = e(V)*(e(N)-`pstar')/e(N)
qui predict `predS' if `touse'
qui gen `residS' = (`ystar' - `predS') if `touse'
mat `b' = e(b)

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
if "`cmd'"=="logit" | "`cmd'"=="logistic" {
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
gen `opSS' = sqrt(`hprimeS')
gen `opSP' = sqrt(`hprimeP') * mi(`y')
gen `opPP' = sqrt(`hprimeP') * !mi(`y')
mat opaccum `BSS' = `xSlist' if `touse', group(`id') opvar(`opSS') `constant'
mat opaccum `BSP0' = `xSlist' `xPlist' if `touse', group(`id') opvar(`opSP') `constant'
mat `BSP' = `BSP0'[`top'..`bottom',`left'..`right']
if `hascons' mat `BSP' = `BSP' \ `BSP0'[`right',`left'..`right']
mat `BSP' = -`BSP'
mat opaccum `BPP' = `xPlist' if `touse', group(`id') opvar(`opPP') `constant'
mat `BPS' = J(`pP',`pS',0)
mat `B' = (`BSS', `BSP' \ `BPS', `BPP') 
`ifdebug' mat list `C', title(C)

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
`ifdebug' di as text _new "*** Estimating effective sample size ***"
local Vpred = det(`Vpred')^(-1/`pS')
local Vdrop = det(`Vdrop')^(-1/`pS')
local Vmiss = det(`Vmiss')^(-1/`pS')
local misswt = (`Vmiss'-`Vdrop') / (`Vpred'-`Vdrop')
if `misswt'<0 local misswt=0 
if `misswt'>1 local misswt=1
if "`debug'"=="debug" di "Vpred = `Vpred', Vmiss = `Vmiss', Vdrop = `Vdrop', misswt = `misswt'"
scalar `neff' = `nobs' + `misswt'*(`ntot'-`nobs')

// DF CORRECTION
local factor = `neff' / (`neff' - `pstar')
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

end
