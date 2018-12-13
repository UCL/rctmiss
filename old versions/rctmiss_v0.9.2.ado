*! version 0.9.3 IRW 7jan2017
/*******************************************************************************
TO DO:
	better error code for (or even handle) collinear covariates
	delete pmm_glm in due course
	should we have neff for SM? e.g. sum(w_i^2) / ave(w_i)^2 ?
	give better error message when pmm_glm3 run without missing values

HISTORY
version 0.9.3 7jan2017
	speeded up neff calculation by cleverer matrix coding
version 0.9.2 3jan2017
    actually SM was correct before: auxiliary mustn't be in numerator of SW
version 0.9.1 29dec2016
    corrected error in SM: stabilised weights partly ignored auxiliary
version 0.9 16dec2016
	new auxiliary option and new call to pmm_glm3 (line 212)
	non-integer dof used: now CI has half-width exactly =invttail(e(neff)-e(pstar),.025)*_se[alloc]
version 0.8.1   29aug2013
    bug fix - arm labelling in legend of sens graph was wrong when rand wasn't 0/1
version 0.8   7jan2013 -- ON BSU WEBSITE
    Corrected dof errors:
        - pstar returned wrongly by pmm_glm
        - dof wrongly set in non-regress sensitivity analysis 
    Situation now is that CIs use 
        - t-distribution after regress 
        - Normal distribution after other regcmds
    Sytnax changed to senstype(equal|unequal|all) though the ambiguous both|one still work.
    Help file improved.
version 0.7.1   25may2012
    neff, pstar returned as scalars by subroutines
    new ereturn of n*, pstar, method
    nmissmin option becomes min suboption of basemiss (documented)
    pmm_reg: 
        new calculation of neff, based on ratio of small to large sample variance
        now the two variances are added with no scaling
    output formatted in columns
    testscript - OK
version 0.7   15may2012
    Getting CIs right too, by posting dof
    Trying to return neff and pstar - so far done only for mean score linear regression?
version 0.6.5    8may2012
    New sandwich option forces use of sandwich variance
version 0.6.4   30jan2012 -- ON BSU WEBSITE
    level() enabled (either as prefix option or as regression command option)
version 0.6.3   13sep2011
    new undocumented mmstore() option stores results of missingness model
version 0.6.2   30aug2011
    change IM_exp to exp|numlist, [log base(#)]
version 0.6.1   30aug2011
    cistyle(line) becomes ciband
    saving() becomes savedta() so that saving() applies to graph
    replace option moved to suboption of savedta()
version 0.6   26aug2011
    horizontal axis is on same scale as requested (can change using xscale(log))
    delta(log 0(0.1)1) now works because log(0) is taken as -999
    better error capture for misspecified pmmdelta() or smdelta()
version 0.5   27jun2011
version 0.4   3may2011  renamed mnar_mml.ado as rctmiss_smlik.ado + added to package; default changed from smdelta(0) to pmmdelta(0); nosw option; effective sample size added & used in small-sample correction; dfcorrection removed; leaves no matrices in memory
version 0.3.2   12dec2010  only tidied up comments
version 0.3.1   3nov2010  clear works with nograph
version 0.3    30jun2010  new basemiss and nmissmin options; noconstant works properly; doesn't replay after sens() option; new options cistyle(line) savewt() senstype(); log tidied up
version 0.2.4  24jun2010  eform works properly; allow log as first element of smdelta or pmmdelta (is everything labeled correctly? no, output data & list aren't) ***
version 0.2.3  14jun2010  deleted unused call to index() that crashed v10.0
version 0.2.2   1jun2010  all files in one
version 0.2.1  21may2010  drops use of dicmd
version 0.2    16mar2010  rand() changed to sens(); gphoptions now added `loose'; new options debug robust lpattern() nograph; now calls rctmiss_*.ado not mnar_*.ado; various bug fixes

Test script:
    H:\missing\sensanal\RCTmiss\rctmiss_testscript.do
********************************************************************************/

prog def rctmiss, eclass
version 10

*** SEPARATE PREFIX AND REGRESSION COMMANDS ***
gettoken prefix command : 0, parse(":")
local command : subinstr local command ":" ""
local prefix : subinstr local prefix ":" ""

*** REPLAY ***
if "`command'"=="" {
    if "`e(cmd)'"!="rctmiss" {
        di as error "last estimates not found"
        exit 301
    }
    ereturn display `prefix'
    exit
}

*** PARSE REGRESSION COMMAND ***
gettoken regcmd restofcommand : command
unabcmd `regcmd'
local regcmd = r(cmd)

local 0 `restofcommand'
syntax varlist [if] [in] [fweight aweight iweight pweight], [level(passthru) *]
marksample touse, novarlist
gettoken yvar xvars: varlist
if "`weight'"!="" local weightexp [`weight'`exp']
local regifinwt `if' `in' `weightexp'
local regopts `options'
local level1 `level'

*** PARSE PREFIX COMMAND ***
local 0 `prefix'
syntax, [ ///
    sens(varname) SMDelta(string) PMMDelta(string) nosw AUXiliary(varlist) /// model options
    basemiss(string)                                        /// missing baseline options
    stagger(real -1) COLors(string) LWidth(passthru)        /// graph options
    LPATterns(string) nograph ciband *                      /// graph options
    senstype(string) list LISTOPTions(string) clear         /// output and display options
    savewt(string) eform(string) savedta(string) level(passthru) /// output and display options
    debug mmstore(passthru) SANDwich keepmat(passthru) dicmd neff(string) /// undocumented options
    ]
local gphoptions `options'
local level2 `level'

if !mi("`level1'") & !mi("`level2'") {
    di as error "Please specify level() only once"
    exit 198
}
local 0 , `level1' `level2'
syntax, [level(cilevel)]

if "`regcmd'"=="logistic" {
    if "`eform'"=="" local eform Odds ratio // exponentiate graph
}
if "`eform'"!="" local eformopt eform(`eform')
if "`eform'"!="" local bname "`eform'"
else local bname "Coefficient"

if "`smdelta'"!="" & "`pmmdelta'"!="" {
    di as error "Please specify only one of smdelta() and pmmdelta()"
    exit 198
}
if "`smdelta'"=="" & "`pmmdelta'"=="" {
    if "`sens'"=="" {
        di as error "Assuming pmmdelta(0)"
        local pmmdelta 0
    }
    else {
        di as error "smdelta(numlist) or pmmdelta(numlist) must be specified with sens()"
        exit 198
    }
}

if "`sens'"=="" { // Check for syntax errors
    local 0 ,`gphoptions'
    syntax, [sens *]
    if "`sens'"=="sens" {
        di as error "Syntax for option sens is sens(varname)"
        exit 198
    }
    foreach thing in savedta clear list listoptions gphoptions colors lwidth {
        if "``thing''"!="" & "`badthings'"!="" local s s
        if "``thing''"!="" local badthings `badthings' ``thing''
    }
    if "`badthings'"!="" di as error "sens() not specified, ignoring option`s': `badthings'"
}

local 0 `basemiss'
syntax [anything], [min(int 3)]
local basemissmethod = cond("`anything'"=="", "mean", "`anything'")
if !inlist("`basemissmethod'","mean","mim") {
    di as error "Syntax: basemiss(mean|mim, [min(#)])"
    exit 198
}
local basemissmin `min'

* for output
local col as result _col(26)

if !mi("`neff'") local neffopt neffvalue(`neff')

// IPWs
if "`savewt'"!="" {
    confirm new variable `savewt'
    local savewtopt savewt(`savewt')
}

if "`debug'"=="" local ifdebug *

*** HANDLE INCOMPLETE BASELINES AND AUXILIARIES ***
preserve
tempname orig
if "`xvars'"!="" {
    foreach xvar of varlist `xvars' `auxiliary' {
        qui count if mi(`xvar') & `touse'
        local basemissn = r(N)
        if `basemissn'>0 {
            rename `xvar' `orig'`xvar'
            if "`xvarschanged'"=="" di as text "Incomplete covariate:" _c
            local xvarschanged `xvarschanged' `xvar'
            di `col' "`xvar'" as text " has " as result r(N) as text " missing values"
            di `col' as text "  - imputed with the mean" _c
            if "`basemissmethod'"=="mim" {
                if `basemissn'>=`basemissmin' {
                    di `col' as text " + indicator " as result "M`xvar'"
                    gen M`xvar'=mi(`orig'`xvar') if `touse'
                    local mvars `mvars' M`xvar'
                }
                else {
                    di _new `col' as text "  - no indicator because <`basemissmin' missing values"
                }
            }
            else di
            qui summ `orig'`xvar' if `touse', meanonly
            qui gen `xvar' = cond(mi(`orig'`xvar'), r(mean), `orig'`xvar') if `touse'
        }
    }
}

*** SET UP COMMANDS ***
if !mi("`auxiliary'") local auxopt auxiliary(`auxiliary')
local command `regcmd' `yvar' `xvars' `mvars' `regifinwt', `regopts' `auxopt'
local restofcommand `yvar' `xvars' `mvars' `regifinwt', `regopts' `auxopt'
tempname b V wt neff pstar
if "`smdelta'"!="" {
    * SELECTION MODEL / IPW METHOD
    local maincmd sm_ipw `command' b(`b') v(`V') neff(`neff') pstar(`pstar') `savewtopt' `sw' `mmstore' `neffopt'
    local deltaname1 SM 
}
if "`pmmdelta'"!="" {
    * PATTERN MIXTURE MODEL / MEAN SCORE METHOD
    if "`regcmd'"!="regress" | "`sandwich'"=="sandwich" | !mi("`auxiliary'") {
        local maincmd pmm_glm3 `command' b(`b') v(`V') neff(`neff') pstar(`pstar') `keepmat'
    }
    else {
        local maincmd pmm_reg `restofcommand' b(`b') v(`V') neff(`neff') pstar(`pstar') 
    }
    local deltaname1 PMM 
}
* PARSE DELTA
local 0 `smdelta'`pmmdelta'
syntax anything, [log base(string)]
local delta `anything'
local deltaname2 = cond("`log'"=="log", "exp(delta)", "delta")
if "`base'"=="" local deltabase = 0
else {
    confirm number `base'
    if "`log'"=="log" & "`base'"=="0" local deltabase -999
    else local deltabase = `log'(`base')
}
local deltaname `deltaname1' `deltaname2' 
local deltaparm = lower("`deltaname1'`deltaname2'")
local assumption MNAR
if "`sens'"=="" {
    cap assert `delta'==0
    if !_rc local assumption MAR
}

// COUNT & REPORT OBS
qui count if `touse'
local ntot = r(N)
qui count if `touse' & !mi(`yvar')
local nobs = r(N)
local nmis = `ntot'-`nobs'    
local missfate = cond("`assumption'" == "MNAR", "used in analysis", "ignored")
di as text "Missing data assumption: " `col' "`assumption'"
di as text "Observed outcomes:" `col' `nobs'
di as text "Unobserved outcomes:" `col' `nmis' " (`missfate')"
if !mi("`auxiliary'") di as text "Auxiliary variables:" `col' "`auxiliary'"

*** ANALYSIS ***
if "`sens'"=="" {
    * EXPRESSION SPECIFIED: SINGLE RUN
    tempvar deltavble
    cap gen `deltavble'=`delta' if `touse'
    if _rc {
        di as error "Syntax without sens(): pmmdelta(expression) or smdelta(expression)"
        exit 198
    }
    if "`log'"=="log" {
        qui count if `delta'<0
        if r(N)>0 {
            di as error r(N) " individuals have negative delta - not allowed with log option"
            exit 498
        }
        qui replace `deltavble' = log(`deltavble')
        qui replace `deltavble'=-999 if `delta'==0 
    }
    * catch missing values of delta
    qui count if mi(`deltavble') & `touse'
    if r(N)>0 {
        cap assert mi(`deltavble') if `touse'
        if _rc di as error "`deltaparm' could not be computed for " r(N) " observations"
        else di as error "`deltaparm' could not be computed"
        exit 498
    }
    * display delta
    qui summ `deltavble' if `touse'
    if r(sd)==0 di as text "`deltaname': " `col' r(mean)
    else di as text "`deltaname': " `col' "`log' `delta'"
      
    `ifdebug' di as text `"Running command: `maincmd' delta(`deltavble') `debug'"'
    `dicmd' `maincmd' delta(`deltavble') `debug'
    if "`regcmd'" == "regress" {
        local dof = `neff' - `pstar'
        local dof dof(`dof')
    }
    `ifdebug' di as text "Running command: ereturn post `b' `V', depname(`yvar') obs(`ntot') esample(`touse') `dof'"
    ereturn post `b' `V', depname(`yvar') esample(`touse') `dof'
    foreach stat in neff pstar ntot nobs nmis {
        ereturn scalar `stat' = ``stat''
    }
    ereturn local method = word("`maincmd'",1)
    if "`smdelta'"!="" & "`savewt'"!="" ereturn local IPW `savewt'
    `ifdebug' di as text "*** Final results ***"
    di as text "Effective sample size: " `col' `neff'  
    ereturn display, `eformopt' level(`level')
    ereturn local cmd rctmiss
    foreach xvar in `xvarschanged' {
        drop `xvar'
        rename `orig'`xvar' `xvar'
    }
    if "`mvars'"!="" drop `mvars'
    restore, not
}
else {
    * SENSITIVITY ANALYSIS
    cap numlist "`delta'"
    if _rc {
        di as error "Syntax with sens(): pmmdelta(numlist [,log base(#)]) or smdelta(numlist [,log base(#)])"
        exit 198
    }
    if wordcount(r(numlist))==1 di as error "Warning: only one value in delta: graph will look weird"
    // check some output is requested
    if "`graph'"=="nograph" & "`savedta'"=="" & "`clear'"=="" & "`list'"=="" {
        di as error "Nograph option, please specify one or more of: list, savedta(), clear"
        exit 498
    }
    qui levelsof `sens' if `touse', local(randlevels)
    if wordcount("`randlevels'")>2 {
        di as error "Sorry, I can only handle two-arm trials at present"
        exit 498
    }
    if wordcount("`randlevels'")<2 {
        di as error "`sens' does not vary"
        exit 498
    }
    di as text "Performing sensitivity analyses" _c
    local randcon = word("`randlevels'",1)
    local randint = word("`randlevels'",2)
    local randlab0 : label (`sens') `randcon'
    local randlab1 : label (`sens') `randint'
    tempname post
    if "`savedta'"=="" tempfile savedtafile
    else {
        parse "`savedta'", parse(",")
        local savedtafile `1'
        local savedtareplace `3'
    }
    postfile `post' type delta b se dof neff using `savedtafile', `savedtareplace'
    if inlist("`senstype'","equal","both") local typelist 2
    else if inlist("`senstype'","unequal","one") local typelist 1 3
    else local typelist 1 2 3
    foreach del of numlist `delta' {
        di "." _c
        foreach type in `typelist' {
            local logdel = cond("`log'"=="log" & `del'==0, -999, `log'(`del'))
            if `type'==1 local deltavar cond(`sens'==`randint',`logdel',`deltabase')
            if `type'==2 local deltavar `logdel'
            if `type'==3 local deltavar cond(`sens'==`randcon',`logdel',`deltabase')
            `ifdebug' di as input _new "delta=`logdel', type=`type'"
            `ifdebug' di as input "`maincmd' delta(`deltavar') `debug'"
            `dicmd' qui `maincmd' delta(`deltavar') `debug'
            if "`regcmd'" == "regress" local dof = `neff' - `pstar'
            else local dof .
            ereturn post `b' `V'
            post `post' (`type') (`logdel') (_b[`sens']) (_se[`sens']) (`dof') (scalar(`neff'))
        }
    }
    di
    postclose `post'

    use `savedtafile', clear
    label def type 1 "`randlab1' only" 2 "both arms" 3 "`randlab0' only"
    label val type type
    * sort out x-variable
    if "`log'"=="log" {
        * want delta output and graphed on exp-scale
        gen exp_delta = exp(delta)
        gen deltagraph = exp(delta)
        label var exp_delta "exp(delta)"
        label var deltagraph "exp(delta), staggered for graph"
        local dlistvar exp_delta
    }
    else {
        gen deltagraph = delta
        label var deltagraph "delta, staggered for graph"
        local dlistvar delta
    }
    * sort out y-variable
    gen zcrit = cond(`dof' == ., invnorm(.5+`level'/200), invttail(dof,.5-`level'/200))
    if "`eform'"!="" {
        gen exp_b = exp(b)
        gen exp_b_low = exp(b-zcrit*se)
        gen exp_b_upp = exp(b+zcrit*se)
        local gphoptions yscale(log) `gphoptions'
        local blistvars exp_b exp_b_low exp_b_upp
        local bvar exp_b
    }
    else {
        gen b_low = b-zcrit*se
        gen b_upp = b+zcrit*se
        local blistvars b se
        local bvar b
    }
    
    if "`list'"=="list" {    
        if "`listoptions'"=="" local listoptions sepby(delta) abb(10)
        list type `dlistvar' `blistvars', `listoptions'
    }
    
    if "`graph'"!="nograph" {
        *** DRAW A GRAPH
        di "Drawing graph..."
        local col1 = word("`colors'",1)
        local col2 = word("`colors'",2)
        local col3 = word("`colors'",3)
        if "`col1'"=="" local col1 blue
        if "`col2'"=="" local col2 purple
        if "`col3'"=="" local col3 red
        local lpattern1 = word("`lpatterns'",1)
        local lpattern2 = word("`lpatterns'",2)
        local lpattern3 = word("`lpatterns'",3)
        if "`ciband'"=="" { // confidence limits as rspikes
            if `stagger'<0 {
                qui sum deltagraph, meanonly
                local stagger = (r(max)-r(min))/100
            }
            qui replace deltagraph=deltagraph-`stagger' if type==1
            qui replace deltagraph=deltagraph+`stagger' if type==3
            if "`lpattern1'"!="" local lpattern1 lpattern(`lpattern1')
            if "`lpattern2'"!="" local lpattern2 lpattern(`lpattern2')
            if "`lpattern3'"!="" local lpattern3 lpattern(`lpattern3')
            local legendboth label(3 "both arms") 
            local legendone label(1 "`randlab1' only") label(5 "`randlab0' only")
            if "`senstype'"=="both" local legendopt legend(order(3) `legendboth' rows(1))
            else if "`senstype'"=="one" local legendopt legend(order(1 5) `legendone' rows(1))
            else local legendopt legend(order(1 3 5) `legendboth' `legendone' rows(1))
            #delimit ;
            local graphcmd twoway
                (line `bvar' deltagraph if type==1, lcol(`col1') `lwidth' `lpattern1') 
                (rspike `bvar'_low `bvar'_upp deltagraph if type==1, lcol(`col1') `lwidth' `lpattern1')
                (line `bvar' deltagraph if type==2, lcol(`col2') `lwidth' `lpattern2') 
                (rspike `bvar'_low `bvar'_upp deltagraph if type==2, lcol(`col2') `lwidth' `lpattern2')
                (line `bvar' deltagraph if type==3, lcol(`col3') `lwidth' `lpattern3') 
                (rspike `bvar'_low `bvar'_upp deltagraph if type==3, lcol(`col3') `lwidth' `lpattern3')
                ,
                `legendopt'
                ytitle("`bname' for `sens' (`level'% CI)")
                xtitle(`deltaname' in specified arm(s))
                `gphoptions'
            ;
            #delimit cr
        }
        else { // confidence limits as lines
            if "`lpattern1'"=="" local lpattern1 solid
            if "`lpattern2'"=="" local lpattern2 dash
            local lpattern lpattern(`lpattern1' `lpattern2' `lpattern2')
            local legendboth label(4 "both arms") 
            local legendone label(1 "`randlab1' only") label(7 "`randlab0' only")
            if "`senstype'"=="both" local legendopt legend(order(4) `legendboth' rows(1))
            else if "`senstype'"=="one" local legendopt legend(order(1 7) `legendone' rows(1))
            else local legendopt legend(order(1 4 7) `legendboth' `legendone' rows(1))
            #delimit ;
            local graphcmd twoway
                (line `bvar' `bvar'_low `bvar'_upp deltagraph if type==1, lcol(`col1' `col1' `col1') `lwidth' `lpattern') 
                (line `bvar' `bvar'_low `bvar'_upp deltagraph if type==2, lcol(`col2' `col2' `col2') `lwidth' `lpattern') 
                (line `bvar' `bvar'_low `bvar'_upp deltagraph if type==3, lcol(`col3' `col3' `col3') `lwidth' `lpattern') 
                ,
                `legendopt'
                ytitle(`bname' for `sens')
                xtitle(`deltaname' in specified arm(s))
                `gphoptions'
            ;
            #delimit cr
        }
        `ifdebug' di as text `"*** Running: `graphcmd'"'
        `graphcmd'
        if "`clear'"!="" {
             global F9 `graphcmd'
             di as text "Graph command stored in F9"
        }
    }
    if "`clear'"!="" {
         restore, not
    }
    ereturn clear // Nothing sensible to ereturn
}
end

*************************** END OF RCTMISS PROGRAM *******************************************

prog def pmm_reg
version 10
syntax varlist(min=1) [if] [in], [delta(string) robust debug noCONStant b(string) V(string) neff(string) pstar(string)]

// PARSE
marksample touse, novarlist
gettoken y xlist : varlist
if "`debug'"=="" local ifdebug qui
di as text "Method:" _col(26) as result "two linear regressions"

tempname bI vI vIlarge bD vD vDlarge vlarge

// IMPUTATION MODEL
`ifdebug' di as text "*** Imputation model ***"
`ifdebug' reg `y' `xlist' if `touse', `robust' `constant'
local factorI = e(N) / e(df_r)
mat `bI' = e(b)
mat `vI' = e(V)
mat `vIlarge' = e(V) * e(df_r) / e(N) // Correction factor is (nobs-pstar)/nobs

// CORRECTION MODEL 
// fitted to all obs
tempvar mdz
qui gen `mdz' =  mi(`y') * `delta' if `touse'
`ifdebug' di as text "*** Analysis model ***"
`ifdebug' reg `mdz' `xlist' if `touse', robust `constant'
mat `bD' = e(b)
mat `vD' = e(V)
mat `vDlarge' = e(V) * e(df_r) / e(N) // Correction factor is (ntot-pstar)/ntot

// COMBINED
mat `b' = `bI' + `bD'
mat `v' = `vI' + `vD'
mat `vlarge' = `vIlarge' + `vDlarge'
scalar `pstar' = colsof(`b')
local dfcorr = (det(`vlarge')/det(`v'))^(1/`pstar') // Estimates (neff-pstar)/neff
scalar `neff' = `pstar'  / (1 - `dfcorr' )

* sureg fails because it requires the same obs for both regns (and the same weights)
* but I verified that the residuals are exactly uncorrelated

end

******************************** END OF PMM_REG PROGRAM ************************************

prog def pmm_glm
version 10

/*******************************************
NOTES
    Avoid -predict, residual- which uses unexpected formulae.
*******************************************/

// PARSE
syntax anything [if] [in], delta(string) [debug noCONStant b(string) V(string) neff(string) pstar(string) keepmat(string)]
gettoken cmd varlist : anything
unabcmd `cmd'
local cmd = r(cmd)
gettoken y xlist : varlist
if "`debug'"=="" local ifdebug qui
di as text "Method:" _col(26) as result "mean score + joint sandwich variance"

// SET UP
marksample touse, novarlist
tempvar rowmiss id ri pi ra pa ystar offsetvar rai r2
tempname bi Vi Va Vpred Vdrop Vmiss
tempname fullV Binv B BIA BII BAI BAA C CAI CAA CII
gen `id'=_n

// COUNT OBS
qui count if `touse'
local ntot = r(N)
qui count if `touse' & !mi(`y')
local nobs = r(N)

// IMPUTATION MODEL
`ifdebug' di as text _new "*** Imputation model ***"
if "`cmd'" != "regress" {
    qui gen `offsetvar' = missing(`y')*`delta' if `touse'
    local offsetopt offset(`offsetvar')
}
`ifdebug' `cmd' `y' `xlist' if `touse', `offsetopt' `constant'
qui predict `pi' if `touse'
if "`cmd'" == "regress" {
    qui replace `pi' = `pi' + missing(`y')*`delta' if `touse'
}
qui gen `ri' = cond(mi(`y'), 0, `y'-`pi') if `touse'
mat `bi' = e(b)
local p = colsof(`bi')

// ANALYSIS MODEL
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
* CC analysis (only for calculating effective sample size)
`ifdebug' di as text _n "*** CC analysis (only for calculating effective sample size) ***"
`ifdebug' `cmd2' `y' `xlist' if `touse', `opts' `constant' robust
scalar `pstar' = cond("`cmd'"=="regress",`p',1)
mat `Vdrop' = e(V)*(e(N)-`pstar')/e(N)
* main analysis
qui gen `ystar' = cond(missing(`y'),`pi',`y') if `touse'
`ifdebug' di as text _new "*** Analysis model ***"
`ifdebug' `cmd2' `ystar' `xlist' if `touse', `opts' `constant' robust
mat `Vpred' = e(V)*(e(N)-`pstar')/e(N)
qui predict `pa' if `touse'
qui gen `ra' = (`ystar' - `pa') if `touse'
mat `b' = e(b)

// CONSTRUCT C MATRIX
sort `id'
qui gen `rai'=sqrt(`ri'*`ra') if `touse'
mat opaccum `CII' = `xlist' if `touse', group(`id') opvar(`ri') `constant'
mat opaccum `CAA' = `xlist' if `touse', group(`id') opvar(`ra') `constant'
mat opaccum `CAI' = `xlist' if `touse', group(`id') opvar(`rai') `constant'
mat `C' = (`CAA', `CAI' \ `CAI', `CII')

// CONSTRUCT B MATRIX 
tempvar hprimeA hprimeI opAA opAI opII
if "`cmd'"=="logit" | "`cmd'"=="logistic" {
    qui gen `hprimeA' = `pa'*(1-`pa')
    qui gen `hprimeI' = `pi'*(1-`pi')
}
else if "`cmd'"=="regress" {  
    qui gen `hprimeA' = 1
    qui gen `hprimeI' = 1
}
else if "`cmd'"=="poisson" {
    qui gen `hprimeA' = `pa'
    qui gen `hprimeI' = `pi'
}
gen `opAA' = sqrt(`hprimeA')
gen `opAI' = sqrt(`hprimeI') * mi(`y')
gen `opII' = sqrt(`hprimeI') * !mi(`y')
mat opaccum `BAA' = `xlist' if `touse', group(`id') opvar(`opAA') `constant'
mat opaccum `BAI' = `xlist' if `touse', group(`id') opvar(`opAI') `constant'
mat `BAI' = -`BAI'
mat opaccum `BII' = `xlist' if `touse', group(`id') opvar(`opII') `constant'
mat `BIA' = J(`p',`p',0)
mat `B' = (`BAA', `BAI' \ `BIA', `BII') 

if "`debug'"=="debug" {
    mat l `B', title("B")
    mat l `C', title("C")
}

// CALCULATE V MATRIX
mat `Binv' = inv(`B')
mat `fullV' = `Binv' * `C' * `Binv'' 
mat `v' = `fullV'[1..`p',1..`p']
if "`debug'"=="debug" {
    mat l `Binv', title("Binv")
    mat l `fullV', title("fullV")
}
mat `Vmiss' = `v'

// EFFECTIVE SAMPLE SIZE
local Ipred = det(`Vpred')^(-1/`p')
local Idrop = det(`Vdrop')^(-1/`p')
local Imiss = det(`Vmiss')^(-1/`p')
local misswt = (`Imiss'-`Idrop') / (`Ipred'-`Idrop')
if `misswt'<0 local misswt=0 
if `misswt'>1 local misswt=1
if "`debug'"=="debug" di as text "det(V)^(-1/`p') = " as result `Imiss' ///
    as text ", det(Vdrop)^(-1/`p') = " as result `Idrop' ///
    as text ", det(Vpred)^(-1/`p') = " as result `Ipred' ///
    as text ", misswt = " as result `misswt'
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

*************************** END OF PMM_GLM PROGRAM ****************************************

prog def sm_ipw
version 10
syntax anything [if] [in], [delta(string) robust debug noSUMwt savewt(string) ///
	AUXiliary(varlist) noCONStant nosw ///
	b(string) V(string) neff(string) pstar(string) mmstore(string)]

// PARSE
marksample touse, novarlist
gettoken cmd varlist : anything
unabcmd `cmd'
local cmd = r(cmd)
gettoken y xlist : varlist
if "`debug'"=="" local ifdebug qui
di as text "Method:" _col(26) as result "inverse probability weighting"

qui count if `touse'
local ntot = r(N)
qui count if `touse' & !mi(`y')
local nobs = r(N)
local nmis = `ntot'-`nobs'

if `nmis'==0 {
    di as error "No incomplete observations: no weights used"
    `cmd' `varlist', robust `constant'
    exit
}
local col as result _col(26)

// FIT MISSINGNESS MODEL
`ifdebug' di _new as text "*** Fitting missingness model ***"
tempvar miss offset lp1 lp2 weight
qui gen `miss' = mi(`y') if `touse'
qui gen `offset' = cond(`miss',0,`delta'*`y')
qui ml model lf rctmiss_smlik (`miss' = `xlist' `auxiliary', offset(`offset') `constant'), vce(robust)
`ifdebug' ml maximize
if "`mmstore'"!="" est store `mmstore'
qui predict `lp1'

// COMPUTE WEIGHTS
qui gen `weight' = 1 + exp(`lp1') if `touse'
if "`sw'"!="nosw" {
    // FIT MAR MISSINGNESS MODEL
    `ifdebug' di _new as text "*** Fitting MAR missingness model ***"
	* without auxiliary!
    qui ml model lf rctmiss_smlik (`miss' = `xlist', `constant'), vce(robust)
    `ifdebug' ml maximize
    qui predict `lp2'
    qui replace `weight' = `weight' / (1 + exp(`lp2')) if `touse'
    di as text "Weights: " `col' "stabilised"
}
else di as text "Weights: " `col' "not stabilised"

// SUMMARISE WEIGHTS
if "`sumwt'" != "nosumwt" {
	local col2 _col(45)
    qui summ `weight' if `touse' & !mi(`y')
    di as text "Summary of weights:" `col' as text "CV  = " as result r(sd)/r(mean) `col2' as text "Max/min = " as result r(max)
    di `col' as text "Max = " as result r(max)/r(min) `col2' as text "Min     = " as result r(min) 
    local wts = r(N)    
    qui count if `weight'==0 & `touse' & !mi(`y')
    local wt0 = r(N)
    di `col' as text ">0  = " as result `wts'-`wt0' `col2' as text "Zero    = " as result `wt0'
}

// FIT WEIGHTED ANALYSIS MODEL
`ifdebug' di _new as text "*** Fitting weighted analysis model ***" 
`ifdebug' `cmd' `varlist' if `touse' [pw=`weight'], `constant'
mat `b'=e(b)
mat `v'=e(V)

// OPTIONALLY SAVE WEIGHT
if "`savewt'"!="" rename `weight' `savewt'

// Compute neff, pstar???

scalar `neff' = e(N)
scalar `pstar' = colsof(`b')
`ifdebug' di as text "SM_IPW completed successfully"
`ifdebug' scalar dir
end

**************************** END OF SM_IPW PROGRAM ***************************************

prog def dicmd
noi di as input `"`0'"'
`0'
end

**************************** START OF PMM_GLM3 PROGRAM ***************************************

prog def pmm_glm3
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
		tidy up the code - it currently leaves matrices around
15dec2016 v3 - using mata - works nicely
	grafting back into rctmiss	
*******************************************/

// PARSE
syntax anything [if] [in] [iweight/], ///
	delta(string) [AUXiliary(varlist) /// model specification
	b(string) V(string) neff(string) pstar(string) /// returned values
	noCONStant /// optional model specification
	NEFFValue(real 0) /// optional analysis specification
	keepmat(string) /// optional returned values
	debug display INFluence(string)] // output settings

gettoken cmd vars : anything
unabcmd `cmd'
local cmd = r(cmd)
if "`cmd'"=="logistic" local cmd logit
if !inlist("`cmd'","regress","logit","poisson") {
    di as error "Sorry, command `cmd' is not yet supported"
    exit 498
}
gettoken y xlist : vars
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
unab xPlist : `xlist' `auxiliary'
local nxPlist : word count `xPlist'
unab xSlist : `xlist'
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
if "`debug'"=="debug" {
	mac list _pS
	mac list _pP
	foreach thing in BSS BSP BPS BPP b bCC {
		mat list ``thing'', title(`thing')
	}
}
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
	qui gen `influence'obs = .
	qui gen `influence'full = .

	mata: residS = st_data(., "`residS'")
	mata: xS = st_data(., "`xSlist' `one'")
	mata: residP = st_data(., "`residP'")
	mata: xP = st_data(., "`xPlist' `one'")

	mata: Binv=st_matrix("`Binv'")
	mata: v=st_matrix("`v'")
	mata: BSS=st_matrix("`BSS'")

	mata: U = (residS:*xS, residP:*xP) // loose Hadamard product
	mata: A = (I(`pS'), J(`pS',`pP',0))
	mata: ABinvU = A*Binv*U'
	mata: vinv = invsym(v)
	*mata: inf = diagonal(ABinvU'*vinv*ABinvU) // slow
	mata: inf = rowsum((ABinvU'*vinv):*ABinvU') // much faster
	mata: st_store(., "`influence'obs", inf)

	mata: BSSinv=luinv(BSS)
	*mata: inf2 = diagonal(xS*BSSinv'*vinv*BSSinv*xS') // slow 
	mata: xSBSSinv=xS*BSSinv'
	mata: inf2 = rowsum((xSBSSinv*vinv):*xSBSSinv) // much faster 
	mata: st_store(.,"`influence'full",inf2)

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

**************************** END OF PMM_GLM3 PROGRAM ***************************************

