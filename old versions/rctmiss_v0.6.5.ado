*! version 0.6.5   Ian White    8may2012
/*******************************************************************************
HISTORY
version 0.6.5    8may2012
    New sandwich option forces use of sandwich variance
version 0.6.4   30jan2012
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

STILL TO DO
- fix df too to match MAR? see ereturn post,... dof()
- SM: noconstant drops the constant from the missingness model too - is this what we want?
- Alho option? e.g. attempts(varname)
- Correct df for CIs
- return error if user specifies sens(14/14)! - or more helpful graph
- optionally give CIs as dotted lines or spikes
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
    sens(varname) SMDelta(string) PMMDelta(string) nosw /// model options
    basemiss(string) nmissmin(int 3) /// missing baseline options
    stagger(real -1) COLors(string) LWidth(passthru) LPATterns(string) nograph ciband * /// graph options
    senstype(string) list LISTOPTions(string) clear savewt(string) eform(string) savedta(string) level(passthru) /// output and display options
    debug mmstore(passthru) sandwich keepmat(passthru) /// undocumented options
    ]
local gphoptions `options'
local level2 `level'

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

if "`basemiss'"=="" local basemiss mean
else if !inlist("`basemiss'","mean","mim") {
    di as error "Syntax: basemiss(mean|mim)"
    exit 198
}

// IPWs
if "`savewt'"!="" {
    confirm new variable `savewt'
    local savewtopt savewt(`savewt')
}

if "`debug'"=="" local ifdebug *

*** HANDLE INCOMPLETE BASELINES ***
preserve
tempname orig
if "`xvars'"!="" {
    foreach xvar of varlist `xvars' {
        qui count if mi(`xvar') & `touse'
        local nmiss = r(N)
        if `nmiss'>0 {
            if inlist("`basemiss'","mean","mim") {
                rename `xvar' `orig'`xvar'
                local xvarschanged `xvarschanged' `xvar'
                di as text "Imputing " as result r(N) as text " missing values of " as result "`xvar'" as text " with the observed mean"
                if "`basemiss'"=="mim" {
                    if `nmiss'>=`nmissmin' {
                        di as text "  - and including missing indicator " as result "M`xvar'"
                        gen M`xvar'=mi(`orig'`xvar') if `touse'
                        local mvars `mvars' M`xvar'
                    }
                    else {
                        di as text "  - not including missing indicator because <`nmissmin' missing values (nmissmin option)"
                    }
                }
                qui summ `orig'`xvar' if `touse', meanonly
                qui gen `xvar' = cond(mi(`orig'`xvar'), r(mean), `orig'`xvar') if `touse'
            }
        }
    }
}

*** SET UP COMMANDS ***
local command `regcmd' `yvar' `xvars' `mvars' `regifinwt', `regopts'
local restofcommand `yvar' `xvars' `mvars' `regifinwt', `regopts'
tempname b V wt
if "`smdelta'"!="" {
    * SELECTION MODEL / IPW METHOD
    local maincmd sm_ipw `command' b(`b') v(`V') `savewtopt' `sw' `mmstore'
    local deltaname1 SM 
}
if "`pmmdelta'"!="" {
    * PATTERN MIXTURE MODEL / MEAN SCORE METHOD
    if "`regcmd'"!="regress" | "`sandwich'"=="sandwich" {
        local maincmd pmm_glm `command' b(`b') v(`V') `keepmat'
    }
    else {
        local maincmd pmm_reg `restofcommand' b(`b') v(`V')
    }
    local deltaname1 PMM 
}
* PARSE DELTA
local 0 `smdelta'`pmmdelta'
syntax anything, [log base(real 0)]
local delta `anything'
local deltaname2 = cond("`log'"=="log", "exp(delta)", "delta")
local deltabase = `log'(`base')
if "`log'"=="log" & `base'==0 local deltabase -999
local deltaname `deltaname1' `deltaname2' 
local deltaparm = lower("`deltaname1'`deltaname2'")
local mnar 1
if "`sens'"=="" {
    cap assert `delta'==0
    if !_rc local mnar 0
}

// COUNT & REPORT OBS
qui count if `touse'
local n = r(N)
qui count if `touse' & !mi(`yvar')
local ncomplete = r(N)
local nincomplete = `n'-`ncomplete'    
if `mnar' {
    di as text "Using " as result `ncomplete' as text " observed outcomes and " as result `nincomplete' as text " unobserved outcomes"
    di as text "Results allowing for MNAR"
}
else {
    di as text "Using " as result `ncomplete' as text " observed outcomes but ignoring " as result `nincomplete' as text " unobserved outcomes"
    di as text "Results assuming MAR"
}

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
    if r(sd)==0 di as text "`deltaname': " as result r(mean)
    else di as text "`deltaname': " as result "`log' `delta'"
    * options
    if "`gphoptions'`colors'`lwidth'"!="" di as error "Options `gphoptions' `colors' `lwidth' ignored"
      
    `ifdebug' di as text `"Running command: `maincmd' delta(`deltavble') `debug'"'
    `maincmd' delta(`deltavble') `debug'
    ereturn post `b' `V', depname(`yvar') obs(`n') esample(`touse')
    if "`smdelta'"!="" & "`savewt'"!="" ereturn local IPW `savewt'
    `ifdebug' di as text "*** Final results ***"
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
    // check some output is requested
    if "`graph'"=="nograph" & "`savedta'"=="" & "`clear'"=="" & "`list'"=="" {
        di as error "Please specify either list or savedta() or clear with nograph"
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
    local randlab0 : label (`sens') 0
    local randlab1 : label (`sens') 1
    tempname post
    if "`savedta'"=="" tempfile savedtafile
    else {
        parse "`savedta'", parse(",")
        local savedtafile `1'
        local savedtareplace `3'
    }
    postfile `post' type delta b se using `savedtafile', `savedtareplace'
    if "`senstype'"=="both" local typelist 2
    else if "`senstype'"=="one" local typelist 1 3
    else local typelist 1 2 3
    foreach del of numlist `delta' {
        di "." _c
        foreach type in `typelist' {
            local logdel = cond("`log'"=="log" & `del'==0, -999, `log'(`del'))
            if `type'==1 local deltavar cond(`sens'==`randint',`logdel',`deltabase')
            if `type'==2 local deltavar `logdel'
            if `type'==3 local deltavar cond(`sens'==`randcon',`logdel',`deltabase')
            qui `maincmd' delta(`deltavar') `debug'
            ereturn post `b' `V'
            post `post' (`type') (`logdel') (_b[`sens']) (_se[`sens'])
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
    local zcrit = invnorm(.5+`level'/200)
    if "`eform'"!="" {
        gen exp_b = exp(b)
        gen exp_b_low = exp(b-`zcrit'*se)
        gen exp_b_upp = exp(b+`zcrit'*se)
        local gphoptions yscale(log) `gphoptions'
        local blistvars exp_b exp_b_low exp_b_upp
        local bvar exp_b
    }
    else {
        gen b_low = b-`zcrit'*se
        gen b_upp = b+`zcrit'*se
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
                (rspike `bvar'_low `bvar'_upp deltagraph if type==1, lcol(`col1') `lwidth')
                (line `bvar' deltagraph if type==2, lcol(`col2') `lwidth' `lpattern2') 
                (rspike `bvar'_low `bvar'_upp deltagraph if type==2, lcol(`col2') `lwidth')
                (line `bvar' deltagraph if type==3, lcol(`col3') `lwidth' `lpattern3') 
                (rspike `bvar'_low `bvar'_upp deltagraph if type==3, lcol(`col3') `lwidth')
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
        `ifdebug' "*** Running: `graphcmd'"
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

*************************** END OF MAIN PROGRAM *******************************************

prog def pmm_reg
version 10
syntax varlist(min=1) [if] [in], [delta(string) robust debug noCONStant b(string) V(string)]

// PARSE
marksample touse, novarlist
gettoken y xlist : varlist
if "`debug'"=="" local ifdebug qui
di as text "Method: two linear regressions"

// IMPUTATION MODEL
`ifdebug' di as text "*** Imputation model ***"
`ifdebug' reg `y' `xlist' if `touse', `robust' `constant'
global df_m = e(df_m)
local factorI = e(N) / e(df_r)
mat `b'=e(b)
mat `v'=e(V)

// ANALYSIS MODEL
tempvar mdz
if "`delta'"!="" {
    qui gen `mdz' =  mi(`y') * `delta' if `touse'
    `ifdebug' di as text "*** Analysis model ***"
    `ifdebug' reg `mdz' `xlist' if `touse', robust `constant'
    local factorA = e(N) / e(df_r)
    mat `b' = `b' + e(b)
    mat `v' = `v' + e(V)*`factorI'/`factorA'
}

* sureg fails because it requires the same obs for both regns (and the same weights)
* but I verified that the residuals are exactly uncorrelated

end

******************************** END OF PMM_REG **************************************

prog def pmm_glm
version 10

/*******************************************
NOTES
    Avoid -predict, residual- which uses unexpected formulae.
*******************************************/

// PARSE
syntax anything [if] [in], delta(string) [debug noCONStant b(string) V(string) keepmat(string)]
gettoken cmd varlist : anything
unabcmd `cmd'
local cmd = r(cmd)
gettoken y xlist : varlist
if "`debug'"=="" local ifdebug qui
di as text "Method: mean score + joint sandwich variance"

// SET UP
marksample touse, novarlist
tempvar rowmiss id ri pi ra pa ystar offsetvar rai r2
tempname bi Vi Va Vpred Vdrop Vmiss
tempname fullV Binv B BIA BII BAI BAA C CAI CAA CII
gen `id'=_n

// COUNT OBS
qui count if `touse'
local n = r(N)
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
`ifdebug' `cmd2' `y' `xlist' if `touse', `opts' `constant' robust
local dfcorrection = cond("`cmd'"=="regress",`p',1)
mat `Vdrop' = e(V)*(e(N)-`dfcorrection')/e(N)
* main analysis
qui gen `ystar' = cond(missing(`y'),`pi',`y') if `touse'
`ifdebug' di as text _new "*** Analysis model ***"
`ifdebug' `cmd2' `ystar' `xlist' if `touse', `opts' `constant' robust
mat `Vpred' = e(V)*(e(N)-`dfcorrection')/e(N)
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
local Vpred = det(`Vpred')^(-1/`p')
local Vdrop = det(`Vdrop')^(-1/`p')
local Vmiss = det(`Vmiss')^(-1/`p')
local misswt = (`Vmiss'-`Vdrop') / (`Vpred'-`Vdrop')
local neff = `nobs' + `misswt'*(`n'-`nobs')
if "`debug'"=="debug" di "Vpred = `Vpred', Vmiss = `Vmiss', Vdrop = `Vdrop'
di as text "Effective sample size: " as result `neff'

// DF CORRECTION
local factor = `neff' / (`neff' - `dfcorrection')
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

*************************** END OF PMM_GLM *******************************************

prog def sm_ipw
version 10
syntax anything [if] [in], [delta(string) robust debug noSUMwt savewt(string) noCONStant b(string) V(string) nosw mmstore(string)]

// PARSE
marksample touse, novarlist
gettoken cmd varlist : anything
unabcmd `cmd'
local cmd = r(cmd)
gettoken y xlist : varlist
if "`debug'"=="" local ifdebug qui
di as text "Method: inverse probability weighting"

qui count if `touse'
local n = r(N)
qui count if `touse' & !mi(`y')
local ncomplete = r(N)
local nincomplete = `n'-`ncomplete'

if `nincomplete'==0 {
    di as error "No incomplete observations: no weights used"
    `cmd' `varlist', robust `constant'
    exit
}

// FIT MISSINGNESS MODEL
`ifdebug' di _new as text "*** Fitting missingness model ***"
tempvar miss offset lp1 lp2 weight
qui gen `miss' = mi(`y') if `touse'
qui gen `offset' = cond(`miss',0,`delta'*`y')
qui ml model lf rctmiss_smlik (`miss' = `xlist', offset(`offset') `constant'), vce(robust)
`ifdebug' ml maximize
if "`mmstore'"!="" est store `mmstore'
qui predict `lp1'

// COMPUTE WEIGHTS
qui gen `weight' = 1 + exp(`lp1') if `touse'
if "`sw'"!="nosw" {
    // FIT MAR MISSINGNESS MODEL
    `ifdebug' di _new as text "*** Fitting MAR missingness model ***"
    qui ml model lf rctmiss_smlik (`miss' = `xlist', `constant'), vce(robust)
    `ifdebug' ml maximize
    qui predict `lp2'
    qui replace `weight' = `weight' / (1 + exp(`lp2')) if `touse'
    di as text "Weights: " as result "stabilised"
}
else di as text "Weights: " as result "not stabilised"

// SUMMARISE WEIGHTS
if "`sumwt'" != "nosumwt" {
    qui summ `weight' if `touse' & !mi(`y')
    di _new as text "Summary of weights: CV  = " as result r(sd)/r(mean) _col(40) as text " Max/min = " as result r(max)
    di      as text "                    Max = " as result r(max)/r(min) _col(40) as text " Min     = " as result r(min) 
    local wts = r(N)    
    qui count if `weight'==0 & `touse' & !mi(`y')
    local wt0 = r(N)
    di      as text "                    >0  = " as result `wts'-`wt0' _col(40) as text " Zero    = " as result `wt0'
}

// FIT WEIGHTED ANALYSIS MODEL
`ifdebug' di _new as text "*** Fitting weighted analysis model ***" 
qui `cmd' `varlist' if `touse' [pw=`weight'], `constant'
mat `b'=e(b)
mat `v'=e(V)

// OPTIONALLY SAVE WEIGHT
if "`savewt'"!="" rename `weight' `savewt'
end

**************************** END OF SM_IPW ******************************************

