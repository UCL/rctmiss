*! version 0.3 28-30jun2010  new basemiss and nmissmin options; noconstant works properly; doesn't replay after sens() option; new options cistyle(line) savewt() senstype(); log tidied up
* version 0.2.4 24jun2010  eform works properly; allow log as first element of smdelta or pmmdelta (is everything labeled correctly? no, output data & list aren't) ***
* version 0.2.3 14jun2010  deleted unused call to index() that crashed v10.0
* version 0.2.2   1jun2010  all files in one
* version 0.2.1  21may2010  drops use of dicmd
* version 0.2    16mar2010  rand() changed to sens(); gphoptions now added `loose'; new options debug robust lpattern() nograph; now calls rctmiss_*.ado not mnar_*.ado; various bug fixes
prog def rctmiss, eclass
version 10
/*************************************************************************
rctmiss: analyse a RCT allowing for informatively missing outcome data.

Notes:
    missing values of baseline covariates are imputed with the mean of their observed values (White & Thompson, 2005).
    user can choose a fixed analysis or a sensitivity analysis
    user can express departures from MAR in two ways
    - as the coefficient of y in a model for m on x and y (selection model): the coefficient is the increase in the log odds of missingness for a 1-unit increase in y (the log imor)
    - as the coefficient of m in a model for y on x and m (pattern mixture model): for linear regression, the coefficient is the difference between mean unobserved outcome and the mean observed outcome; for logistic regression, the coefficient is the log odds ratio between the outcome and missingness (the log imor) 

    After smdelta, the same results can be produced by
        reg varlist [pw=`e(IPW)']

SYNTAX
    rctmiss, [sens(varname)] [smdelta(IM_expression)|pmmdelta(IM_expression)] [more_options] [graph_options]: regression_command
    
where more_options can be
    stagger(real -1) colors(string) saving(string) clear replace list LISTOPTions(string) nopreserve

IM_expression can be
    expression -> gives a single analysis
    [log] numlist -> gives a sensitivity analysis; sens(varname) is required

STILL TO DO
    SM: noconstant drops the constant from the missingness model too - is this what we want?
    Alho option e.g. attempts(varname)
    Correct df for CIs
    help file
        return error if user specifies sens(14/14)! - or more helpful graph
        optionally give CIs as dotted lines or spikes

LIMITATIONS
    only 2 arms
    not (st)cox
*************************************************************************/

*** SEPARATE PREFIX AND REGRESSION COMMANDS ***
local strpos = strpos("`0'",":")
if `strpos'==0 {
    local prefix `0'
}
else {
    local prefix = substr("`0'",1,`strpos'-1)
    local command = substr("`0'",`strpos'+1,.)
}

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
syntax varlist [if] [in] [fweight aweight iweight pweight], [*]
marksample touse, novarlist
gettoken yvar xvars: varlist
if "`weight'"!="" local weightexp [`weight'`exp']
local regifinwt `if' `in' `weightexp'
local regopts `options'

*** PARSE PREFIX COMMAND ***
local 0 `prefix'
syntax, [ ///
    sens(varname) SMDelta(string) PMMDelta(string) /// model options
    basemiss(string) nmissmin(int 3) /// missing baseline options
    stagger(real -1) COLors(string) LWidth(passthru) LPATterns(string) nograph CIstyle(string) * /// graph options
    senstype(string) list LISTOPTions(string) clear savewt(string) eform(passthru) saving(string) replace /// output and display options
    debug /// undocumented options
    ]
local gphoptions `options'

if "`smdelta'"!="" & "`pmmdelta'"!="" {
    di as error "Please specify only one of smdelta() and pmmdelta()"
    exit 498
}
if "`smdelta'"=="" & "`pmmdelta'"=="" {
    if "`sens'"=="" {
        di as error "Assuming smdelta(0)"
        local smdelta 0
    }
    else {
        di as error "smdelta(numlist) or pmmdelta(numlist) must be specified with sens()"
        exit 498
    }
}

if "`basemiss'"!="mim" local basemiss mean

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
    local maincmd sm_ipw `command' b(`b') v(`V') `savewtopt'
    local deltaname1 SM 
}
if "`pmmdelta'"!="" {
    * PATTERN MIXTURE MODEL / MEAN SCORE METHOD
    if "`regcmd'"!="regress" {
        local maincmd pmm_glm `command' b(`b') v(`V')
    }
    else {
        local maincmd pmm_reg `restofcommand' b(`b') v(`V')
    }
    local deltaname1 PMM 
}
local delta `smdelta'`pmmdelta'
if word("`delta'",1)=="log" {
    gettoken log delta : delta
    local deltaname2 exp(delta)
}
else local deltaname2 delta
local deltaname `deltaname1' `deltaname2' 
local deltaparm = lower("`deltaname1'`deltaname2'")

// COUNT & REPORT OBS
qui count if `touse'
local n = r(N)
qui count if `touse' & !mi(`yvar')
local ncomplete = r(N)
local nincomplete = `n'-`ncomplete'    
if "`delta'"!="" {
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
    qui gen `deltavble'=`delta' if `touse'
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
    if r(sd)==0 di as text "`deltaname' = " as result r(mean)
    else di as text "`deltaname' = " as result "`delta'"
    * options
    if "`gphoptions'`colors'`lwidth'"!="" di as error "Options `gphoptions' `colors' `lwidth' ignored"
    `ifdebug' di as text `"Running command: `maincmd' delta(`log'(`delta')) `debug'"'
    `maincmd' delta(`log'(`delta')) `debug'
    ereturn post `b' `V', depname(`yvar') obs(`n') esample(`touse')
    if "`smdelta'"!="" & "`savewt'"!="" ereturn local IPW `savewt'
    `ifdebug' di as text "*** Final results ***"
    ereturn display, `eform'
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
    // check some output is requested
    if "`graph'"=="nograph" & "`saving'"=="" & "`clear'"=="" & "`list'"=="" {
        di as error "Please specify either list or saving() or clear with nograph"
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
    if "`saving'"=="" tempfile saving
    postfile `post' type delta b se using `saving', `replace'
    if "`senstype'"=="both" local typelist 2
    else if "`senstype'"=="one" local typelist 1 3
    else local typelist 1 2 3
    foreach del of numlist `delta' {
        di "." _c
        foreach type in `typelist' {
            if `type'==1 local deltavar cond(`sens'==`randint',`log'(`del'),0)
            if `type'==2 local deltavar `log'(`del')
            if `type'==3 local deltavar cond(`sens'==`randcon',`log'(`del'),0)
            qui `maincmd' delta(`deltavar') `debug'
            ereturn post `b' `V'
            post `post' (`type') (`log'(`del')) (_b[`sens']) (_se[`sens'])
        }
    }
    di
    postclose `post'

    use `saving', clear
    label def type 1 "`randlab1' only" 2 "both arms" 3 "`randlab0' only"
    label val type type
    if "`list'"=="list" list, `listoptions'
    if "`clear'"!="" {
        cap drop _deltagraph
        local deltagraph _deltagraph
    }
    else tempvar deltagraph
    if "`log'"=="log" | "`eform'"!="" {
        * want results output and graphed on e-scale
        gen `deltagraph'=exp(delta)
    }
    else gen `deltagraph'=delta
    label var `deltagraph' "Delta, staggered for graph"
    
    if "`graph'"!="nograph" {
        *** DRAW A GRAPH
        di "Drawing graph..."
        gen low = b-1.96*se
        gen upp = b+1.96*se
        local col1 = word("`colors'",1)
        local col2 = word("`colors'",2)
        local col3 = word("`colors'",3)
        if "`col1'"=="" local col1 blue
        if "`col2'"=="" local col2 purple
        if "`col3'"=="" local col3 red
        local lpattern1 = word("`lpatterns'",1)
        local lpattern2 = word("`lpatterns'",2)
        local lpattern3 = word("`lpatterns'",3)
        if "`cistyle'"!="line" { // confidence limits as rspikes
            if `stagger'<0 {
                qui sum `deltagraph', meanonly
                local stagger = (r(max)-r(min))/100
            }
            qui replace `deltagraph'=`deltagraph'-`stagger' if type==1
            qui replace `deltagraph'=`deltagraph'+`stagger' if type==3
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
                (line b `deltagraph' if type==1, lcol(`col1') `lwidth' `lpattern1') 
                (rspike low upp `deltagraph' if type==1, lcol(`col1') `lwidth')
                (line b `deltagraph' if type==2, lcol(`col2') `lwidth' `lpattern2') 
                (rspike low upp `deltagraph' if type==2, lcol(`col2') `lwidth')
                (line b `deltagraph' if type==3, lcol(`col3') `lwidth' `lpattern3') 
                (rspike low upp `deltagraph' if type==3, lcol(`col3') `lwidth')
                ,
                `legendopt'
                ytitle(Coefficient of `sens')
                xtitle(`deltaname' in specified arm(s))
                xscale(`log')
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
                (line b low upp `deltagraph' if type==1, lcol(`col1' `col1' `col1') `lwidth' `lpattern') 
                (line b low upp `deltagraph' if type==2, lcol(`col2' `col2' `col2') `lwidth' `lpattern') 
                (line b low upp `deltagraph' if type==3, lcol(`col3' `col3' `col3') `lwidth' `lpattern') 
                ,
                `legendopt'
                ytitle(Coefficient of `sens')
                xtitle(`deltaname' in specified arm(s))
                xscale(`log')
                `gphoptions'
            ;
            #delimit cr
        }
        `ifdebug' "*** Running: `graphcmd'"
        `graphcmd'
        if "`clear'"!="" {
             global F9 `graphcmd'
             di as text "Graph command stored in F9"
             restore, not
        }
    }
    ereturn clear // Nothing sensible to ereturn
}
end

**********************************************************************

prog def pmm_reg
version 10
syntax varlist(min=1) [if] [in], [delta(string) robust debug noCONStant b(string) V(string)]

// PARSE
marksample touse, novarlist
gettoken y xlist : varlist
if "`debug'"=="" local qui qui

// IMPUTATION MODEL
`qui' di as text "*** Imputation model ***"
`qui' reg `y' `xlist' if `touse', `robust' `constant'
local factorI = e(N) / e(df_r)
mat `b'=e(b)
mat `v'=e(V)

// ANALYSIS MODEL
tempvar mdz
if "`delta'"!="" {
    qui gen `mdz' =  mi(`y') * `delta' if `touse'
    `qui' di as text "*** Analysis model ***"
    `qui' reg `mdz' `xlist' if `touse', robust `constant'
    local factorA = e(N) / e(df_r)
    mat `b' = `b' + e(b)
    mat `v' = `v' + e(V)*`factorI'/`factorA'
}

* sureg fails because it requires the same obs for both regns (and the same weights)
* but I verified that the residuals are exactly uncorrelated

end

**********************************************************************

prog def pmm_glm
version 10

/*******************************************
NOTES

Avoid -predict, residual- which uses unexpected formulae.

Stata's robust se's use a fudge factor of n/(n-p) for regress, n/(n-1) for ML commands.
I've done this with n = complete cases. 
*******************************************/

// PARSE
syntax anything [if] [in], delta(string) [debug DFCORRection(int -999) noCONStant b(string) V(string)]
gettoken cmd varlist : anything
unabcmd `cmd'
local cmd = r(cmd)
gettoken y xlist : varlist
if "`debug'"=="" local qui qui
if `dfcorrection'<0 & `dfcorrection'!=-999 di as error "Warning: negative dfcorrection"

// SET UP
marksample touse, novarlist
tempvar rowmiss id ri pi ra pa ystar offsetvar rai r2
tempname bi Vi Va
gen `id'=_n

// COUNT OBS
qui count if `touse'
local n = r(N)

// IMPUTATION MODEL
`qui' di as text _new "*** Imputation model ***"
if "`cmd'" != "regress" {
    qui gen `offsetvar' = missing(`y')*`delta' if `touse'
    local offsetopt offset(`offsetvar')
}
`qui' `cmd' `y' `xlist' if `touse', `offsetopt' `constant'
qui predict `pi' if `touse'
if "`cmd'" == "regress" {
    qui replace `pi' = `pi' + missing(`y')*`delta' if `touse'
}
qui gen `ri' = cond(mi(`y'), 0, `y'-`pi') if `touse'
mat `bi' = e(b)
local p = colsof(`bi')

// DF CORRECTION
if `dfcorrection'==-999 local dfcorrection = cond("`cmd'"=="regress",`p',1)
local factor = e(N) / (e(N) - `dfcorrection')

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
qui gen `ystar' = cond(missing(`y'),`pi',`y') if `touse'
`qui' di as text _new "*** Analysis model ***"
`qui' `cmd2' `ystar' `xlist' if `touse', `opts' `constant'
qui predict `pa' if `touse'
qui gen `ra' = (`ystar' - `pa') if `touse'
mat `b' = e(b)


// CONSTRUCT C MATRIX
sort `id'
qui gen `rai'=sqrt(`ri'*`ra') if `touse'
mat opaccum CII = `xlist' if `touse', group(`id') opvar(`ri') `constant'
mat opaccum CAA = `xlist' if `touse', group(`id') opvar(`ra') `constant'
mat opaccum CAI = `xlist' if `touse', group(`id') opvar(`rai') `constant'
mat C = (CAA, CAI \ CAI, CII)

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
mat opaccum BAA = `xlist' if `touse', group(`id') opvar(`opAA') `constant'
mat opaccum BAI = `xlist' if `touse', group(`id') opvar(`opAI') `constant'
mat BAI = -BAI
mat opaccum BII = `xlist' if `touse', group(`id') opvar(`opII') `constant'
mat BIA = J(`p',`p',0)
mat B = (BAA, BAI \ BIA, BII) 

if "`debug'"=="debug" {
    mat l B
    mat l C
}

// CALCULATE V MATRIX
mat Binv = inv(B)
mat fullV = Binv * C * Binv' * `factor'
mat `v' = fullV[1..`p',1..`p']
if "`debug'"=="debug" {
    mat l Binv
    mat l fullV
}
mat coleq `b' = ""
mat rownames `b' = `y'
mat roweq `v' = ""
mat coleq `v' = ""

end

**********************************************************************

prog def sm_ipw
version 10
syntax anything [if] [in], [delta(string) robust debug noSUMwt savewt(string) noCONStant b(string) V(string)]

// PARSE
marksample touse, novarlist
gettoken cmd varlist : anything
unabcmd `cmd'
local cmd = r(cmd)
gettoken y xlist : varlist
if "`debug'"=="" local qui qui

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
`qui' di _new as text "*** Fitting missingness model ***"
tempvar miss offset lp1 lp2 weight
qui gen `miss' = mi(`y') if `touse'
qui gen `offset' = cond(`miss',0,`delta'*`y')
qui ml model lf mnar_mml (`miss' = `xlist', offset(`offset') `constant'), vce(robust)
`qui' ml maximize
qui predict `lp1'

// FIT MAR MISSINGNESS MODEL
`qui' di _new as text "*** Fitting MAR missingness model ***"
qui ml model lf mnar_mml (`miss' = `xlist', `constant'), vce(robust)
`qui' ml maximize
qui predict `lp2'

// GET WEIGHTS
qui gen `weight' = (1 + exp(`lp1')) / (1 + exp(`lp2')) if `touse'
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
`qui' di _new as text "*** Fitting weighted analysis model ***" 
qui `cmd' `varlist' if `touse' [pw=`weight'], `constant'
mat `b'=e(b)
mat `v'=e(V)

// OPTIONALLY SAVE WEIGHT
if "`savewt'"!="" rename `weight' `savewt'
end

**********************************************************************

