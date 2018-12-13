*! version 0.2   16mar2010     rand() changed to sens(); gphoptions now added `loose'; new options debug robust lpattern() nograph; now calls rctmiss_*.ado not mnar_*.ado; various bug fixes
prog def rctmiss
version 10
/*************************************************************************
rctmiss: analyse a RCT allowing for informatively missing outcome data.

Notes:
    missing values of baseline covariates are imputed with the mean of their observed values (White & Thompson, 2005).
    user can choose a fixed analysis or a sensitivity analysis
    user can express departures from MAR in two ways
    - as the coefficient of y in a model for m on x and y (selection model): the coefficient is the increase in the log odds of missingness for a 1-unit increase in y (the log imor)
    - as the coefficient of m in a model for y on x and m (pattern mixture model): for linear regression, the coefficient is the difference between mean unobserved outcome and the mean observed outcome; for logistic regression, the coefficient is the log odds ratio between the outcome and missingness (the log imor) 

SYNTAX
    rctmiss, [sens(varname)] [smdelta(IM_expression)|pmmdelta(IM_expression)] [more_options] [graph_options]: regression_command
    
where more_options can be
    stagger(real -1) colors(string) saving(string) replace list LISTOPTions(string) nopreserve

IM_expression can be
    expression -> gives a single analysis
    numlist -> gives a sensitivity analysis; sens(varname) is required

STILL TO DO
    Allow smdelta() to be expresed via imor() or logimor()
    Alho option e.g. attempts(varname)
    Alternative handling of missing baselines: MIM, reg imputation
    Correct df for CIs
    help file

LIMITATIONS
    only 2 arms
    not (st)cox
*************************************************************************/

*** SEPARATE PREFIX AND REGRESSION COMMANDS ***

tokenize "`0'", parse(":")
if "`1'"==":" {
    local prefix
    local command `2'
}
else if "`2'"==":" {
    local prefix `1'
    local command `3'
}
else exit 198

*** PARSE REGRESSION COMMAND ***

gettoken regcmd restofcommand : command
unabcmd `regcmd'
local regcmd = r(cmd)

local 0 `restofcommand'
syntax varlist [if] [in] [fweight aweight iweight pweight], [robust *]
marksample touse, novarlist
gettoken yvar xvars: varlist
if index("`command'",",")==0 local comma ,
if "`robust'"=="robust" local robust1 robust
if "`weight'"!="" local weightexp [`weight'`exp']
local command `regcmd' `varlist' `if' `in' `weightexp', `options'
local restofcommand `varlist' `if' `in' `weightexp', `options'

*** PARSE PREFIX COMMAND ***
local 0 `prefix'
syntax, [sens(varname) SMDelta(string) PMMDelta(string) stagger(real -1) COLors(string) LWidth(passthru) LPATterns(string) saving(string) replace list LISTOPTions(string) nograph nopreserve debug robust *]
local gphoptions `options'
if "`robust'"=="" local robust `robust1'

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
        di as error "smdelta() or pmmdelta() must be specified if sens() is not specified"
        exit 498
    }
}
if "`smdelta'"!="" {
    * SELECTION MODEL / IPW METHOD
    local maincmd rctmiss_ipw `command'
    local deltaname SM delta
}
if "`pmmdelta'"!="" {
    * PATTERN MIXTURE MODEL / MEAN SCORE METHOD
    if "`regcmd'"!="regress" | "`robust'"=="robust" {
        local maincmd rctmiss_pmm `command'
    }
    else {
        local maincmd rctmiss_reg `restofcommand'
    }
    local deltaname PMM delta
}
local delta `smdelta'`pmmdelta'

// check some output is requested
if "`graph'"=="nograph" & "`saving'"=="" & "`list'"=="" {
    di as error "Please specify either list or saving() with nograph"
    exit 498
}
if "`preserve'"!="nopreserve" preserve


*** HANDLE INCOMPLETE BASELINES ***
if "`xvars'"!="" {
    foreach xvar of varlist `xvars' {
        qui count if mi(`xvar') & `touse'
        if r(N)>0 di as text "Imputing " as result r(N) as text " missing values of " as result "`xvar'" as text " with the observed mean"
        qui summ `xvar' if `touse', meanonly
        qui replace `xvar' = r(mean) if mi(`xvar') & `touse'
    }
}

// COUNT & REPORT OBS
qui count if `touse'
local n = r(N)
qui count if `touse' & !mi(`yvar')
local ncomplete = r(N)
local nincomplete = `n'-`ncomplete'    
if "`delta'"!="" {
	di as text "Results allowing for MNAR"
    di as text "Using " as result `ncomplete' as text " observed outcomes and " as result `nincomplete' as text " unobserved outcomes"
}
else {
	di as text "Results assuming MAR"
    di as text "Using " as result `ncomplete' as text " observed outcomes but ignoring " as result `nincomplete' as text " unobserved outcomes"
}

*** ANALYSIS ***
if "`sens'"=="" {
    * EXPRESSION SPECIFIED: SINGLE RUN
    if "`gphoptions'`colors'`lwidth'"!="" di as error "Options `gphoptions' `colors' `lwidth' ignored"
    dicmd `maincmd' delta(`delta') `debug'
}
else {
    * SENSITIVITY ANALYSIS
    di as text "Performing sensitivity analysis..."
    qui levelsof `sens' if `touse', local(randlevels)
    if wordcount("`randlevels'")>2 {
        di as error "Sorry, I can only handle two-arm trials at present"
        exit 498
    }
    if wordcount("`randlevels'")<2 {
        di as error "`sens' does not vary"
        exit 498
    }
    local randcon = word("`randlevels'",1)
    local randint = word("`randlevels'",2)
    local randlab0 : label (`sens') 0
    local randlab1 : label (`sens') 1
    tempname post
    if "`saving'"=="" tempfile saving
    postfile `post' type delta b se using `saving', `replace'
    forvalues type=1/3 {
        foreach del of numlist `delta' {
            if `type'==1 local deltavar cond(`sens'==`randint',`del',0)
            if `type'==2 local deltavar `del'
            if `type'==3 local deltavar cond(`sens'==`randcon',`del',0)
            qui `maincmd' delta(`deltavar') `debug'
            post `post' (`type') (`del') (_b[`sens']) (_se[`sens'])
        }
    }
    postclose `post'
    use `saving', clear
    label def type 1 "`randlab1' only" 2 "both arms" 3 "`randlab0' only"
    label val type type
    if "`list'"=="list" list, `listoptions'
    
    if "`graph'"!="nograph" {
        *** DRAW A GRAPH
        if `stagger'<0 {
            qui sum delta, meanonly
            local stagger = (r(max)-r(min))/100
        }
        qui replace delta=delta-`stagger' if type==1
        qui replace delta=delta+`stagger' if type==3
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
        if "`lpattern1'"!="" local lpattern1 lpattern(`lpattern1')
        if "`lpattern2'"!="" local lpattern2 lpattern(`lpattern2')
        if "`lpattern3'"!="" local lpattern3 lpattern(`lpattern3')
        #delimit ;
        twoway
            (line b delta if type==1, lcol(`col1') `lwidth' `lpattern1') 
            (rspike low upp delta if type==1, lcol(`col1') `lwidth')
            (line b delta if type==2, lcol(`col2') `lwidth' `lpattern2') 
            (rspike low upp delta if type==2, lcol(`col2') `lwidth')
            (line b delta if type==3, lcol(`col3') `lwidth' `lpattern3') 
            (rspike low upp delta if type==3, lcol(`col3') `lwidth')
            ,
            legend(order(1 3 5) label(1 "`randlab1' only") label(3 "both arms") label(5 "`randlab0' only") rows(1))
            ytitle(Coefficient of `sens')
            xtitle(`deltaname' in specified arm(s))
            `gphoptions'
        ;
        #delimit cr
    }
}
end

