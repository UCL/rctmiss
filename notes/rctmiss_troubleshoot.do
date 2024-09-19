use "http://www.homepages.ucl.ac.uk/~rmjwiww/stata/missing/smoke.dta", clear

tab rand quit, miss
gen quit2 = quit
replace quit2 = 0 if missing(quit)
*MAR
logistic quit2 rand

*Sens analysis - same as MAR
rctmiss, pmmdelta(0, expdelta): logistic quit rand

*Sensitivity analysis based around missing=smoking, delta equal in both arms:
local n = 1
forvalues d = 0.1(0.1)1 {
	xi: rctmiss, pmmdelta(`d', expdelta): logistic quit i.rand
	*RENAME FACTOR VARIABLE
	local colnms: coln e(b)
	di "`colnms'"

	local colnms = subinstr("`colnms'","_Irand_1", "1.rand", .)
	mat b= e(b)
	mat colnames b= `colnms'

	erepost b= b, rename
	
	margins i.rand, expression(exp(predict(xb))/(1+exp(predict(xb))))
	mat e = r(table)
	local d`n'e_b = string(100*e[1,1], "%4.1f")
	local d`n'e_ll = string(100*e[5,1], "%4.1f")
	local d`n'e_ul = string(100*e[6,1], "%4.1f")
	di "Delta: `d': `d`n'e_b' (95% CI: `d`n'e_ll', `d`n'e_ul')"
	local n = `n'+1
}

*Sensitivity analysis based around missing=smoking, unequal delta by study arm:
local n = 1
forvalues d = 0.1(0.1)1 {
	xi: rctmiss, sens(i.rand, list nograph unequal) pmmdelta(0, expdelta): logistic quit i.rand
	*RENAME FACTOR VARIABLE
	local colnms: coln e(b)
	di "`colnms'"

	local colnms = subinstr("`colnms'","_Irand_1", "1.rand", .)
	mat b= e(b)
	mat colnames b= `colnms'

	erepost b= b, rename
	/*
	margins i.rand, expression(exp(predict(xb))/(1+exp(predict(xb)))) // produces error "last estimates not found"
	mat e = r(table)
	local d`n'e_b = string(100*e[1,1], "%4.1f")
	local d`n'e_ll = string(100*e[5,1], "%4.1f")
	local d`n'e_ul = string(100*e[6,1], "%4.1f")
	di "Delta: `d': `d`n'e_b' (95% CI: `d`n'e_ll', `d`n'e_ul')"
	*/
	local n = `n'+1
}
