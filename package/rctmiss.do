/* 
rctmiss.do: example of -rctmiss-, following help file

History:
25apr2019 syntax corrected - sens options become suboptions
12jan2017 created

Notes:
Start in C:\ado\ian\rctmiss\package
*/

// UK500 data (quantitative outcome)

use UK500, clear

* Analysis assuming MAR, dropping missing baselines:

reg sat96 rand sat94 i.centreid

* Analysis assuming MAR, with mean imputation for missing baselines
gen sat94fill = sat94
summ sat94
replace sat94fill = r(mean) if mi(sat94)
reg sat96 rand sat94fill i.centreid

* Same using rctmiss

xi: rctmiss, pmmdelta(0): reg sat96 rand sat94 i.centreid

* Single MNAR analysis, assuming missing values are 5 units lower than observed values in both arms:

xi: rctmiss, pmmdelta(-5): reg sat96 rand sat94 i.centreid

* Sensitivity analysis, assuming missing values are from 0 to 10 units 
* lower than observed values, in one arm or in both arms:

xi: rctmiss, sens(rand) pmmdelta(-10/0): reg sat96 rand sat94 i.centreid

* Improving appearance:

xi: rctmiss, pmmdelta(-10/0) ///
	sens(rand, legend(rows(3)) stagger(0.05) list(sepby(delta)) ///
		title(Sensitivity analysis for UK500 data)): ///
	reg sat96 rand sat94 i.centreid


// Smoking data (binary outcome)

use smoke, clear

tab rand quit, miss

* Analysis assuming missing = smoking:

gen quit2 = quit

replace quit2 = 0 if missing(quit)

logistic quit2 rand

* Same analysis using rctmiss:

rctmiss, pmmdelta(0, expdelta): logistic quit rand

* Sensitivity analysis based around missing=smoking:

rctmiss, sens(rand) pmmdelta(0(0.1)1, expdelta base(0)): logistic quit rand

