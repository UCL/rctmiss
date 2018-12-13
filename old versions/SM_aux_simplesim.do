* SM_aux_simplesim: simple simulation with auxiliary variable in selection model
* should be easy to build into code
* IW 4nov2016

adopath ++ C:\Users\ian\Dropbox\ado\rctmiss

* Simulate data with trt effect of 1 fully mediated by w
clear
set seed 47160
set obs 100000
gen z = runiform()<0.5
gen w = z+rnormal()
gen ytrue = w+rnormal()
gen r = runiform()<invlogit(w-1)
gen y = ytrue if r

* analysis before data deletion
reg ytrue z

* CCA
reg y z 

* SM with incorrect weights - IMP in wrong selection model (W omitted)
rctmiss, smdelta(-1): reg y z

* SM with correct weights in wrong analysis 
rctmiss, smdelta(-1) savewt(wt): reg y w z

* SM with correct weights in right analysis 
reg y z [pw=wt]
