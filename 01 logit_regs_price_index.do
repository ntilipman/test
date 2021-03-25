clear
cd "C:\Users\tilipman\Box Sync\Colorectal Cancer\Data and Documentation\Intermediate"
global outpath "C:\Users\tilipman\Box Sync\Colorectal Cancer\Output"
global matpath "C:\Users\tilipman\Box Sync\Colorectal Cancer\Code\Matlab"
set more off

/*##################################################################
Author: Nicholas Tilipman
Last Modified: 12/15/20
Purpose: This computes the logit price index
####################################################### */


capture program drop logit_index
program define logit_index 


local date = "200421"

global reg_list_eff1 eff1 side4 tablet second pfizer i.mkt

*Made this following data
use logit_log_regression.dta, clear

*Original order of observations
gen orig = _n


/* Defining variables */
replace price = price/1000


capture drop second
gen second = 0
replace second = 1 if eff1==9.5 | eff3==1.5 | (eff3>4.0 & eff3<4.1)

capture drop eff2second
gen eff2second = eff2*second

capture drop side5second
gen side5second = side5*second


egen prod = group(eff1 eff2 eff3 side1 side2 side3 side4 side5 tablet)

/*Fixes using sean's new numbers */

/* Cetiximab Fix */
replace eff1 = 6.9 if prod==1
replace eff1 = 8.6 if prod==2
/* Oxaliplatin + 5 FU LV */
replace eff1 = 19.7 if prod==9
replace eff2 = 46.8 if prod==9
/* Oxaliplatin + capecitabine */
replace eff1 = 19.7 if prod==8
/* Irinotecan */
replace eff1 = 9.6 if prod==3
replace eff2 = 14.3 if prod==3

replace side4 = 19 if prod==12


sort orig
outsheet second using "$matpath\second.out", replace
outsheet eff1 using "$matpath\eff1_corrected.out", replace
outsheet eff2 using "$matpath\eff2_corrected.out", replace



capture drop firm
gen firm = 1 if prod==4 
replace firm = 2 if prod==6 | prod==3
replace firm = 3 if prod==5 | prod==7
replace firm = 4 if prod==9 | prod==8
replace firm = 5 if prod==1 | prod==2
replace firm = 6 if prod==10 | prod==12 | prod==11

capture drop fu
gen fu = 0
replace fu = 1 if prod==4 | prod==6 | prod==9 | prod==10 | prod==11

capture drop pfizer
gen pfizer = 0
replace pfizer = 1 if prod==6 | prod==3 | prod==7 | prod==2 | prod==10


capture drop roche
gen roche = 0
replace roche=1 if prod==5 | prod==7 | prod==8 | prod==12


capture drop sanofi
gen sanofi = 0
replace sanofi = 1 if prod==9 | prod==8 | prod==11 | prod==12

capture drop imclone
gen imclone = 0
replace imclone = 1 if prod==1 | prod==2

capture drop genetech
gen genetech = 0
replace genetech = 1 if prod==10 | prod==12 | prod==11

sort obs
outsheet pfizer roche sanofi imclone genetech using "$matpath\firm_new.out", replace

/* Number of products */

gen obs = 1
capture drop numprods
bysort mkt: egen numprods = sum(obs)

/* Alternate Instruments */

bysort mkt firm: egen numprods_firm = sum(obs)

gen iv3 = numprods - numprods_firm

capture drop obs

tsset prod mkt

gen iv4 = l1.logprice

gen othermonths = l1.eff1
*replace othermonths = 0 if othermonths==.
bysort mkt: egen summonths = sum(othermonths)
gen diffmonths = summonths-othermonths

gen iv5 = diffmonths/(numprods-1)

bysort mkt firm: egen sumfirmmonths = sum(othermonths)
gen difffirmmonths = summonths-sumfirmmonths

gen iv6 = difffirmmonths/iv3

char mkt[omit] 51


capture drop const
gen const = 1
*OLS regression
xi: reg y logprice ${reg_list_eff1} const, noconst
estimates save "$outpath\ols_logprice_tab_second", replace

*IV regression
xi: ivreg y ${reg_list_eff1} const (logprice = iv1 iv2), noconst
estimates save "$outpath\logit_logprice_tab_second", replace

*generating residuals
predict resid, r
foreach var of varlist _I* {
	gen new_`var' = _b[`var']
	replace `var' = 0
	}
capture drop delta
*Y is the original mean utility (with time dummies)
gen delta = y
*yhat here is without time dummy (for outsidech below)
predict yhat, xb
*summing over delta (the orig mean util with time dummies, not the yhat)
gen expyhat = exp(delta)
bysort mkt: egen sumyhat = sum(expyhat)
gen shr2 = expyhat/(1+sumyhat)
*Saving the results file for drawing the index later
save tmp_logit.dta, replace


sort prod mkt
bysort prod: gen prod_obs = _n

*First residual observed
gen r_first = resid if prod_obs==1
bysort prod: egen max_r = max(r_first)

capture drop delta
gen delta = y-resid

*If keeping unobserved quality at first period obsevation, replacing residual at its first value
if "`2'"=="product0" {
	replace resid = max_r
	}

*If outside option changes over time, replacing mean utility with mean utility minus time dummy
if "`1'" == "outsidech" {
	replace delta = yhat
	}
	
replace delta = delta+resid
	
/* Time dummies defining */
capture drop dum
gen dum = 0
forvalues i = 1/50 {
	replace dum = new__Imkt_`i' if mkt==`i'
	}


replace price = logprice
*gen price = exp(logprice)/1000

capture drop theta1
gen theta1 = _b[logprice]



gen util = exp(delta)
gen maxutil = delta
gen expdum = exp(-dum)

replace price = price*shr2





collapse (max) maxutil (sum) util price (mean) dum expdum theta1 logprice, by(mkt)

local j = 1
local k = 0
if "`1'"=="outsidech" {
	local j = expdum
	local k = -dum
	}

/*If outside option remains the same, `j'=1 so add 1 to utility cuz normalize to 0. if outside option changes, `j'=exp(-time dummy) b/c this is the value of outside option for this period.  */
replace util = ln(util+`j')

gen cv = (util - util[_n-1])/-theta1

gen mut = cv/(cv + price)

local start1 = 1
local start2 = 2


gen pindex = 100 if _n==`start1'
forvalues i = `start2'/51 {
   replace pindex = pindex[_n-1]*(1-mut) if _n==`i'
   }
   
keep mkt pindex
save pindex_logit_`1'_`2'.dta, replace  

end

logit_index outside0 product0
logit_index outside0 productch
logit_index outsidech product0
logit_index outsidech productch



