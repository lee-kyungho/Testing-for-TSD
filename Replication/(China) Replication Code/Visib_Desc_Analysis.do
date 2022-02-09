clear

use NCDC_panel0820, clear

* Keeping only data after 2011 (following AZ)
keep if year>=2011

* defining quarter dummy
gen quarter = 1 if month == 1 | month == 2 | month == 3
replace quarter = 2 if month == 4 | month == 5 | month == 6
replace quarter = 3 if month == 7 | month == 8 | month == 9
replace quarter = 4 if month == 10 | month == 11 | month == 12

gen yyqq = year * 10 + quarter

* constructing quarterly panel data
collapse (mean) visib temp wdsp prcp year month quarter provgb citygb CNT, by(yyqq stn)

* time index
gen time = yq(year, quarter)
format time % tq
tsset stn time

* save data
save quarterly_visib, replace
export delimited NCDC_panel_quarterly.csv, replace

* data summary
use quarterly_visib, replace

gen phase = 0 if year >= 2011 & year < 2014
replace phase = 1 if year >= 2014 & yyqq < 20163
replace phase = 2 if yyqq >= 20163 

sort CNT phase
bysort CNT : su visib wdsp temp prcp 

**** Plotting ****
use quarterly_visib, clear

sort CNT
codebook stn if CNT == 1
codebook stn if CNT == 0

reg visib wdsp temp prcp i.quarter i.provgb
predict res_visib, res

collapse (mean) visib res_visib (sd) sd_visib = visib sd_res = res_visib (max) max_res = res_visib (min) min_res = res_visib , by(time CNT)
* save quarterly residuals data
save CNT_visib_res, replace
export delimited CNT_visib_res.csv, replace

* Plotting
use CNT_visib_res, replace

grstyle init
grstyle color background white
grstyle color major_grid gs8
grstyle linewidth major_grid thin
grstyle linepattern major_grid dot
grstyle yesno draw_major_hgrid yes
grstyle yesno grid_draw_min yes
grstyle yesno grid_draw_max yes
grstyle anglestyle vertical_tick horizontal
grstyle gsize axis_title_gap tiny
grstyle set legend 2, inside

* Mean Visib
tw (line visib time if CNT == 1) (line  visib time if CNT == 0, lp(dash)) , tline(`=q(2014q1)' `=q(2016q3)', lcolor(gs10)) legend(col(1) lab(1 "Regulated regions") lab(2 "Unregulated regions")) ylabel(6(1)15) ytitle("Visibility")
graph export "visibility-trend.png", replace

* Mean Residuals
tw (line res_visib time if CNT == 1) (line  res_visib time if CNT == 0, lp(dash)), tline(`=q(2014q1)' `=q(2016q3)', lcolor(gs10)) legend(col(1) lab(1 "Regulated regions") lab(2 "Unregulated regions") ) ylabel(-3(1)3) ytitle("Residuals")
graph export "residual-trend.png", replace
