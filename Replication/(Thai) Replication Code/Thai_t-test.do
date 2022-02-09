* t-test results
* authors: Kyungho Lee, Oliver Linton, Yoon-Jae Whang

import delimited ThaiVillage_Panel.csv, clear
tsset case_id year

ttest tc if year == 2002, by(small) 
ttest tc if year == 2003, by(small) 
ttest tc if year == 2004, by(small) 
ttest tc if year == 2005, by(small) 
ttest tc if year == 2006, by(small) 
ttest tc if year == 2007, by(small) 

ttest netinc if year == 2002, by(small) 
ttest netinc if year == 2003, by(small) 
ttest netinc if year == 2004, by(small) 
ttest netinc if year == 2005, by(small) 
ttest netinc if year == 2006, by(small) 
ttest netinc if year == 2007, by(small) 