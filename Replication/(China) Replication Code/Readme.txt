This folder contains MATLAB and Stata codes for replicating results in Section 8.2 of 'Testing for Time Stochastic Dominance' by Kyungho Lee, Oliver Linton, and Yoon-Jae Whang

1. China_replicate.m
- This MATLAB code implements time stochastic dominance testings Section 8.2.

2. Visib_Desc_Analysis.do
- This do file implements descriptive analysis and plotting in Section 8.2.

The rest MATLAB codes are functions for running China_phAll.m

3. TSD_contact_China.m
- Implement Contact-set approach

4. Lambda.m
- Calculate Lambda 

5. contact_set_estimation.m
- Estimate contact set

6. operation.m
- Calculate cumulated EDFs

7. operation_T.m
- Calculate cumulated EDFs of the terminal period

8. PowerSet.m
- Calculate power set

9. test_data.m
- Arrange data to be fit with testing

These are descriptions about data:

1. NCDC_panel0820.csv
- Data from Almond and Zhang (2021)

2. NCDC_panel_quarterly.csv
- Quarterly panel data constructed from NCDC_panel0820.csv

3. quarterly_visib.dta
- dta file of NCDC_panel_quarterly.csv

4. China_CNT_1028.csv
- Testing Results
