# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 16:16:49 2020

@author: Kyungho Lee, Oliver Linton, Yoon-Jae Whang

"""

"""

Plotting
Descriptive Statistics of Thai Panel Data

"""

import pandas as pd
import matplotlib.pyplot as plt

Thai = pd.read_csv("ThaiVillage_Panel.csv")
plt.style.use("seaborn-whitegrid")

###### Group by Descriptive Stat
Thai_yearly = Thai.groupby(["year","small"])
Thai_mean = Thai_yearly.mean()[["tc","netinc"]].reset_index()
Thai_sd = Thai_yearly.std()[["tc","netinc"]].reset_index()
Thai_q1 = Thai_yearly.quantile(0.25)[["tc","netinc"]].reset_index()
Thai_q2 = Thai_yearly.quantile(0.5)[["tc","netinc"]].reset_index()
Thai_q3 = Thai_yearly.quantile(0.75)[["tc","netinc"]].reset_index()


###########################
Thai_mean_0 = Thai_mean[Thai_mean["small"]==0]
Thai_mean_1 = Thai_mean[Thai_mean["small"]==1]
Thai_sd_0 = Thai_sd[Thai_sd["small"]==0]
Thai_sd_1 = Thai_sd[Thai_sd["small"]==1]
Thai_q1_0 = Thai_q1[Thai_q1["small"]==0]
Thai_q1_1 = Thai_q1[Thai_q1["small"]==1]
Thai_q2_0 = Thai_q2[Thai_q2["small"]==0]
Thai_q2_1 = Thai_q2[Thai_q2["small"]==1]
Thai_q3_0 = Thai_q3[Thai_q3["small"]==0]
Thai_q3_1 = Thai_q3[Thai_q3["small"]==1]

# Quantiles of TC
plt.figure(figsize=(7,5))
# plt.title("Total Consumption",fontsize=15)
plt.xlim(1997,2007)
plt.ylim(20000, 160000)
plt.xlabel("year",fontsize=15)
plt.ylabel("Total Consumption (THB)",fontsize=12)
plt.axvline(2001.5,color="black")
plt.fill_between(Thai_q1_0.year,Thai_q1_0.tc, Thai_q3_0.tc, alpha=0.2, color="blue")
plt.plot(Thai_q2_0.year, Thai_q2_0.tc,ls=":",marker="o",label="Big:Median", color = 'blue')
plt.plot(Thai_mean_0.year, Thai_mean_0.tc,marker="o",label="Big:Mean", color = 'blue')
plt.fill_between(Thai_q1_1.year,Thai_q1_1.tc, Thai_q3_1.tc, alpha=0.2,color="red")
plt.plot(Thai_q2_1.year, Thai_q2_1.tc,ls=":",marker = "s",label="Small:Median",color="red")
plt.plot(Thai_mean_1.year, Thai_mean_1.tc,marker="s",label="Small:Mean",color="red")
plt.legend(loc="upper left")
plt.savefig("Figure_ThaiA.png")
plt.show()


# Quantiles of NI
plt.figure(figsize=(7,5))
# plt.title("Net Income",fontsize=15)
plt.xlim(1997,2007)
plt.ylim(20000, 160000)
plt.xlabel("year",fontsize=15)
plt.ylabel("Net Income (THB)",fontsize=12)
plt.axvline(2001.5,color="black")
plt.fill_between(Thai_q1_0.year,Thai_q1_0.netinc, Thai_q3_0.tc, alpha=0.2,color="blue")
plt.plot(Thai_q2_0.year, Thai_q2_0.netinc,ls=":",marker="o",label="Big:Median",color="blue")
plt.plot(Thai_mean_0.year, Thai_mean_0.netinc,marker="o",label="Big:Mean",color="blue")
plt.fill_between(Thai_q1_1.year,Thai_q1_1.netinc, Thai_q3_1.tc, alpha=0.2,color="red")
plt.plot(Thai_q2_1.year, Thai_q2_1.netinc,ls=":",marker = "s",label="Small:Median",color="red")
plt.plot(Thai_mean_1.year, Thai_mean_1.netinc,marker="s",label="Small:Mean",color="red")
plt.legend(loc="upper left")
plt.savefig("Figure_ThaiB.png")
plt.show()


"""

Applying Static stochastic dominance testing on Thai Panel data
We use the python package 'pysdtest' for the testing.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysdtest

Thai = pd.read_csv("Thai/ThaiVillage_Panel.csv")

year_list, variable_list, s_list, sb_pval_list, bs_pval_list = [], [], [], [], []

for year in np.arange(2002,2008):
    for variable in ['netinc','tc']:
        for order in [1,2]:
        
            Thai_Panel_small = np.array(Thai[(Thai['year'] == year) & (Thai['small'] == 1)][variable].dropna())
            Thai_Panel_big = np.array(Thai[(Thai['year'] == year) & (Thai['small'] == 0)][variable].dropna())
            # Scaling
            Thai_Panel_big = Thai_Panel_big/1000000
            Thai_Panel_small = Thai_Panel_small/1000000

            testing_contact_sb = pysdtest.test_sd_contact(Thai_Panel_small, Thai_Panel_big, ngrid = 100, s = order, resampling = 'bootstrap')
            testing_contact_sb.testing()
            sb_pval_list.append(testing_contact_sb.result['pval'])
            
            testing_contact_bs = pysdtest.test_sd_contact(Thai_Panel_big, Thai_Panel_small, ngrid = 100, s = order, resampling = 'bootstrap')
            testing_contact_bs.testing()
            bs_pval_list.append(testing_contact_bs.result['pval'])
            
            year_list.append(year)
            variable_list.append(variable)
            s_list.append(order)
            Testing_result = pd.DataFrame({'year' : year_list,
             'variable': variable_list,
             's': s_list,
             'sb_pval': sb_pval_list,
             'bs_pval': bs_pval_list})
            Testing_result.to_excel("Testing_result.xlsx", index = False)
