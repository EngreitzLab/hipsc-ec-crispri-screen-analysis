#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:04:03 2024

@author: chloemc
"""

import pandas as pd
import numpy as np
import scipy.stats 

####################
# EDIT INPUTS HERE #
####################

# Insert prefix in geneFrequency file used for non-targeting and safe guides
nontargeting_name = "NON-TARGETING"
safetargeting_name = 'SAFE'

# Insert names of the conditions you would like to compare
Condition1 = 'D3pos'
Condition2 = 'D3neg'

# Number of replicates (This code only works if you have two bioreplicates)
Replicates = 2

# Index of columns for each condition/bio-replicate to compare
# Ex. If i want to compare D3pos vs D3neg for bio-replicate 1, its indexes are
# [4, 2] in geneFrequency.csv

columns_to_compare = [[4,2],[5,3]]

# Name for output files
outputName = 'D3EC_vs_D3nonEC'



#########################################
# CALCULATE LOG2FC FOR EACH GENE TARGET #
#########################################

# Load in dataframe containing frequency of each guide across samples
data = pd.read_csv('geneFrequency.csv', index_col=0)

# Subset out safe and non-targeting (aka control guides)
non_df = data.loc[nontargeting_name]
safe_df = data.loc[safetargeting_name]
safe_avg = list(safe_df)
non_avg = list(non_df)

# Calculate fold change in control guides (BR1 D3+ EC vs D3- non-EC)
d3_safe_control_b1 = safe_avg[columns_to_compare[0][0]] / safe_avg[columns_to_compare[0][1]]
d3_safe_control_b2 = safe_avg[columns_to_compare[1][0]] / safe_avg[columns_to_compare[1][1]]

# Calculate fold change in control guides (BR2 D3+ EC vs D3- non-EC)
d3_non_control_b1 = non_avg[columns_to_compare[0][0]] / non_avg[columns_to_compare[0][1]]
d3_non_control_b2 = non_avg[columns_to_compare[1][0]] / non_avg[columns_to_compare[1][1]]

# Find average fold change of control guides within each bioreplicate
d3_b1_control = np.mean([d3_safe_control_b1,d3_non_control_b1])
d3_b2_control = np.mean([d3_safe_control_b2,d3_non_control_b2])

# Calculate fold change for all non-control guides for each bioreplicate (D3+ EC vs D3- non-EC)
D3_B1 = pd.DataFrame()
D3_B2 = pd.DataFrame()
names = [Condition1 + '1', Condition2 + '1',Condition1 + '2', Condition2 + '2']
D3_B1['original'] = np.divide(data[names[0]] , data[names[1]])
D3_B2['original'] = np.divide(data[names[2]] , data[names[3]])

# Normalize fold change to control guides
D3_B1['Ctrl_Normalized'] = D3_B1['original'] / d3_b1_control
D3_B2['Ctrl_Normalized'] = D3_B2['original'] / d3_b2_control

# Calculate log2 fold change 
D3_B1 = D3_B1.astype(float)
D3_B2 = D3_B2.astype(float)
D3_B1['log'] = np.log2(D3_B1['Ctrl_Normalized'])
D3_B2['log'] = np.log2(D3_B2['Ctrl_Normalized'])

# Create dataframe with fold change results
finalLogResults = pd.DataFrame()
finalLogResults['BR1'] = D3_B1['log']
finalLogResults['BR2'] = D3_B2['log']

finalLogResults.index = list(data.index)
finalLogResults.to_csv('Log2FC_byGene' + outputName + '.csv')



###################################
# CALCULATE LOG2FC FOR EACH GUIDE #
###################################

# Load in dataframe containing frequency of each guide across samples
data = pd.read_csv('guideFrequency.csv', index_col=0)

# Subset out safe and non-targeting (aka control guides)
non_df = data.loc[nontargeting_name]
safe_df = data.loc[safetargeting_name]
safe_avg = list(np.mean(safe_df))
non_avg = list(np.mean(non_df))

# Calculate fold change in control guides (BR1 D3+ EC vs D3- non-EC)
d3_safe_control_b1 = safe_avg[columns_to_compare[0][0]] / safe_avg[columns_to_compare[0][1]]
d3_safe_control_b2 = safe_avg[columns_to_compare[1][0]] / safe_avg[columns_to_compare[1][1]]

# Calculate fold change in control guides (BR2 D3+ EC vs D3- non-EC)
d3_non_control_b1 = non_avg[columns_to_compare[0][0]] / non_avg[columns_to_compare[0][1]]
d3_non_control_b2 = non_avg[columns_to_compare[1][0]] / non_avg[columns_to_compare[1][1]]

# Find average fold change of control guides within each bioreplicate
d3_b1_control = np.mean([d3_safe_control_b1,d3_non_control_b1])
d3_b2_control = np.mean([d3_safe_control_b2,d3_non_control_b2])

# Calculate fold change for all non-control guides for each bioreplicate (D3+ EC vs D3- non-EC)
D3_B1 = pd.DataFrame()
D3_B2 = pd.DataFrame()
D3_B1['original'] = np.divide(data[names[0]] , data[names[1]])
D3_B2['original'] = np.divide(data[names[2]] , data[names[3]])

# Normalize fold change to control guides
D3_B1['Ctrl_Normalized'] = D3_B1['original'] / d3_b1_control
D3_B2['Ctrl_Normalized'] = D3_B2['original'] / d3_b2_control

# Calculate log2 fold change 
D3_B1 = D3_B1.astype(float)
D3_B2 = D3_B2.astype(float)
D3_B1['log'] = np.log2(D3_B1['Ctrl_Normalized'])
D3_B2['log'] = np.log2(D3_B2['Ctrl_Normalized'])

# Create dataframe with fold change results
finalLogResults = pd.DataFrame()
finalLogResults['BR1'] = D3_B1['log']
finalLogResults['BR2'] = D3_B2['log']

finalLogResults.index = list(data.index)
finalLogResults.to_csv('Log2FC_byGuide_' + outputName + '.csv')



##############################################
# CALCULATE ADJ P VALUE FOR EACH GENE TARGET #
##############################################

# Load in dataframe containing frequency of each guide across samples
data = pd.read_csv('guideFrequency.csv', index_col=0)

# Get frequency of all control guides, per bioreplicate
a = D3_B1.loc[safetargeting_name]
b = D3_B1.loc[nontargeting_name]
c = pd.concat([a,b])
ctrl_ratios_d3pos_1 = [x for x in c['Ctrl_Normalized']]
a = D3_B2.loc[safetargeting_name]
b = D3_B2.loc[nontargeting_name]
c = pd.concat([a,b])
ctrl_ratios_d3pos_2 = [x for x in c['Ctrl_Normalized']]
non_pos = ctrl_ratios_d3pos_1 + ctrl_ratios_d3pos_2

#Calculate p value for each gene
P_Values = []
genes = list(set(data.index))
for i in genes:
    # Locate all guides for gene "i"
    tempD3posB1 =D3_B1.loc[i]
    tempD3posB2 =D3_B2.loc[i]

    # Locate the normalized frequency for all guides for gene "i"
    d3pos_1 = [x for x in tempD3posB1['Ctrl_Normalized']]
    d3pos_2 = [x for x in tempD3posB2['Ctrl_Normalized']]

    # Combine bioreplicates
    d3pos = d3pos_1 + d3pos_2

    # Calculate p value
    P_Values.append(float(scipy.stats.mannwhitneyu(d3pos,non_pos).pvalue))

    # Calculate Benjamini-Hochberg adjusted p value

p_vals = pd.DataFrame()
p_vals['p'] = P_Values
p_vals['p_adj'] = P_Values
p_vals.index = genes

p_vals.to_csv('PValues_' + outputName + '.csv')
























