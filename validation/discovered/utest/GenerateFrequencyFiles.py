#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 16:29:55 2024

@author: chloemc
"""
# This code takes in the CRISPRi-FlowFish outputs and generates a dataframe
# that creates a csv file for the frequency/count of each guide.


import pandas as pd
import numpy as np

####################
# EDIT INPUTS HERE #
####################

# List all conditions (prefixes used in the byExperimentRep files)
samples = ['D0', 'D3neg', 'D3pos']

# Read in csv file of all guide names in your experiment
guidenamesFileName = 'guideNames.csv'

# Number of replicates
Replicates = 2



######################################################
# CREATE GUIDECOUNTS.CSV AND GUIDEFREQUENCY.CSV FILE #
######################################################
types = ['freq','counts']

# Read in csv file of all guide names in your experiment
guidenames = pd.read_csv(guidenamesFileName, header=None)

samples.insert(0, 'OligoID')

for j in types:
    countbyPCRrep = pd.DataFrame()
    for k in samples: 
        for i in range(1,Replicates + 1):
            
            name = 'byExperimentRep/' + k + '-Rep' + str(i) + '.bin_' + j + '.txt'
            try:
                # Read in each counts/frequency file
                temp = pd.read_table(name)
                
                # Find any missing guides that are not present in the file (interferes with sorting)
                missing = list ( set(guidenames[0]) - set(temp['OligoID']) )
                
                # Create a new row for each missing guide with nan values
                for guide in missing:
                    new = pd.DataFrame([guide, np.nan, np.nan, np.nan]).T
                    new.columns = samples
                    temp = pd.concat([temp, new], axis = 0, ignore_index=True)
                    temp.index = list(temp['OligoID'])
                
                # create dataframe for each gene's guides and sort alphabetically
                temp = temp.sort_values(by=['OligoID'])
                addThis = temp[['OligoID',k]]
                addThis.index = guidenames
                
                # Sort all files and create entire dataframe with all guides for all genes
                addThis.sort_index(inplace=True) 
                addThis.columns = [[k + str(i) + "_OligoID", k + str(i)]]
                countbyPCRrep = pd.concat([countbyPCRrep, addThis[ k + str(i) ]], axis=1)
                
            except:
                pass

    # Add label for the gene that corresponds with each guide
    gene = []
    for h in list(guidenames[0]):
        temp = h.split('_')
        gene.append(temp[0])
    countbyPCRrep.index = gene
    countbyPCRrep['Guide'] = list(guidenames[0])
    
    # Save results to csv file
    if j == 'counts':
        countbyPCRrep.to_csv('guideCounts.csv')
    else:
        countbyPCRrep.to_csv('guideFrequency.csv')
    
    
    
#######################################################
# CREATE GENECOUNTS.CSV AND GENEFREQUENCY.CSV FILE #
#######################################################

GeneCounts = pd.DataFrame()
GeneFreq = pd.DataFrame()

# Get list of all genes
genes = list(set(countbyPCRrep.index))
genes.sort()

totalconditions = (len(samples) - 1) * Replicates

# Find sum of all counts for each condition/bioreplicate
sums = np.sum(countbyPCRrep)[0:totalconditions]

for i in genes:
    # Sum all counts for each gene
    temp =countbyPCRrep.loc[i]
    temp_sum = np.sum(temp)
    temp_sum = temp_sum[0:totalconditions]

    # Find frequency of gene counts
    temp_freq = temp_sum / sums

    # append counts, frequency to dataframe
    GeneCounts[i] = temp_sum
    GeneFreq[i] = temp_freq

# Save results to csv file
GeneCounts = GeneCounts.T
GeneFreq = GeneFreq.T
GeneCounts.to_csv('geneCounts.csv')
GeneFreq.to_csv('geneFrequency.csv')
    
    