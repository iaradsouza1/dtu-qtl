#!/usr/bin/env python
# coding: utf-8

# In[37]:


import os
import glob
import re
from functools import reduce

import pandas as pd
import numpy as np


# In[65]:


# Get all human samples 
ann = pd.read_csv('data/meta/SraRunTable.txt')
samples = ann.loc[  ann['Organism'] == 'Homo sapiens', 'Run' ]

file = open('/home/iaradsouza/samples.csv', 'w')
for i in samples:
    file.write(f'{i}.ann.vcf\n')
file.close()


# In[42]:


# Function to process and filter SNPs in each VCF
def process_vcf(sample, path):
    
    vcf_file = str(path) + str(sample) + '.ann.vcf'
    
    with open(vcf_file) as file:

        res = {
            'chr':[],
            'pos':[],
            'REF':[],
            'ALT':[],
            'ID':[],
            sample:[]
        }

        chrm = [
            'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
            'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
            'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
            'chr22', 'chrX', 'chrY'
        ]

        for line in file:
            
            # Skip header
            if not (line.startswith('##')):
                line = line.split('\t')

                # Get only the positions on the standard chromosomes,
                # Get only SNPs
                # Get only the biallelic sites
                if line[0] in chrm and len(line[3]) == 1 and len(line[4]) == 1 and len(line[4].split(',')) == 1:
                    res['chr'].append(line[0])
                    res['pos'].append(line[1])
                    res['ID'].append(line[2])
                    res['REF'].append(line[3])
                    res['ALT'].append(line[4])
                    
                    # Define the genotype code:
                    # ref homozygous = 0
                    # heterozygous = 1
                    # alt homozygous = 2
                    if line[9].split(':')[0] == '0/0' or line[9].split(':')[0] == '0|0':
                        res[sample].append(int(0))
                    elif line[9].split(':')[0] == '0/1' or line[9].split(':')[0] == '0|1':
                        res[sample].append(int(1))
                    elif line[9].split(':')[0] == '1/1' or line[9].split(':')[0] == '1|1':
                        res[sample].append(int(2))
    
    df = pd.DataFrame(res)
    return(df)


# In[46]:


# Process all VCFs and gather them
res = []
for i in range(len(samples)):
    print(f'Processing sample: {samples[i]} - ({i+1}/{len(samples)})')
    res.append(process_vcf(samples[i], 'data/vcf/'))

df = reduce(lambda x, y: pd.merge(x, y, on=['chr', 'pos', 'REF', 'ALT', 'ID'], how='outer'), res)
df = df.fillna(0)


# In[50]:


# Save
if not os.path.isdir('results'):
    print(f'Creating "results/" directory at {os.getcwd()}')
    os.mkdir('results')
    
df.to_csv('results/genotypes.csv')

