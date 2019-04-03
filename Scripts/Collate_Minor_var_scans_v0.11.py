#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess
import os
import math
import pandas as pd
import numpy as np
import re

from optparse import OptionParser
usage = "Usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-p", "--path",action="store", type="string", dest="inpath",help="Specify path of .bcf.analysed.tsv files", default='./')
parser.add_option("-l", "--list",action="store", type="string", dest="inlist",help="Specify filename of list with sequence names to be collated (line separated)")
parser.add_option("-r", "--reflength",action="store", type="string", dest="reflength",help="Specify length of reference all files are mapped to [%default]",default=1139780)

parser.add_option("-o", "--outprefix",action="store", type="string", dest="outprefix",help="Specify prefix for naming output files [%default]", default='')



parser.add_option("--mincov",action="store", type="string",dest="mincov",help="Specify minimum read coverage to keep 100% variant (minors already called separately) [%default]", default=8)
parser.add_option("--minfreq",action="store", type="string",dest="minfreq",help="Specify minimum percentage to call minor variant [%default]", default=5)
parser.add_option("--biasthreshold",action="store", type="string",dest="biasthreshold",help="Specify p value below which minor variants will fail bias tests [%default]", default=0.01)

parser.add_option("--ratio",action="store_true", dest="freq2ratio", help="Variant frequencies will output as a ratio instead of a percentage", default=False)

(options, args) = parser.parse_args()

# Read in the seqlist
myseqlist = []    
with open(options.inlist) as f:
    for line in f:
        line = line.strip('\n')
        myseqlist.append(line)

myseqlistfull = myseqlist
#myseqlistfull.insert(0, 'POS')

mypath = options.inpath


percratio=[]
if options.freq2ratio == True:
    percratio=100
else:
    percratio=1

reflength = int(options.reflength)

mincov=int(options.mincov)
minfreq=int(options.minfreq)
biasthreshold=int(options.biasthreshold)



print "Reading in sequences from " + options.inlist

# Bring in info for all sequences on whether a minor variant is present/absent at each site
# This is a new version that also brings in the alleles into another df so we can check if it's always the same alt allele

all_minors_freq = pd.DataFrame(range(1,reflength),columns=['POS'])
all_minors_ALTallele = pd.DataFrame(range(1,reflength),columns=['POS'])
all_minors_POS = pd.DataFrame()
all_minors_rawALT = pd.DataFrame(range(1,reflength),columns=['POS']) #####
all_minors_rawFreq = pd.DataFrame(range(1,reflength),columns=['POS']) #####
all_minors_Coverage = pd.DataFrame(range(1,reflength),columns=['POS']) #####



minorvars1 = pd.read_table(mypath + myseqlist[1] + '.bcf.analysed.tsv',sep='\t', header=0)
for currentseq in myseqlist:
    print 'Reading in ' + mypath + currentseq + '.bcf.analysed.tsv'
    minorvars1 = pd.read_table(mypath + currentseq + '.bcf.analysed.tsv',sep='\t', header=0)

# Remove duplicate positions (where there are more than one var calls per position - e.g. complex calls - always removes the later one - may want to improve this)
    minorvars1.drop_duplicates(subset='POS', keep='first', inplace=True)
    minorvars1.reset_index(inplace=True,drop=True)

    # Adjust positions to allow for 'missing' positions (e.g. zero coverage)
    minorvars1_length = pd.DataFrame(range(1,reflength),columns=['POS'])

    minorvars1 = pd.merge(minorvars1, minorvars1_length,on='POS',how='outer')
    minorvars1.sort_values(by=['POS'],inplace=True)
    minorvars1.reset_index(inplace=True,drop=True)

    minorvars1_POS = minorvars1.loc[:,['POS']]
    minorvars1_POS = minorvars1_POS.rename(columns={'POS': currentseq})
    all_minors_POS = pd.concat([all_minors_POS[:], minorvars1_POS[:]], axis=1)


# Check consensus variants first - make either 0 or 100 (REF or ALT)
    minorvars1['PassFreq'] = np.where((minorvars1['INDEL'] == False) &
                                      (minorvars1['MinReadSup'] == "ok") &
                                      (minorvars1['AltPerc'] > 50) &
                                      (minorvars1['DP'] >= mincov), 100,0)


    # Check if variant passes criteria for calling as a minor var - if not, leave as consensus
    minorvars1['PassFreq'] = np.where((minorvars1['StrandBias'] >= biasthreshold) &
                                      (minorvars1['BaseQBias'] >= biasthreshold) &
                                      (minorvars1['MapQBias'] >= biasthreshold) &
                                      (minorvars1['TailDistBias'] >= biasthreshold) &                                
                                      (minorvars1['INDEL'] == False) &
                                      (minorvars1['MinReadSup'] == "ok") &
                                      (minorvars1['MinorPerc'] >= minfreq),minorvars1['AltPerc'] ,minorvars1['PassFreq'])


####    
    minorvars1_freq = pd.DataFrame(minorvars1['POS'],columns=['POS'])
    if percratio == 100:
        minorvars1_freq['PassFreq'] = minorvars1.loc[:,['PassFreq']]/100
    else:
        minorvars1_freq['PassFreq'] = minorvars1.loc[:,['PassFreq']]
    minorvars1_freq = minorvars1_freq.rename(columns={'PassFreq': currentseq})

    all_minors_freq = pd.merge(all_minors_freq, minorvars1_freq, on='POS',how='outer')
    all_minors_freq.sort_values(by='POS',inplace=True)
    all_minors_freq.reset_index(inplace=True,drop=True)
    all_minors_freq = all_minors_freq.replace(np.nan, 0) # Deal with errors/nans


    
    ########### Now do the same for the alleles (i.e. pull out the possible ALT alleles using similar criteria)
   # Check consensus variants first - make either '-' (REF) or ALT (based on consensus)
    minorvars1['ALTallele'] = np.where((minorvars1['INDEL'] == False) &
                                      (minorvars1['MinReadSup'] == "ok") &
                                      (minorvars1['AltPerc'] > 50) &
                                      (minorvars1['DP'] >= mincov), minorvars1['ALT'],'-')

    # If a site has insufficient coverage, make it 'n'
    minorvars1['ALTallele'] = np.where((minorvars1['INDEL'] == False) &                                     
                                      (minorvars1['DP'] < mincov), 'n',minorvars1['ALTallele'])    
    
    # Check if variant passes criteria for calling as a minor var - if not, leave as consensus
    minorvars1['ALTallele'] = np.where((minorvars1['StrandBias'] >= biasthreshold) &
                                      (minorvars1['BaseQBias'] >= biasthreshold) &
                                      (minorvars1['MapQBias'] >= biasthreshold) &
                                      (minorvars1['TailDistBias'] >= biasthreshold) &                                
                                      (minorvars1['INDEL'] == False) &
                                      (minorvars1['MinReadSup'] == "ok") &
                                      (minorvars1['MinorPerc'] >= minfreq),minorvars1['ALT'] ,minorvars1['ALTallele'])  
    
    minorvars1_ALTallele = pd.DataFrame(minorvars1['POS'],columns=['POS'])
    minorvars1_ALTallele['ALTallele'] = minorvars1.loc[:,['ALTallele']]
    minorvars1_ALTallele = minorvars1_ALTallele.rename(columns={'ALTallele': currentseq})
    all_minors_ALTallele = pd.merge(all_minors_ALTallele, minorvars1_ALTallele, on='POS',how='outer')
    all_minors_ALTallele.sort_values(by='POS',inplace=True)
    all_minors_ALTallele.reset_index(inplace=True,drop=True)
        
#######
  # Extract 'raw' AltPerc aswell (to allow recalling of failed sites in a multisample analysis)    
    minorvars1_rawFreq = minorvars1.loc[:,['POS','AltPerc']] ###
    minorvars1_rawFreq = minorvars1_rawFreq.rename(columns={'AltPerc': currentseq}) ###
    all_minors_rawFreq = pd.merge(all_minors_rawFreq, minorvars1_rawFreq, on='POS',how='outer')
    all_minors_rawFreq.sort_values(by='POS',inplace=True)
    all_minors_rawFreq.reset_index(inplace=True,drop=True)   
    all_minors_rawFreq = all_minors_rawFreq.replace(np.nan, '0') # Deal with errors/nans
  # Extract 'raw' ALT allele (uncorrected)  
    minorvars1_rawALT = minorvars1.loc[:,['POS','ALT']] ###
    minorvars1_rawALT = minorvars1_rawALT.rename(columns={'ALT': currentseq}) 
    all_minors_rawALT = pd.merge(all_minors_rawALT, minorvars1_rawALT, on='POS',how='outer')###
    all_minors_rawALT.sort_values(by='POS',inplace=True)
    all_minors_rawALT.reset_index(inplace=True,drop=True)    
    all_minors_rawALT = all_minors_rawALT.replace(np.nan, '-') # Deal with errors/nans
  # Extract 'raw' Coverage/depth (to allow for coverage based analyses)    
    minorvars1_Coverage = minorvars1.loc[:,['POS','DP']] ###
    minorvars1_Coverage = minorvars1_Coverage.rename(columns={'DP': currentseq}) ###
    all_minors_Coverage = pd.merge(all_minors_Coverage, minorvars1_Coverage, on='POS',how='outer')###
    all_minors_Coverage.sort_values(by='POS',inplace=True)
    all_minors_Coverage.reset_index(inplace=True,drop=True)  
    all_minors_Coverage = all_minors_Coverage.replace(np.nan, '0') # Deal with errors/nans
#######
print 'Finished reading in files'


# If generating a 'ratio' instead of percentage, do this now (and create working file for downstream stuff) and round down values to maneagable amount 
freq2ratio = False
if freq2ratio == True:
#if options.freq2ratio == True:
    all_minors_freq = all_minors_freq.round(3) # round percentages down to a sensible number
else:
    all_minors_freq = all_minors_freq.round(1) # round percentages down to a sensible number

print ''
print 'All sequences read and compiled!'


#######
print ''
print 'Identifying invariant sites and generating list for subsetting tables'

# get rid of 'POS' column now we're merging rather than concat-ing
all_minors_freq_trim = all_minors_freq.loc[:,myseqlist]
# Identify sites that have a minor var (in any sequence), plus generate a 'position' column
all_minors_freq['Invariant'] = all_minors_freq_trim.eq(all_minors_freq_trim.iloc[:, 0], axis=0).all(1) # is whole row equal to first element?

mvar_sites = all_minors_freq.loc[all_minors_freq['Invariant']==False]
mvar_site_list = pd.DataFrame(mvar_sites['POS'],columns=['POS'])


print 'Identifying singleton sites'
Singleton = []
for currentrow in list(mvar_sites.index):
    Singleton.append(np.where((list((mvar_sites.loc[currentrow,myseqlist]!=0) & (mvar_sites.loc[currentrow,myseqlist]!=(100/percratio))).count(True))<=1, 'Singleton','Multiple'))
mvar_sites['Singleton'] = Singleton

print 'Identifying sites with only 2 frequences, e.g. 0/100'
# Identify sites where there are only two possible frequencies
numfreqs = []
for currentrow in list(mvar_sites.index):
    numfreqs.append(np.where(len(set((list(mvar_sites.loc[currentrow,myseqlist])))) <= 2, "2Freqs",'MultipleFreqs'))
mvar_sites['NumFreqs'] = numfreqs

####
print 'Identify sites where there are more than one ALT allele'
#all_minors_ALTallele = all_minors_ALTallele.replace('.',np.nan)
mvar_ALTallele = all_minors_ALTallele.loc[all_minors_freq['Invariant']==False]

# Check if there are more than one ALT allele at each site
Singleton_ALT = []
for currentrow in list(mvar_sites.index):
    Singleton_ALT.append(np.where(len(set(mvar_ALTallele.loc[currentrow,myseqlist].dropna()))==2, 'SingleALT','MultipleALT'))
mvar_ALTallele['Singleton_ALT'] = Singleton_ALT

mvar_sites['Singleton_ALT'] = mvar_ALTallele['Singleton_ALT']
####


print 'Identify unique alleles at each loci'
allelelist = []
for currentrow in list(mvar_ALTallele.index):
#    allelelist.append(",".join(set(list(all_minors_ALTallele.loc[currentrow,myseqlist]))))
    currentallele = ",".join(sorted(set(list(mvar_ALTallele.loc[currentrow,myseqlist]))))
    #currentallele = re.sub('n','',currentallele)
    #currentallele = re.sub('-','',currentallele)
    currentallele = re.sub(',$', '', currentallele)
    currentallele = re.sub('^,', '', currentallele)
    currentallele = re.sub(',,', ',', currentallele)
    allelelist.append(currentallele)
mvar_ALTallele['UniqueAlleles'] = allelelist
mvar_sites['UniqueAlleles'] = mvar_ALTallele['UniqueAlleles']


print 'Identify sites that include a minor var site (since code now captures 100% sites aswell)'
nonMinor2 = []
for currentrow in list(mvar_sites.index):
    nonMinor2.append(np.where((list((mvar_sites.loc[currentrow,myseqlist]!=0) & (mvar_sites.loc[currentrow,myseqlist]!=(100/percratio))).count(True))>0, 'Minors','None'))
mvar_sites['nonMinor'] = nonMinor2
mvar_only = mvar_sites.loc[mvar_sites['nonMinor']=='Minors']


myseqlistfull.insert(0, 'POS')
mvar_site_freq = mvar_sites.loc[:, myseqlistfull]
mvar_site_rawFreq = pd.merge(all_minors_rawFreq, mvar_site_list, on='POS',how='right')
#mvar_site_ALTallele = pd.merge(all_minors_ALTallele, mvar_site_list, on='POS',how='right')
mvar_site_rawALT = pd.merge(all_minors_rawALT, mvar_site_list, on='POS',how='right')
mvar_site_Coverage = pd.merge(all_minors_Coverage, mvar_site_list, on='POS',how='right')




print ''
print 'Analysis complete! Generating output files.'



# Output a full collated minor var file, plus one only containing multi-sample mvars
###
print ''
print 'Output full dataframe of sites with any variant present (unformatted): Collated_MinorVar_Any-Variant-Site_Full-FFreq-Table.tsv'
mvar_sites.to_csv(options.outprefix+"Collated_MinorVar_Any-Variant-Site_Full-FFreq-Table.tsv",sep='\t', index=False)
###
print 'Filtered sites and filtered frequencies - sites that contain minor vars (as opposed to having only 100% REF or ALT): Collated_MinorVar_Minority-Variant-Site_Filtered-Freqs-Table.tsv'
mvar_only.to_csv(options.outprefix+"Collated_MinorVar_Minority-Variant-Site_Filtered-Freq-Table.tsv",sep='\t', index=False)
print ''
###
print 'Filtered sites and filtered frequences where there is any type of variant : Collated_MinorVar_Any-Variant-Site_Filtered-Freqs.tsv'
mvar_site_freq.to_csv(options.outprefix+"Collated_MinorVar_Any-Variant-Site_Filtered-Freqs.tsv",sep='\t', index=False)
###
print 'Filtered sites with true variants in them, but with raw frequency percentages (less stringent) : Collated_MinorVar_Any-Variant-Site_Raw-Freqs.tsv'
mvar_site_rawFreq.to_csv(options.outprefix+"Collated_MinorVar_Any-Variant-Site_Raw-Freqs.tsv",sep='\t', index=False)
###
print 'Filtered sites and filtered ALT alleles where there is any type of variant : Collated_MinorVar_Any-Variant-Site_Filtered-ALTalleles.tsv'
#mvar_site_ALTallele.to_csv(options.outprefix+"Collated_MinorVar_Any-Variant-Site_Filtered-ALTalleles.tsv",sep='\t', index=False)
mvar_ALTallele.to_csv(options.outprefix+"Collated_MinorVar_Any-Variant-Site_Filtered-ALTalleles.tsv",sep='\t', index=False)
###
print 'Filtered sites and raw ALT alleles whre there is any type of variant : Collated_MinorVar_Any-Variant-Site_Raw-ALTalleles.tsv'
mvar_site_rawALT.to_csv(options.outprefix+"Collated_MinorVar_Any-Variant-Site_Raw-ALTalleles.tsv",sep='\t', index=False)
###
print 'Filtered sites and raw coverage/depth where there is any type of variant : Collated_MinorVar_Any-Variant-Site_Raw-Coverage.tsv'
mvar_site_Coverage.to_csv(options.outprefix+"Collated_MinorVar_Any-Variant-Site_Raw-Coverage.tsv",sep='\t', index=False)
###
print ''
print 'Finished'


