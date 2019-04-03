#!/usr/bin/env python


import sys
import subprocess
import os
import math
import pandas as pd
import numpy as np
import ggplot as gg

# To generate new annotations for snpEff: -
#	Create folder in ~/programs/snpEff/data/ with name of genome
#	Create entry in ~/programs/snpEff/snpEff.config with details of genome (search for Treponema for examples)
# 	Deposit WGS fasta (labelled sequences.fa) and gff3 (labelled genes.gff) in folder.
# 	Ensure sequence headers inside file correspond to the one used by the target bcf/vcf files (can parse using one liner)
# 	Build database in snpEff, e.g. : java -jar ~/programs/snpEff/snpEff.jar build -gff3 -v Neisseria_gonorrhoeae_FA_1090

from optparse import OptionParser
#def get_options():
usage = "Usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", "--infile",action="store", type="string", dest="vcfinput",help="Specify input FILE")
parser.add_option("-r", "--reference",action="store", type="string", dest="reference",help="Specify name of primary reference [%default]", default='CP001050')
parser.add_option("-m", "--minreads",action="store",default=5, type="int",dest="minreads", help="Minimum number of reads required to call minor variant [%default]")
parser.add_option("-p", "--minperc",action="store", type="int",dest="minperc", default=5, help="Minimum percentage minor variant to accept [%default]")
parser.add_option("-t", "--covthreshold",action="store", type="int",dest="covthreshold", default=100, help="Minimum number of reads at any site to be included in analysis (variant and non variant) [%default]")
parser.add_option("-g", "--graphics", action="store_true", dest="graphics", default=False, help="Generate graphical output using sliding windows [%default]")
parser.add_option("-w", "--windowsize",action="store", type="int",dest="windowsize", default=10000, help="Size of sliding window for graphical output [%default]")
parser.add_option("-a", "--annotation",action="store",dest="annotation", default='None', help="Incorporate annotation using SnpEff (must provide a valid snpEff reference) [%default]")

(options, args) = parser.parse_args()


#	return parser.parse_args()
#	if len(args) < 1:
#        	parser.error("incorrect number of arguments - must specify input using '-i'\n")
#get_options()

####

vcfinput = options.vcfinput

#minreads = 5 # minimum number of reads to support (any) variant
minreads = options.minreads
#minfreq = 5 # minimum frequency (%) of variant to include in minor vars
minfreq = options.minperc
#covthreshold = 100
covthreshold = options.covthreshold
#windowsize = 10000 # size of sliding window for graphics
windowsize = options.windowsize
reference = options.reference
annotation = options.annotation
####




#arguments = str(sys.argv)
#vcfinput = sys.argv[1]


# input name of bcf file - will want to tweak this to take command line arguments and maybe work from a sam/bam aswell

#vcfinput = "/Users/mb29/testing_data/17723_4#7.bcf"
#vcfinput = "/Users/mb29/testing_data/17723_4#7.bcf"

vcfoutput = vcfinput + ".temp"
#print vcfoutput


# Use bcftools to extract the required columns into a tsv file
bcftools_path = 'bcftools-1.3'
snpeff_path = '/nfs/users/nfs_m/mb29/programs/snpEff/'


#bcftools-1.3 view $infile | java -jar ~/programs/snpEff/snpEff.jar Treponema_pallidum_DAL_1_uid87065 -c ~/programs/snpEff/snpEff.config -csvStats -lof -onlyProtein - > $infile.ann.vcf


print '\nExtracting relevant columns from '+ vcfinput +' into '+ vcfoutput

#if options.annotation == 'None':
#	bcfcommand = [bcftools_path,"query","-f",'%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t%INFO/PV4\n',"-o",vcfoutput,vcfinput]
#else:
	


def bcfcommand():
	if options.annotation == 'None':
		bcfrun = [bcftools_path,"query","-f",'%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t%INFO/PV4\n',"-o",vcfoutput,vcfinput]
		subprocess.call(bcfrun)
	else:	
		vcfinputanot=vcfinput +".ann.vcf"
		os.environ["PYTHONPATH"] = snpeff_path
		snpeffjar = snpeff_path + "snpEff.jar"
		snpeffconf = snpeff_path + "snpEff.config"
#		vcfinputanot=vcfinput +".ann.vcf"
		with open(vcfinputanot,'a') as f1:
			bcfview = subprocess.Popen([bcftools_path,'view',vcfinput],stdout=subprocess.PIPE)
#			runsnpeff = subprocess.Popen(['java','-jar',str(snpeffjar),'-c',str(snpeffconf),'-no-downstream','-no-intergenic','-no-intron','-no-upstream','-no-utr',annotation],stdin=bcfview.stdout,stdout=subprocess.PIPE)
			runsnpeff = subprocess.Popen(['java','-jar',str(snpeffjar),'-c',str(snpeffconf),'-no-downstream','-no-intergenic','-no-intron','-no-upstream','-no-utr','-no','INTRAGENIC',annotation],stdin=bcfview.stdout,stdout=subprocess.PIPE)

#		runsnpeff = subprocess.Popen(['java','-jar',str(snpeffjar),'-c',str(snpeffconf),annotation],stdin=bcfview.stdout,stdout=[subprocess.PIPE,vcfan])
#		printsnpeff = subprocess.Popen(['java','-jar',str(snpeffjar),'-c',str(snpeffconf),annotation,'>',str(vcfinputanot) ],stdin=bcfview.stdout)
			
#			bcfview.stdout.close()
#                vcfan = open(vcfinputanot, "w")
#		vcfan = runsnpeff.stdout
			
#			bcfquery = subprocess.Popen([bcftools_path,"query","-f",'%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t%INFO/PV4\t%INFO/ANN\n',"-o",vcfoutput],stdin=runsnpeff.stdout)
			bcfquery = subprocess.Popen([bcftools_path,"query","-f",'%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t%INFO/PV4\t%INFO/ANN\n',"-o",vcfoutput],stdin=subprocess.PIPE,stdout=f1)
			runsnpeff_out=runsnpeff.communicate()[0]
#			runsnpeff.stdout.close()
#			output = bcfquery.communicate()[0]
			bcfquery.communicate(runsnpeff_out)
			runsnpeff.wait()
			bcfquery.wait()
			bcfview.stdout.close()
			runsnpeff.stdout.close()
#			bcfquery.stdout.close()	

#	bcfcommand = [bcftools_path,"view",vcfinput,"|","java","-jar",snpeff_path,"snpEff.jar",options.annotation,"-c",snpeff_path,"snpEff.config","-csvStats","-onlyProtein","-","|",bcftools_path,"query","-f",'%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\t%INFO/PV4\t%INFO/ANN\n',"-o",vcfoutput]

#print bcfcommand 
#subprocess.call(bcfcommand)




# need to restore this to make the script work again
bcfcommand()




#for temporary work on a smaller file
#vcfoutput = '/Users/mb29/testing_data/17723_4#7_3000.tsv'


### Now pull the tsv into the python script
vcfdf = pd.read_table(vcfoutput,sep='\t', header=None,names=['CHROM','POS','REF','ALT','DP','DP4','PV4','ANN'])
### Get rid of temporary output
#subprocess.call(["rm",vcfoutput])


#binsize = 2 # size of frequency bin (not currently used)

print 'Analysing '+ vcfoutput +' using reference ' +reference
print 'Calling minor variants >=' + str(minfreq) + '% with >=' + str(minreads) + ' supporting reads at sites with >=' + str(covthreshold) + ' reads'


# Label positions where there is a potential variant
vcfdf['PossVariant'] = np.where(vcfdf['ALT']=='.', 'None', 'Variant')


# Identify INDELS
def isINDEL(myrows):
    INDEL = []
    for currentrow in myrows:
        reflength = len(vcfdf.loc[currentrow]['REF'])
        altlength = len(vcfdf.loc[currentrow]['ALT'])
        if reflength >1 or altlength > 1:
            INDEL.append('TRUE')
        else:
            INDEL.append('FALSE')
    return INDEL
vcfdf['INDEL'] = isINDEL(range(len(vcfdf)))


# Expand 'DP4' columnn into 4 separate columns, then add back in
DP4 = pd.DataFrame(vcfdf['DP4'].map(eval).tolist())
DP4.columns = ['RF','RR','AF','AR']
# Remerge with vcfdf
vcfdf = pd.concat([vcfdf[:], DP4[:]], axis=1)
#vcfdf


# Split up the PV4 (tests of bias)
# Deal with columns without any values by generating '1' four times, else use PV4 values
vcfdf['PV4'] = np.where(vcfdf['PV4']=='.', '1,1,1,1', vcfdf['PV4'])

# Now split into columns and remerge (as previously)
PV4 = pd.DataFrame(vcfdf['PV4'].map(eval).tolist())
PV4.columns = ["StrandBias","BaseQBias","MapQBias","TailDistBias"]

vcfdf = pd.concat([vcfdf[:], PV4[:]], axis=1)


# Calculate supporting reads from each side
vcfdf['RefReads'] = vcfdf[['RF','RR']].sum(axis=1)
vcfdf['AltReads'] = vcfdf[['AF','AR']].sum(axis=1)

# Also mark variants where count is below read threshold and coverage threshold
vcfdf['MinReadSup'] = np.where(vcfdf[['RefReads','AltReads']].sum(axis=1) >= minreads, 'ok','low')
vcfdf['CoverageThreshold'] = np.where(vcfdf[['RefReads','AltReads']].sum(axis=1) >= covthreshold, 'ok','low')

# If count is below read threshold, blank down to 0
vcfdf['RefReads'] = np.where(vcfdf['RefReads']<= minreads, 0, vcfdf['RefReads'])
vcfdf['AltReads'] = np.where(vcfdf['AltReads']<= minreads, 0, vcfdf['AltReads'])

# get variant percentages
vcfdf['RefPerc'] = (vcfdf['RefReads'] / vcfdf[['RefReads','AltReads']].sum(axis=1))*100
vcfdf['AltPerc'] = (vcfdf['AltReads'] / vcfdf[['RefReads','AltReads']].sum(axis=1))*100

# Work out the minor variant allele percentage
vcfdf['MinorPerc'] = vcfdf[['RefPerc','AltPerc']].min(axis=1)

# Determine variants that fail the bias tests
vcfdf['TestBias'] = 'Fail'
vcfdf.loc[(vcfdf['StrandBias'] >= 0.05) & 
          (vcfdf['BaseQBias'] >=0.05) &
          (vcfdf['MapQBias'] >=0.05) &
          (vcfdf['TailDistBias'] >=0.05), 'TestBias'] = 'Pass'

# Determine Minor SNPs
vcfdf['MinorSNP'] = 'None'
vcfdf.loc[(vcfdf['MinorPerc'] >= minfreq) & 
          (vcfdf['MinReadSup'] =="ok") &
          (vcfdf['INDEL'] == 'FALSE') &
          (vcfdf['PossVariant'] == 'Variant') &
          (vcfdf['StrandBias'] >= 0.05) & 
          (vcfdf['BaseQBias'] >=0.05) &
          (vcfdf['MapQBias'] >=0.05) &
          (vcfdf['TailDistBias'] >=0.05), 'MinorSNP'] = 'Minor'

              
print "Calculating intra-host diversity..."
# Estimate nucleotide diversity (Pi)
# Pi = (D/(D -1)) x (1-(p^2 + (p-1)^2) 
#where 'D' is the sequencing depth at that position and 'p' is the minorvar frequency (not percentage!!)

def getpi(myrow):
    D = vcfdf['RefReads'].loc[myrow] + vcfdf['AltReads'].loc[myrow]
    p = vcfdf['MinorPerc'].loc[myrow]/100
    Pi = ((D/(D-1)) * (1-((p**2) + ((p-1)**2))))
    return Pi
    
vcfdf['Pi'] = getpi(range(len(vcfdf)))


# Estimate Shannon's Entropy (Sn)
# Sn = (-(p * ln(p)) + ((1 - p) * ln(1-p))) / ln(4)
# where p' is the minorvar frequency (not percentage!!)

def getShannon(myrows):
    Snout = []
    for currentrow in myrows:
        p = (vcfdf['MinorPerc'].loc[currentrow])/100
        if (p != 0):
            Sn = (-(p * np.log(p)) + ((1 - p) * np.log(1 - p))) / np.log(4)
            Snout.append(Sn)
        else:
            Sn = 0
            Snout.append(Sn)
    return Snout
    

vcfdf['Sn'] = getShannon(range(len(vcfdf)))


# Get Mean Pi
gw_Pi = np.mean(vcfdf.loc[(vcfdf['CoverageThreshold']=='ok') &
          (vcfdf['MinReadSup']=='ok') &
          (vcfdf['INDEL']=='FALSE') &
          (vcfdf['TestBias']=='Pass') &
          (vcfdf['CHROM']==reference) ]['Pi'])

# and Mean Entropy
gw_Sn = np.mean(vcfdf.loc[(vcfdf['CoverageThreshold']=='ok') &
          (vcfdf['MinReadSup']=='ok') &
          (vcfdf['INDEL']=='FALSE') &
          (vcfdf['TestBias']=='Pass') &
          (vcfdf['CHROM']==reference) ]['Sn'])

#print gw_Pi
#print gw_Sn


#print 'Genome-wide Mean Nucleotide Diversity (Pi; ' + u"\u03c0" + ') = ' + str(round(gw_Pi,6))
print 'Genome-wide Mean Nucleotide Diversity (Pi; ' + str(round(gw_Pi,6))
print 'Genome-wide Mean Shannon\'s Entropy (Sn) = ' + str(round(gw_Sn,6))

with open(vcfinput +'.Stats.txt', "w") as text_file:
#    text_file.write("Pi: %s \nSn: %f" % (gw_Pi, gw_Sn))
    text_file.write("Pi\tSn\n%s\t%f\n" % (gw_Pi, gw_Sn))

print 'Writing mean Pi and Sn to ' + vcfinput + '.Stats.txt'




################################
# include Annotation data and analysis
#	Define methods first:

# Label each annotation "ANN0, ANN1, ANN2..."
def renamecols(nAnns):
    mycolnames = []
    for col in nAnns:
        mycolnames.append('ANN' + str(col))
    return mycolnames

# Deal with 'None' values in empty splits
def ann_none(myrows,mycol):
    myanns = []
    for currentrow in myrows:
        if myANN[mycol].loc[currentrow] is None:
            myanns.append('-')
        else:
            myanns.append(myANN[mycol].loc[currentrow])
    return myanns

# Extract key SNP interpretations (can be used for syn/nonsyn)
def synsnp(myrows,mycol):
    mysnps = []
    for currentrow in myrows:
        if myANN[mycol].loc[currentrow] is None:
            mysnps.append('-')
        elif 'missense_variant' in myANN[mycol].loc[currentrow]:
            mysnps.append('nonsyn')
        elif 'synonymous_variant' in myANN[mycol].loc[currentrow]:
            mysnps.append('syn')
        elif 'stop_lost' in myANN[mycol].loc[currentrow]:
            mysnps.append('stop_lost')
        elif 'stop_gained' in myANN[mycol].loc[currentrow]:
            mysnps.append('stop_gained')
        elif 'start_lost' in myANN[mycol].loc[currentrow]:
            mysnps.append('start_lost')
        elif 'stop_retained_variant' in myANN[mycol].loc[currentrow]:
            mysnps.append('stop_retained')
        elif 'start_retained' in myANN[mycol].loc[currentrow]:    
            mysnps.append('start_retained')
        else:
            mysnps.append('-')
    return mysnps


# Extract key annotations (gene, nucleotide, amino acid change)
def getkeyannot(myrows,mycol):
    myanns = []
    for currentrow in myrows:
        if myANN[mycol].loc[currentrow] is None:
            myanns.append('-')
        elif myANN[mycol].loc[currentrow] == '-':
            myanns.append('-')    
        elif len(myANN[mycol].loc[currentrow]) >=3:
            annot_expand = myANN[mycol].loc[[currentrow]].str.split(pat='|', n=-1, expand=True)
#            gene = annot_expand[3]
	    gene = annot_expand.loc[currentrow,3]
#            nucl = annot_expand[9]
	    nucl = annot_expand.loc[currentrow,9]
#            amino = annot_expand[10]
	    amino = annot_expand.loc[currentrow,10]
            key_annot = str(gene + '|' + nucl + '|' + amino)
            myanns.append(key_annot)
	else:
	    myanns.append('-')
    return myanns


#	Then run each method:
if options.annotation == 'None':
	print 'No annotations - proceeding to next step'
	vcfdf.drop(['PV4','DP4'],axis=1,inplace=True)
else:
# Each variant has multiple annotations (per line) separated by commas
# Split each Annotation up into different columns, creating empty 'None' columns to fill gaps
	print 'Analysing annotations...'
	myANN = vcfdf['ANN'].to_dict()
	myANN = pd.DataFrame(list(myANN.items()), columns = ['position', 'ANN'])
	myANN = myANN['ANN'].str.split(pat=',', n=-1, expand=True)
	
	myANN.columns = renamecols((list(myANN.columns.values)))

# Get list of annotations (i.e. how many annotation columns are present)
        Annot_headers = list(myANN.columns.values)
        print Annot_headers


#	for col in (list(myANN.columns.values)):
        for col in Annot_headers:
	    myANN[col] = ann_none(range(len(myANN)),col)

	myANN = myANN.replace(np.nan, '-', regex=True)

# Determine effect of variant
#	for col in (list(myANN.columns.values)):
        for col in Annot_headers:
	    myANN[str(col + '_syn')] = synsnp(range(len(myANN)),col)

# Get key elements from annotation
	for col in Annot_headers:
	    myANN[str(col + '_var')] = getkeyannot(range(len(myANN)),col)

# Get rid of the original annotation - just keep the pruned bit    
	myANN.drop(Annot_headers,axis=1,inplace=True)

# Merge back into main spreadsheet
	vcfdf.drop(['ANN','PV4','DP4'],axis=1,inplace=True)
	vcfdf = pd.concat([vcfdf[:],myANN[:]], axis=1)
	#vcfdf


#print 'Genome-wide Mean Nucleotide Diversity (Pi; ' + u"\u03c0" + ') = ' + str(round(gw_Pi,6))
#print 'Genome-wide Mean Nucleotide Diversity (Pi; ' + str(round(gw_Pi,6))
#print 'Genome-wide Mean Shannon\'s Entropy (Sn) = ' + str(round(gw_Sn,6))

#with open(vcfinput +'.Stats.txt', "w") as text_file:
#    text_file.write("Pi: %s \nSn: %f" % (gw_Pi, gw_Sn))
#    text_file.write("Pi\tSn\n%s\t%f\n" % (gw_Pi, gw_Sn))	

#print 'Writing mean Pi and Sn to ' + vcfinput + '.Stats.txt'



# Define non-overlapping windows in genome
nwindows = len(vcfdf)/windowsize #determine number of windows in subset


def makewindow(myrows):
    window = []
    for currentrow in myrows:
        window.append((math.trunc(vcfdf['POS'][currentrow] / windowsize) *windowsize) +(windowsize/2))
    return window

vcfdf['window'] = makewindow(range(len(vcfdf)))


# Extract minor variant positions
#minorvar = vcfdf.loc[np.where(vcfdf['MinorSNP']=='Minor')]
minorvar = vcfdf.loc[(vcfdf['MinorSNP']=='Minor') & (vcfdf['CoverageThreshold']=='ok') & (vcfdf['CHROM']==reference)]



# Define function to determine Max minor variant per window
def windowMax(mywindows):
    testwindows=[]
    for currentwindow in mywindows:
        testwindows.append(vcfdf.loc[(vcfdf['window']==currentwindow) &
                                    (vcfdf['CoverageThreshold']=='ok') &
                                    (vcfdf['MinReadSup']=='ok') &
                                    (vcfdf['INDEL']=='FALSE') &
                                    (vcfdf['TestBias']=='Pass') &
                                    (vcfdf['CHROM']==reference) ]['MinorPerc'].max(axis=0))
    return testwindows

# Define function to determine mean Pi per window
def windowPi(mywindows):
    testwindows=[]
    for currentwindow in mywindows:
        testwindows.append(np.mean(vcfdf.loc[(vcfdf['window']==currentwindow) &
                                             (vcfdf['CoverageThreshold']=='ok') &
                                             (vcfdf['MinReadSup']=='ok') &
                                             (vcfdf['INDEL']=='FALSE') &
                                             (vcfdf['TestBias']=='Pass') &
                                             (vcfdf['CHROM']==reference) ]['Pi']))
    return testwindows



# Generate new dataframe with analyses performed per window
if options.graphics == True:
	print "Analysing by "+ str(windowsize) +"sliding windows and generating plots"
	windowed_df = pd.DataFrame({'window':sorted(list(set(vcfdf['window']))),
        	                   'MaxMinor':windowMax(sorted(list(set(vcfdf['window'])))),
                	           'Pi':windowPi(sorted(list(set(vcfdf['window']))))})


# Now try and plot graph
	p_MaxMinor = gg.ggplot(gg.aes('window', 'MaxMinor'),data=windowed_df) +gg.geom_point() +gg.theme_bw() +gg.labs(x="Genome Position (bp; windowsize="+ str(windowsize) +")", y="Minor Variant Frequency (%)") +gg.ggtitle(vcfoutput + "\n Valid Minor Variant Sites :" + str(len(minorvar))) 


# Plot Nucleotide Diversity (Pi) along genome 
	p_pi =gg.ggplot(gg.aes('window', 'Pi'),data=windowed_df) +gg.geom_point() +gg.theme_bw() +gg.labs(x="Genome Position (bp; windowsize="+ str(windowsize) +")", y="Mean nucleotide diversity (" + u"\u03c0" +")") +gg.scale_y_continuous(expand=(0,0),limits=(0, windowed_df['Pi'].max(axis=0)+0.001)) +gg.ggtitle(vcfoutput + "\n Genome-wide Mean Nucleotide Diversity (" +u"\u03c0"+ ") :" +str(round(gw_Pi,6))) 

#p_pi

# Facetted plot (still not sorted y axes labels yet)
	windowed_df_melt = pd.melt(windowed_df, id_vars=['window'])
	p_combi = gg.ggplot(gg.aes('window', 'value',colour='variable'),data=windowed_df_melt)
	p_combi = p_combi + gg.geom_point(colour='variable') + gg.facet_grid('variable',scales='free_y')+gg.theme_bw() +gg.labs(x="Genome Position (bp; windowsize="+ str(windowsize) +")")

# Print graphs to .png
	p_combi.save(vcfinput + ".MinorVar_combo.png")
	p_MaxMinor.save(vcfinput + ".MinorVar.png")
	p_pi.save(vcfinput + ".Pi-diversity.png")



# Print full dataframe and minor vars only to separate tab delimited files
vcfdf.to_csv(vcfinput + ".analysed.tsv",sep='\t', index=False)
minorvar.to_csv(vcfinput + ".minorvars.tsv",sep='\t', index=False)







