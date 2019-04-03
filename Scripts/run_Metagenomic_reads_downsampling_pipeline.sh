#!/bin/bash


# Usage: run_Metagenomic_reads_downsampling_pipeline.sh [options: [-t] [-c]] <sampleid> <xxx_1.fastq.gz> <xxx_2.fastq.gz> <resultsdir>

# Takes raw fastq files and performs read assignment with kraken, then bins the reads according to defined criteria (e.g. 'Treponema'), followed by downsampling to a specific upper count (note that samples with fewer reads will be labelled as downsampled even if the read counts remain the same). Note that the kraken database is ~30GB, but is not complete - the pipeline will not work well for samples with limited sequence information, highly diverse genomes, or for genomes with novel gene content"

# Submits jobs directly to LSF - does not need to be bsubbed.




##########################################################
# Specify default options for flags (this was originally written for Treponema)'
readcount=2500000
binning_term='Treponema'

# Capture command line options and positional arguments
while getopts ":ht:c:" opt; do
  case "${opt}" in
    h) echo ""
       echo "Uses kraken to classify and bin raw PE fastqs from a species group, then trims and downsamples to a manageable number for assembly/mapping."
       echo "Note: Jobs are submitted to LSF automatically, this does not need to be Bsub-ed"
       echo ""
       echo "Usage: "
       echo "    cmd [options: [-t] [-c]] <sampleid> <xxx_1.fastq.gz> <xxx_2.fastq.gz> <resultsdir>" 
       echo "Options: "
       echo "    -t  -  Specify species to bin reads by [Treponema]"
       echo "    -c  -  Specify max number of reads to retain after binning [2500000]"
       echo ""
       exit 0
      ;;
    t) binning_term="${OPTARG}" ;;
    c) readcount="${OPTARG}" ;;
    \?)echo "Invalid Option: -$OPTARG" 1>&2
       exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# All 4 positional arguments are mandatory - kill script if any are absent
if [ $# -ne 4 ]
then
 echo "Invalid number of positional arguments "
 echo "cmd [options: [-t] [-c]] <sampleprefix> <xxx_1.fastq.gz> <xxx_2.fastq.gz> <resultsdir>"
 exit 1
fi

# Print optional flags to screen
echo "Binning reads classified as $binning_term"
echo "Downsampling binned reads to $readcount reads or less"
#########################################################

#  Does not need to be bsub-ed

sampleid=$1
read1=$2
read2=$3
resultsdir=$4

#binning_term=Treponema
#binning_term=Chlamydia
#binning_term=Mycoplasma

#readcount=2500000



# run kraken on reads (submits 2 jobs)
dokrakenR1=$(bsub -J kraken_$sampleid\_1 -o $resultsdir/k_$sampleid\_1.%J.o -e $resultsdir/k_$sampleid\_1.%J.e -q normal -n 8 -M32000 -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]' "/nfs/users/nfs_m/mb29/bsub_scripts/kraken.sh $sampleid\_1 $read1 $resultsdir/kraken/ 8")
newdokrakenR1=$(echo $dokrakenR1 |perl -pe 's/\>.+$//g'|perl -pe 's/^.+\<//g')
echo "Kraken job submitted as $newdokrakenR1 for $read1"

dokrakenR2=$(bsub -J kraken_$sampleid\_2 -o $resultsdir/k_$sampleid\_2.%J.o -e $resultsdir/k_$sampleid\_2.%J.e -q normal -n 8 -M32000 -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]' "/nfs/users/nfs_m/mb29/bsub_scripts/kraken.sh $sampleid\_2 $read1 $resultsdir/kraken/ 8")
newdokrakenR2=$(echo $dokrakenR2 |perl -pe 's/\>.+$//g'|perl -pe 's/^.+\<//g')
echo "Kraken job submitted as $newdokrakenR2 for $read2"




# Combine the kraken binning and downsampling into same job (reduces number of jobs needed) and also include a read trimming step after
# This job will submit with the others, but will remain pending until completion of the kraken runs
#dokbin=$(bsub -w "done($newdokrakenR1)&&done($newdokrakenR2)"

dokbin=$(bsub -w "ended($newdokrakenR1)&&ended($newdokrakenR2)" -J kbin_$sampleid -o $resultsdir/kbin_$sampleid.%J.o -e $resultsdir/kbin_$sampleid.%J.e -q normal -n 2 -M32000 -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]' "/nfs/users/nfs_m/mb29/scripts/kraken_bin_fastqs.sh $sampleid $resultsdir $binning_term $read1 $read2 $resultsdir/kraken/$sampleid\_1.kraken.labels $resultsdir/kraken/$sampleid\_2.kraken.labels ; java -Xmx32g -jar /software/pathogen/external/apps/usr/bin/trimmomatic-0.33.jar PE -phred33 -threads 2 $resultsdir/$sampleid\.$binning_term\.binned_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.binned_2.fastq.gz $resultsdir/$sampleid\.$binning_term\.paired_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.unpaired_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.paired_2.fastq.gz $resultsdir/$sampleid\.$binning_term\.unpaired_2.fastq.gz ILLUMINACLIP:/nfs/users/nfs_m/mb29/references/AdaptorTrimming/illumina-adaptors.2.fasta:2:30:10 LEADING:3 MINLEN:40 ; /nfs/users/nfs_m/mb29/scripts/downsample_fastq.sh $resultsdir/$sampleid $resultsdir/$sampleid\.$binning_term\.paired_1.fastq.gz $resultsdir/$sampleid\.$binning_term\.paired_2.fastq.gz $readcount ; rm $resultsdir/$sampleid\.$binning_term\.unpaired_*.fastq.gz")

newdokbin=$(echo $dokbin |perl -pe 's/\>.+$//g'|perl -pe 's/^.+\<//g')
echo "Kraken read binning and downsampling submitted as $newdokbin - will start after kraken runs finish"



