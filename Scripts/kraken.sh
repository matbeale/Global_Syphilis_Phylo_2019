#!/bin/bash


# usage for the command
if [ $# -ne 4 ]
then
 echo "Usage: kraken.sh <runid> <input.fastq.gz> </path/to/results> <number-of-threads>"
 exit 1
fi


runid=$1
inputfile=$2
results_dir=$3
cpus=$4

mkdir -pv $results_dir/

kraken_Db=/lustre/scratch118/infgen/pathogen/pathpipe/kraken/pi_qc_2015521/
kraken=/software/pathogen/external/apps/usr/bin/kraken

# Run Kraken
$kraken --threads $cpus --fastq-input --gzip-compressed --classified-out $results_dir/$runid\.classified.fastq --unclassified-out $results_dir/$runid\.unclassified.fastq --output $results_dir/$runid\.kraken.out --db $kraken_Db $inputfile

#$kraken --threads $cpus --fastq-input --classified-out $results_dir/$runid\.classified.txt --unclassified-out $results_dir/$runid\.unclassified.txt --output $results_dir/$runid\.kraken.out --db $kraken_Db $inputfile


# Label the output with taxonomy data
$kraken\-translate --mpa-format --db $kraken_Db $results_dir/$runid\.kraken.out > $results_dir/$runid\.kraken.labels

# Get a summary report of sample contents
$kraken\-report --db $kraken_Db $results_dir/$runid\.kraken.out > $results_dir/$runid\.kraken.report

# Get a summary report in mpa format (metaphlan style) - can also be run on multiple samples
$kraken\-mpa-report --db $kraken_Db $results_dir/$runid\.kraken.out > $results_dir/$runid\.kraken.mpa-report


# gzip up some of the files (can be big for fastqs)
gzip $results_dir/$runid\.kraken.out
gzip $results_dir/$runid\.unclassified.fastq
gzip $results_dir/$runid\.classified.fastq
# get rid of the different iteratations - takes up tonnes of space
rm -r $results_dir/K*/ 


