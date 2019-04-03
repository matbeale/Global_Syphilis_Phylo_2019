#!/bin/bash

echo "Takes two 'kraken.labels' files and two PE fastqs and extracts sequences based on a search term, e.g. 'Homo_sapiens'\n"
if [ $# -ne 7 ]
then
 echo "Usage kraken_bin_fastqs.sh <sampleid> <resultsdir> <searchterm> <read1.fastq> <read2.fastq> <kraken.read1> <kraken.read2>\n"
 exit 1
fi

sampleid=$1
resultsdir=$2
searchterm=$3
read1=$4
read2=$5
kraken1=$6
kraken2=$7

mkdir -pv $resultsdir/

# Extract names of reads with match from each kraken search into single list
grep "$searchterm" $kraken1 | perl -pe 's/\t.+$//g' | perl -pe 's/\/[1,2]{1}//g' >> $resultsdir/$sampleid\.list.tmp
grep "$searchterm" $kraken2 | perl -pe 's/\t.+$//g' | perl -pe 's/\/[1,2]{1}//g' >> $resultsdir/$sampleid\.list.tmp

# To get all paired hits, need to use intersect of both R1 and R2, so get uniques and relabel
cat $resultsdir/$sampleid\.list.tmp | sort | uniq > $resultsdir/$sampleid\.unique-list.tmp
cat $resultsdir/$sampleid\.unique-list.tmp | perl -pe 's/$/\/1/g' | perl -pe 's/^\/1//g' > $resultsdir/$sampleid\.$searchterm\.R1.list
cat $resultsdir/$sampleid\.unique-list.tmp | perl -pe 's/$/\/2/g' | perl -pe 's/^\/2//g' > $resultsdir/$sampleid\.$searchterm\.R2.list
rm $resultsdir/$sampleid\.list.tmp
rm $resultsdir/$sampleid\.unique-list.tmp

# Now use Seqtk to extract the reads from each fastq
seqtk subseq $read1 $resultsdir/$sampleid\.$searchterm\.R1.list | gzip > $resultsdir/$sampleid\.$searchterm\.binned_1.fastq.gz
seqtk subseq $read2 $resultsdir/$sampleid\.$searchterm\.R2.list | gzip > $resultsdir/$sampleid\.$searchterm\.binned_2.fastq.gz



