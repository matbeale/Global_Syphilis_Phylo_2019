#!/bin/bash

# For converting SRA/ENA fastq headers to locally compatible - removes characters and whitespace between read header and fwd/rev identifier (/1 /2)
# "Usage: fix-fastqs-headers.sh <fastq-prefix>"


if [ $# -ne 1 ]
then
 echo "For converting SRA/ENA fastq headers to locally compatible - removes characters and whitespace between read header and fwd/rev identifier (/1 /2)"
 echo "Usage: fix-fastqs-headers.sh <fastq-prefix>"
 exit 1
fi


#filelist=$1
line=$1


#for i in {1..2} ; do while read line ; do cp $line\_$i\.fastq.gz $line\.f_$i\.fastq.gz ; zcat $line\.f_$i\.fastq.gz | perl -pwe 'if (/^@/) { s/\.[0-9]+\K\ .*\//\// }' | gzip -c > $line\_$i\.fastq.gz ; rm $line\.f_$i\.fastq.gz ; done < $filelist ; done

for i in {1..2} ; do cp $line\_$i\.fastq.gz $line\.f_$i\.fastq.gz ; zcat $line\.f_$i\.fastq.gz | perl -pwe "if (/^@/) { s/\.[0-9]+\K\ .*\//\// }" | gzip -c > $line\_$i\.fastq.gz ; rm $line\.f_$i\.fastq.gz  ; done


