######

#!/bin/bash


index=../INDEX

FILES=../FASTQ/OUT/*.gz

for f in $FILES
do
echo $f
base=$(basename $f | cut -c 1-9)
echo $base
STAR --runThreadN 3 --genomeDir $index --readFilesIn ../FASTQ/OUT/${base}_1P.fq.gz ../FASTQ/OUT/${base}_2P.fq.gz --outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix ./$base"_"
done
