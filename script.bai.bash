#!/bin/sh

for f in *.bam
do 
samtools index $f
done


#convert do bed
convert2bed -i gtf < input.gtf > output.bed
