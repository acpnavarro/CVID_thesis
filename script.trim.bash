#!/bin/sh
adaptor=$(readlink -f /media/anna/4TB_2/Anna_Navarro/adaptor.fa) 
mkdir OUT

forw=($ls *_R1.fastq.gz) 
rev=($ls *_R2.fastq.gz) 
na=($(ls *_R1.fastq.gz | cut -c 1-9)) 
let i=0
for sample in *_R2.fastq.gz
do
echo $i 
main=$(echo "${forw[i]}") 
mate=$(echo "${rev[i]}") 
names1=$(echo "${na[i]}") 

trimmomatic PE $main $mate -baseout OUT/${names1}.fq ILLUMINACLIP:$adaptor:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25
let i=$i+1
done






