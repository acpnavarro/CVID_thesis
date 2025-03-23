
#!/bin/sh

mkdir Index

STAR --runMode genomeGenerate --genomeDir Index --genomeFastaFiles Ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
--sjdbGTFfile Ref/Homo_sapiens.GRCh38.110.gtf --runThreadN 6 --sjdbOverhang 75



wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz







#convert do bed
convert2bed -i gtf < input.gtf > output.bed
