#!/bin/bash
set -e

#inputs:
# $1 - genome
# $2 - gff3
# $3 - gff3 header to use (TAIR 'gene_id=')
# $4 - promoter length
# $5 onwards - arguments for MEME-LAB (also takes $4)

#input argument order (for MEME-LAB itself):
# - promoter sequence file (generated in-script here)
# - promoter length
# - cluster membership file
# - number of motifs to generate per cluster
# - minimum motif width
# - maximum motif width

#start off by filtering the .gff3 to gene lines only
if [ ! -f annot.gff3 ]
then
	cp $2 annot.gff3
fi
grep -P '\tgene\t' annot.gff3 > genelines.gff3

#strip the potential FASTA line breaks. creates genome_stripped.fa
if [ ! -f genome.fa ]
then
	cp $1 genome.fa
fi
python3 /scripts/strip_newlines.py

#create the .genome file
samtools faidx genome_stripped.fa
cut -f 1-2 genome_stripped.fa.fai > bedgenome.genome

#parse up the .bed for promoter extraction
python3 /scripts/parse_genelines.py $3
#the python script takes the genelines.gff3 file and makes a genelines.bed out of it
#prepare promoter region information
bedtools flank -l $4 -r 0 -s -i genelines.bed -g bedgenome.genome > promoters.bed
bedtools getfasta -fi genome_stripped.fa -bed promoters.bed -s -fo promoters_rough.fa -name
#no longer needed - -name does this
#python3 /scripts/parse_promoters.py

#meme_lab is picky and demands a .fasta file specifically for some reason
mv promoters.fa promoters.fasta

#okay, now we can run meme-lab
perl /scripts/web_service_meme_lab.pl promoters.fasta "${@:4}"

#cleanup
rm bedgenome.genome
rm genelines.bed
rm genelines.gff3
rm genome.fa
rm genome_stripped.fa
rm genome_stripped.fa.fai
rm promoters.bed
rm promoters.fasta
rm promoters_rough.fa