#! /usr/bin/env bash

#BSUB -J hw
#BSUB -q normal
#BSUB -o %J.%I.out
#BSUB -e %J.%I.err
#BSUB -R "rusage[mem=50]"
#BSUB -n 12
#BSUB -m compute14

# unzippage
gunzip factorx.chr1.fq.gz
gunzip hg19.chr1.fa.gz

# make some indices after unzipping. bowtie2-build really hates getting
# piped an unzipped file for some reason so don't do that.
bowtie2-build hg19.chr1.fa hg19.chr1

# given the index files already built, run Bowtie alignment for
# dataset, sort it and output a bam file
bowtie2 -x hg19.chr1 -U factorx.chr1.fq | samtools sort > output.bam

# use bedtools to compute factorx coverage along the genome
bedtools genomecov -ibam output.bam -bg > output.bg

# use Macs2 to call peaks
macs2 callpeak -t output.bam -n factorx

# slop some 50bp windows for summit searching
bedtools slop -i factorx_summits.bed -g hg19.chrom.sizes -b 25 > slopped_summits.bed

# sample 1000 peaks so that getting reads doesn't take forever
shuf slopped_summits.bed | head -n 1000 > peaks.rand.1000.bed

# take the narrow peak file and getfasta to get the reads at the peaks
bedtools getfasta -fi hg19.chr1.fa -bed  peaks.rand.1000.bed -fo factorx.fa

# use meme to find the motifs. remember to tell it it's DNA
# this takes forever, specify window size to not get garbage out
# may want to either trim this or start from the summits and awk out a bigger
# window / either way want about a 50bp window
meme factorx.fa -dna -maxsize 10000000 -nsites 5 -maxw 20 -minw 8 -p 12

# get the motif
meme-get-motif -id 1 < meme.txt > out

# put it back into TOMTOM and it should match CTCF.
