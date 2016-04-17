#! /usr/bin/env bash

#BSUB -J hw
#BSUB -q normal
#BSUB -o %J.%I.out
#BSUB -e %J.%I.err
#BSUB -R "rusage[mem=50]"
#BSUB -m compute14

# SAME
# unzippage
# gunzip factorx.chr1.fq.gz
# gunzip hg19.chr1.fa.gz

# SAME
# make some indices after unzipping. bowtie2-build really hates getting
# piped an unzipped file for some reason so don't do that.
bowtie2-build hg19.chr1.fa hg19.chr1

# CHANGED TO SAME 
# given the index files already built, run Bowtie alignment for
# dataset, sort it and output a bam file
bowtie2 -x hg19.chr1 -U factorx.chr1.fq | samtools sort -o output.bam

# CHANGED: to specify genome file 
# even though again, this can't be the problem b/c it doesn't get used for peakcalls
# use bedtools to compute factorx coverage along the genome
bedtools genomecov -ibam output.bam -bg -g hg19.chrom.sizes > output.bg

# CHANGED to specify -f BAM (instead of using auto format parameter)
# use Macs2 to call peaks
macs2 callpeak -t output.bam -f BAM -n factorx

# CHANGED to slop 50bp in each direction like the other script even though come on
# slop some 50bp windows for summit searching
bedtools slop -i factorx_summits.bed -g hg19.chrom.sizes -b 50 > slopped_summits.bed

# SAME
# sample 1000 peaks so that getting reads doesn't take forever
shuf slopped_summits.bed | head -n 1000 > peaks.rand.1000.bed

# SAME
# take the narrow peak file and getfasta to get the reads at the peaks
bedtools getfasta -fi hg19.chr1.fa -bed  peaks.rand.1000.bed -fo factorx.fa

# CHANGED: maxsize (down 10x) and removed nsites and removed  multithreading
# use meme to find the motifs. remember to tell it it's DNA
# this takes forever, specify window size to not get garbage out
# may want to either trim this or start from the summits and awk out a bigger
# window / either way want about a 50bp window
meme factorx.fa -dna -maxsize 1000000 -maxw 20 -minw 8

# SAME
# get the motif
meme-get-motif -id 1 < meme.txt > out

# put it back into TOMTOM and it should match CTCF.
# maybe the 17th time is the charm??
