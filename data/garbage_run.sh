#! /usr/bin/env bash

#BSUB -J hw
#BSUB -q normal
#BSUB -o %J.%I.out
#BSUB -e %J.%I.err
#BSUB -R "rusage[mem=50]"
#BSUB -n 12
#BSUB -m compute14

# SAME
# unzippage
gunzip factorx.chr1.fq.gz
gunzip hg19.chr1.fa.gz

# SAME
# make some indices after unzipping. bowtie2-build really hates getting
# piped an unzipped file for some reason so don't do that.
bowtie2-build hg19.chr1.fa hg19.chr1

# DIFFERENCE: I get my samtools sort output with > instead of -o
# given the index files already built, run Bowtie alignment for
# dataset, sort it and output a bam file
bowtie2 -x hg19.chr1 -U factorx.chr1.fq | samtools sort > output.bam

# DIFFERENCE: I don't give genomecov the chrom.sizes file with -g
# THIS DIFFERENCE ONLY AFFECTS MAKING THE BIGWIG TRACK, SO NOT THE PROBLEM
# use bedtools to compute factorx coverage along the genome
bedtools genomecov -ibam output.bam -bg > output.bg

# DIFFERENCE: other script specifies -f BAM (instead of using auto format parameter)
# DOES THIS MATTER? NOT SURE. YOU'D THINK NOT IF IT HAS AUTO DETECTION, BUT...
# use Macs2 to call peaks
macs2 callpeak -t output.bam -n factorx

# DIFFERENCE: other script slops 50bp in each direction, but I tried that too so NOPE.
# slop some 50bp windows for summit searching
bedtools slop -i factorx_summits.bed -g hg19.chrom.sizes -b 25 > slopped_summits.bed

# SAME
# sample 1000 peaks so that getting reads doesn't take forever
shuf slopped_summits.bed | head -n 1000 > peaks.rand.1000.bed

# SAME
# take the narrow peak file and getfasta to get the reads at the peaks
bedtools getfasta -fi hg19.chr1.fa -bed  peaks.rand.1000.bed -fo factorx.fa

# DIFF: I specify a # of motif sites (5) and have 10x larger maxsize
# DIFF: plus I give Tesla 12 threads to speed this thing up
# use meme to find the motifs. remember to tell it it's DNA
# this takes forever, specify window size to not get garbage out
# may want to either trim this or start from the summits and awk out a bigger
# window / either way want about a 50bp window
meme factorx.fa -dna -maxsize 10000000 -nsites 5 -maxw 20 -minw 8 -p 12

# SAME
# get the motif
meme-get-motif -id 1 < meme.txt > out

# put it back into TOMTOM and it should match CTCF.
# BUT IT NEVER DOES
