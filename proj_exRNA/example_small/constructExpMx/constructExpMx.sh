#!/bin/bash

##---------##
## LEVEL I ##
##---------##
# example files are stored in path:
cd /home/younglee/projects/bioinfoTraining/constructExpMx/

## 1 ##
##  Build bedGraph files

## 1.1 use BEDtools and UCSC Kent Utilities
# sort and convert sam to bam file.
samtools sort NC_1.miRNA.rsem.bam > NC_1.miRNA.sorted.bam
# build bam index
samtools index NC_1.miRNA.sorted.bam 
# create bedGraph
bedtools genomecov -ibam NC_1.miRNA.sorted.bam -bga -split -scale 1.0 | sort -k1,1 -k2,2n > NC_1.miRNA.sorted.bedGraph
# convert bedGraph to bigWig
bedGraphToBigWig NC_1.miRNA.sorted.bedGraph /BioII/lulab_b/shared/genomes/human_hg38/sequence/hg38.chrom.sizes NC_1.miRNA.sorted.bw 

## 1.2 use homer 
# create tag directories
makeTagDirectory NC_1.miRNA.tagsDir/ NC_1.miRNA.sorted.bam
# make bedGraph visualization files
makeUCSCfile NC_1.miRNA.tagsDir/ -fragLength given -o auto


## 2 ##
## calculate raw counts/rpkm/rpm
# quatify each gene in Gencode annotations

# 2.1 using HTSeq
samtools view NC_1.miRNA.rsem.bam > NC_1.miRNA.rsem.sam 
htseq-count -m intersection-strict --idattr=Name --type=miRNA_primary_transcript NC_1.miRNA.rsem.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff > NC_1.miRNA.htseq.counts

# 2.2 using featureCounts
featureCounts -t miRNA_primary_transcript -g Name -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff -o NC_1.miRNA.featureCounts.counts NC_1.miRNA.sorted.bam

# 2.3 using homer
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d NC_1.miRNA.tagsDir/ -noadj > NC_1.miRNA.homer.counts
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d NC_1.miRNA.tagsDir/ -rpkm > NC_1.miRNA.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d NC_1.miRNA.tagsDir/ -norm 1e7 > NC_1.miRNA.homer.rpm


## 3 ##
## merge to expression matrix (miRNA)
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 \
-d NC_1.miRNA.tagsDir/ NC_2.miRNA.tagsDir/ NC_3.miRNA.tagsDir/ \
BeforeSurgery_1.miRNA.tagsDir/ BeforeSurgery_2.miRNA.tagsDir/ BeforeSurgery_3.miRNA.tagsDir/ \
AfterSurgery_1.miRNA.tagsDir/ AfterSurgery_2.miRNA.tagsDir/ AfterSurgery_3.miRNA.tagsDir/ \
-noadj > hcc_example.miRNA.homer.ct.tsv

cut -f 1,9- $path0/04.counts/hcc_example.homer.ct.tsv $path0/04.counts/hcc_example.homer.ct.mx
# manually modify the header
# geneID    NC_1    NC_2    NC_3    BeforeSurgery_1 BeforeSurgery_2 BeforeSurgery_3 AfterSurgery_1  AfterSurgery_2  AfterSurgery_3
# hsa-mir-1295a   0.000   1.000   0.000   6.000   11.000  5.000   1.000   18.000  3.000
# hsa-mir-1248    6.000   15.000  5.000   49.000  4.000   63.000  118.000 22.000  81.000
# hsa-mir-130b    127.000 118.000 146.000 148.000 336.000 167.000 92.000  129.000 352.000
# hsa-mir-3142    1.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
# hsa-mir-541     2.000   0.000   4.000   1.000   21.000  0.000   0.000   0.000   2.000


## 4 ##
## filter samples and genes
Rscript bin/filter_exp_mx.R






