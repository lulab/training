#!/bin/bash

##---------##
## LEVEL I ##
##---------##
# example files are stored in path:
cd /home/younglee/projects/hcc_example/simple_examples

## 1 ##
##  Build bedGraph files

## 1.1 use BEDtools and UCSC Kent Utilities
# sort and convert sam to bam file.
samtools sort NC_1.hg38.sam | samtools view -Sb - > NC_1.hg38.sorted.bam
# build bam index
samtools index NC_1.hg38.sorted.bam 
# create bedGraph
bedtools genomecov -ibam NC_1.hg38.sorted.bam -bga -split -scale 1.0 | sort -k1,1 -k2,2n > NC_1.hg38.sorted.bedGraph
# convert bedGraph to bigWig
bedGraphToBigWig NC_1.hg38.sorted.bedGraph /BioII/lulab_b/shared/genomes/human_hg38/sequence/hg38.chrom.sizes NC_1.hg38.sorted.bw 

## 1.2 use homer 
# create tag directories
makeTagDirectory NC_1.hg38.tagsDir/ NC_1.hg38.sorted.bam
# make bedGraph visualization files
makeUCSCfile NC_1.hg38.tagsDir/ -fragLength given -o auto


## 2 ##
## calculate raw counts/rpkm/rpm
# quatify each gene in Gencode annotations

# 2.1 using HTSeq
htseq-count -m intersection-strict NC_1.hg38.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gff3 > NC_1.hg38.htseq.counts

# 2.2 using featureCounts
featureCounts -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gff3 -o NC_1.hg38.featureCounts.counts NC_1.hg38.sorted.bam

# 2.3 using homer
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gtf hg38 -count exons -d NC_1.hg38.tagsDir/ -noadj > NC_1.hg38.homer.counts
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gtf hg38 -count exons -d NC_1.hg38.tagsDir/ -rpkm > NC_1.hg38.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gtf hg38 -count exons -d NC_1.hg38.tagsDir/ -norm 1e7 > NC_1.hg38.homer.rpm


## 3 ##
## merge to expression matrix (miRNA)
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gtf hg38 \
-d NC_1.hg38.tagsDir/ NC_2.hg38.tagsDir/ NC_3.hg38.tagsDir/ \
BeforeSurgery_1.hg38.tagsDir/ BeforeSurgery_2.hg38.tagsDir/ BeforeSurgery_3.hg38.tagsDir/ \
AfterSurgery_1.hg38.tagsDir/ AfterSurgery_2.hg38.tagsDir/ AfterSurgery_3.hg38.tagsDir/ \
-noadj > $path0/04.counts/hcc_example.hg38.homer.ct.tsv


##----------##
## LEVEL II ## 
##----------##
path0=/home/younglee/projects/hcc_example
## 1 ##
# Build bedGraph files
# use BEDtools and UCSC Kent Utilities
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
samtools sort $path0/02.mapping/$j/hg38/$j.hg38.clean.bam > $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bam
samtools index $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bam
bedtools genomecov -ibam $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bam -bga -split | LC_ALL=C sort -k1,1 -k2,2n > $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bedGraph
bedGraphToBigWig $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bedGraph /BioII/lulab_b/shared/genomes/human_hg38/sequence/hg38.chrom.sizes $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bw
for k in miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon;
do echo $k;
samtools sort $path0/02.mapping/$j/$k/$j.$k.rsem.clean.bam > $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bam
samtools index $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bam
bedtools genomecov -ibam $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bam -bga -split | LC_ALL=C sort -k1,1 -k2,2n > $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bedGraph
bedGraphToBigWig $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bedGraph /BioII/lulab_b/shared/genomes/human_hg38/sequence/hg38.chrom.sizes $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bw
done;
samtools sort $path0/02.mapping/$j/hg38other/$j.hg38other.clean.bam > $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bam
samtools index $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bam
bedtools genomecov -ibam $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bam -bga -split | LC_ALL=C sort -k1,1 -k2,2n > $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bedGraph
bedGraphToBigWig $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bedGraph /BioII/lulab_b/shared/genomes/human_hg38/sequence/hg38.chrom.sizes $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bw
done;

# use homer
mkdir $path0/03.tags/
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
mkdir $path0/03.tags/$j $path0/03.tags/$j/hg38
# Make tag directories for each experiment
makeTagDirectory $path0/03.tags/$j/hg38 $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bam
# Make bedGraph visualization files for each tag directory
makeUCSCfile $path0/03.tags/$j/hg38 -fragLength given -o auto
for k in miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon;
do echo $k;
mkdir $path0/03.tags/$j/$k
makeTagDirectory $path0/03.tags/$j/$k $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bam
makeUCSCfile $path0/03.tags/$j/$k -fragLength given -o auto
done;
mkdir $path0/03.tags/$j/hg38other
makeTagDirectory $path0/03.tags/$j/hg38other $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bam
makeUCSCfile $path0/03.tags/$j/hg38other -fragLength given -o auto
done;


## 2 ##
## calculate raw counts/rpkm/rpm
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
samtools view -h $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bam > $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.sam
samtools view -h $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bam > $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.sam
for k in miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon;
do echo $k;
samtools view -h $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.bam > $path0/02.mapping/$j/$k/$j.$k.rsem.clean.sorted.sam
done;
done;

# merge alignments of RNA types into one single sam/bam file;
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
samtools merge $path0/02.mapping/$j/$j.merged.bam `find $path0/02.mapping/$j -name "*rsem.clean.sorted.bam"` $path0/02.mapping/$j/hg38other/$j.hg38other.clean.sorted.bam
samtools view -h $path0/02.mapping/$j/$j.merged.bam > $path0/02.mapping/$j/$j.merged.sam
makeTagDirectory $path0/03.tags/$j/merged $path0/02.mapping/$j/$j.merged.bam
done;

# count for all the genes in gencode v27
mkdir $path0/04.counts
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
mkdir $path0/04.counts/$j
# HTSeq
htseq-count -m intersection-strict $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gff3 > $path0/04.counts/$j/$j.gencodev27.htseq.ct
# featureCounts
featureCounts -s 1 -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gff3 -o $path0/04.counts/$j/$j.gencodev27.featureCounts.ct $path0/02.mapping/$j/hg38/$j.hg38.clean.sorted.bam
# homer
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gff3 hg38 -count exons -d $path0/03.tags/$j/hg38/ -gid -noadj > $path0/04.counts/$j/$j.gencodev27.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gff3 hg38 -count exons -d $path0/03.tags/$j/hg38/ -gid -rpkm > $path0/04.counts/$j/$j.gencodev27.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/gencode.v27.annotation.gff3 hg38 -count exons -d $path0/03.tags/$j/hg38/ -gid -norm 1e7 > $path0/04.counts/$j/$j.gencodev27.homer.rpm
done;

# count for all miRNAs
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=ID --type=miRNA_primary_transcript $path0/02.mapping/$j/miRNA/$j.miRNA.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff > $path0/04.counts/$j/$j.miRNA.htseq.ct
# featureCounts
featureCounts -t miRNA_primary_transcript -g ID -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff -o $path0/04.counts/$j/$j.miRNA.featureCounts.counts $path0/02.mapping/$j/miRNA/$j.miRNA.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d $path0/03.tags/$j/miRNA/ -gid -noadj > $path0/04.counts/$j/$j.miRNA.homer.ct 
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.miRNA.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.miRNA.homer.rpm
done;


# count for all piRNAs
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=Parent --type=exon $path0/02.mapping/$j/piRNA/$j.piRNA.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/piRNA.gff > $path0/04.counts/$j/$j.piRNA.htseq.ct
# featureCounts
featureCounts -t exon -g Parent -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/piRNA.gff -o $path0/04.counts/$j/$j.piRNA.featureCounts.counts $path0/02.mapping/$j/piRNA/$j.piRNA.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/piRNA.gtf hg38 -count exons -d $path0/03.tags/$j/piRNA/ -gid -noadj > $path0/04.counts/$j/$j.piRNA.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/piRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.piRNA.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/piRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.piRNA.homer.rpm
done;


# count for all Y_RNAs
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=Parent --type=exon $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/Y_RNA_exon.gff > $path0/04.counts/$j/$j.Y_RNA_exon.htseq.ct
# featureCounts
featureCounts -t exon -g Parent -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/Y_RNA_exon.gff -o $path0/04.counts/$j/$j.Y_RNA_exon.featureCounts.counts $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/Y_RNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/Y_RNA_exon/ -gid -noadj > $path0/04.counts/$j/$j.Y_RNA_exon.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/Y_RNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.Y_RNA_exon.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/Y_RNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.Y_RNA_exon.homer.rpm
done;


# count for all snRNAs
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=gene_id --type=snRNA $path0/02.mapping/$j/snRNA/$j.snRNA.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/snRNA.gff > $path0/04.counts/$j/$j.snRNA.htseq.ct
# featureCounts
featureCounts -t snRNA -g gene_id -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/snRNA.gff -o $path0/04.counts/$j/$j.snRNA.featureCounts.counts $path0/02.mapping/$j/snRNA/$j.snRNA.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/snRNA.gtf hg38 -count exons -d $path0/03.tags/$j/snRNA/ -gid -noadj > $path0/04.counts/$j/$j.snRNA.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/snRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.snRNA.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/snRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.snRNA.homer.rpm
done;


# count for all srpRNA_exon
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=Parent --type=exon $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/srpRNA_exon.gff > $path0/04.counts/$j/$j.srpRNA_exon.htseq.ct
# featureCounts
featureCounts -t exon -g Parent -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/srpRNA_exon.gff -o $path0/04.counts/$j/$j.srpRNA_exon.featureCounts.counts $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/srpRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/srpRNA_exon/ -gid -noadj > $path0/04.counts/$j/$j.srpRNA_exon.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/srpRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.srpRNA_exon.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/srpRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.srpRNA_exon.homer.rpm
done;


# count for all tRNAs
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=gene_id --type=tRNA $path0/02.mapping/$j/tRNA/$j.tRNA.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/tRNA.gff > $path0/04.counts/$j/$j.tRNA.htseq.ct
# featureCounts
featureCounts -t tRNA -g gene_id -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/tRNA.gff -o $path0/04.counts/$j/$j.tRNA.featureCounts.counts $path0/02.mapping/$j/tRNA/$j.tRNA.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/tRNA.gtf hg38 -count exons -d $path0/03.tags/$j/tRNA/ -gid -noadj > $path0/04.counts/$j/$j.tRNA.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/tRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.tRNA.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/tRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.tRNA.homer.rpm
done;


# count for all lncRNA_exons
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=gene_id --type=exon $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/lncRNA_exon.gff > $path0/04.counts/$j/$j.lncRNA_exon.htseq.ct
# featureCounts
featureCounts -t exon -g gene_id -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/lncRNA_exon.gff -o $path0/04.counts/$j/$j.lncRNA_exon.featureCounts.counts $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/lncRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/lncRNA_exon/ -gid -noadj > $path0/04.counts/$j/$j.lncRNA_exon.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/lncRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.lncRNA_exon.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/lncRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.lncRNA_exon.homer.rpm
done;


# count for all mRNA_exons
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# HTSeq
htseq-count -m intersection-strict --idattr=gene_id --type=exon $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/mRNA_exon.gff > $path0/04.counts/$j/$j.mRNA_exon.htseq.ct
# featureCounts
featureCounts -t exon -g gene_id -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/mRNA_exon.gff -o $path0/04.counts/$j/$j.mRNA_exon.featureCounts.counts $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.rsem.clean.sorted.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/mRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/mRNA_exon/ -gid -noadj > $path0/04.counts/$j/$j.mRNA_exon.homer.ct
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/mRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -rpkm > $path0/04.counts/$j/$j.mRNA_exon.homer.rpkm
# calculate rpm/cpm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/mRNA_exon.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.mRNA_exon.homer.rpm
done;


## 3 ##
## merge to expression matrix (miRNA)
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 \
-d $path0/03.tags/NC_1/miRNA/ $path0/03.tags/NC_2/miRNA/ $path0/03.tags/NC_3/miRNA/ \
$path0/03.tags/BeforeSurgery_1/miRNA/ $path0/03.tags/BeforeSurgery_2/miRNA/ $path0/03.tags/BeforeSurgery_3/miRNA/ \
$path0/03.tags/AfterSurgery_1/miRNA/ $path0/03.tags/AfterSurgery_2/miRNA/ $path0/03.tags/AfterSurgery_3/miRNA/ \
-gid -noadj > $path0/04.counts/hcc_example.miRNA.homer.ct.tsv

cut -f 1,9- $path0/04.counts/hcc_example.miRNA.homer.ct.tsv > $path0/04.counts/hcc_example.miRNA.homer.ct.mx
# manually modify the header
# geneID    NC_1    NC_2    NC_3    BeforeSurgery_1 BeforeSurgery_2 BeforeSurgery_3 AfterSurgery_1  AfterSurgery_2  AfterSurgery_3
# hsa-mir-1295a   0.000   1.000   0.000   6.000   11.000  5.000   1.000   18.000  3.000
# hsa-mir-1248    6.000   15.000  5.000   49.000  4.000   63.000  118.000 22.000  81.000
# hsa-mir-130b    127.000 118.000 146.000 148.000 336.000 167.000 92.000  129.000 352.000
# hsa-mir-3142    1.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
# hsa-mir-541     2.000   0.000   4.000   1.000   21.000  0.000   0.000   0.000   2.000

for k in miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon;
do echo $k
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/$k.gtf hg38 \
-d $path0/03.tags/NC_1/$k/ $path0/03.tags/NC_2/$k/ $path0/03.tags/NC_3/$k/ \
$path0/03.tags/BeforeSurgery_1/$k/ $path0/03.tags/BeforeSurgery_2/$k/ $path0/03.tags/BeforeSurgery_3/$k/ \
$path0/03.tags/AfterSurgery_1/$k/ $path0/03.tags/AfterSurgery_2/$k/ $path0/03.tags/AfterSurgery_3/$k/ \
-gid -noadj > $path0/04.counts/hcc_example.$k.homer.ct.tsv

analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/$k.gtf hg38 \
-d $path0/03.tags/NC_1/merged/ $path0/03.tags/NC_2/merged/ $path0/03.tags/NC_3/merged/ \
$path0/03.tags/BeforeSurgery_1/merged/ $path0/03.tags/BeforeSurgery_2/merged/ $path0/03.tags/BeforeSurgery_3/merged/ \
$path0/03.tags/AfterSurgery_1/merged/ $path0/03.tags/AfterSurgery_2/merged/ $path0/03.tags/AfterSurgery_3/merged/ \
-gid -rpkm > $path0/04.counts/hcc_example.$k.homer.rpkm.tsv

analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/$k.gtf hg38 \
-d $path0/03.tags/NC_1/merged/ $path0/03.tags/NC_2/merged/ $path0/03.tags/NC_3/merged/ \
$path0/03.tags/BeforeSurgery_1/merged/ $path0/03.tags/BeforeSurgery_2/merged/ $path0/03.tags/BeforeSurgery_3/merged/ \
$path0/03.tags/AfterSurgery_1/merged/ $path0/03.tags/AfterSurgery_2/merged/ $path0/03.tags/AfterSurgery_3/merged/ \
-gid -norm 1e7 > $path0/04.counts/hcc_example.$k.homer.rpm.tsv
done;


## 4 ##
## filter samples and genes
cd $path0/04.counts/
Rscript filter_exp_mx.R






