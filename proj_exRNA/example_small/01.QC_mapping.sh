#!/bin/bash
path0=/home/younglee/projects/hcc_example

## 00 ##
## fetch sequence based on gtf, build index for each RNA type.
for i in rRNA_exon miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon;
do echo $i;
gffread -w $path0/src/$i.fa -g /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p12.genome.fa /BioII/lulab_b/shared/genomes/human_hg38/gtf/$i.gff
# bedtools getfasta -fi /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p12.genome.fa -bed /BioII/lulab_b/shared/genomes/human_hg38/gtf/$i.gff -fo $path0/src/$i.fa
# bowtie2-build $path0/src/$i.fa $path0/src/$i > $path0/log/$i.index.log
rsem-prepare-reference --gtf /BioII/lulab_b/shared/genomes/human_hg38/gtf/$i.gtf --bowtie2 /BioII/lulab_b/shared/genomes/human_hg38/sequence/GRCh38.p12.genome.fa $path0/src/${i}
done;

## 01 ##
## fastqc
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do echo $i;
fastqc -o $path0/00.rawdata/ $path0/00.rawdata/$i;
done

## 02 ##
## cut out low-quality nucleotide and remove adapters
mkdir $path0/01.preprocess
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
# fastx_trimmer -v -Q 33 -f 6 -l 75 -i $path0/00.rawdata/$j.fastq -o $path0/00.rawdata/$j.trim.fastq >> $path0/log/$j.trim.cutadapt.log 2>> $path0/log/$j.trim.cutadapt.err
cutadapt -u -100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 16 --trim-n --too-short-output=$path0/01.preprocess/$j.tooShort.fastq -o $path0/01.preprocess/$j.cutadapt.fastq $path0/00.rawdata/$j.fastq >> $path0/log/$j.cutadapt.log 2>> $path0/log/$j.cutadapt.err;
fastqc -o $path0/01.preprocess/ $path0/01.preprocess/$j.cutadapt.fastq
done;

## 03 ##
## remove rRNA mapped reads
mkdir $path0/02.mapping
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
mkdir $path0/02.mapping/$j $path0/02.mapping/$j/no_rRNA
bowtie2 -p 4 --sensitive-local --no-unal --un $path0/02.mapping/$j/no_rRNA/$j.no_rRNA.fq -x $path0/src/rRNA_exon $path0/01.preprocess/$j.cutadapt.fastq -S $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.sam > $path0/log/$j.rm_rRNA.log 2> $path0/log/$j.rm_rRNA.err
samtools view -b -f 16 $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.sam > $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.sam > $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.clean.bam
samtools fastq $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.reverseMap.bam > $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.reverseMap.fastq
cat $path0/02.mapping/$j/no_rRNA/$j.no_rRNA.fq $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.reverseMap.fastq > $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.unmapped.fastq
done;

## mapping to whole human genome
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
mkdir $path0/02.mapping/$j/hg38
bowtie2 -p 4 --sensitive-local --no-unal --un $path0/02.mapping/$j/hg38/$j.unAligned.fq -x /BioII/lulab_b/shared/genomes/human_hg38/index/bowtie2_hg38_index/GRCh38_p10 $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.unmapped.fastq -S $path0/02.mapping/$j/hg38/$j.hg38.sam > $path0/log/$j.map2hg38.log 2> $path0/log/$j.map2hg38.err
samtools view -S -b  $path0/02.mapping/$j/hg38/$j.hg38.sam > $path0/02.mapping/$j/hg38/$j.hg38.clean.bam 
ln -s $path0/02.mapping/$j/hg38/$j.unAligned.fq $path0/02.mapping/$j/hg38/$j.hg38.unmapped.fastq
done;

## sequential map to each RNA type 
## rRNA_exon miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# map to miRNA
mkdir $path0/02.mapping/$j/miRNA
echo "start map to miRNA" >> $path0/log/$j.miRNA.log 2>> $path0/log/$j.miRNA.err
bowtie2 -p 4 --sensitive-local --no-unal --un $path0/02.mapping/$j/miRNA/$j.unAligned.fq -x $path0/src/miRNA $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.unmapped.fastq -S $path0/02.mapping/$j/miRNA/$j.miRNA.sam >> $path0/log/$j.miRNA.log 2>> $path0/log/$j.miRNA.err
samtools view -b -f 16 $path0/02.mapping/$j/miRNA/$j.miRNA.sam > $path0/02.mapping/$j/miRNA/$j.miRNA.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/miRNA/$j.miRNA.sam > $path0/02.mapping/$j/miRNA/$j.miRNA.clean.bam
rsem-tbam2gbam $path0/src/miRNA $path0/02.mapping/$j/miRNA/$j.miRNA.clean.bam $path0/02.mapping/$j/miRNA/$j.miRNA.rsem.clean.bam > $path0/log/$j.miRNA.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/miRNA/$j.miRNA.reverseMap.bam > $path0/02.mapping/$j/miRNA/$j.miRNA.reverseMap.fastq
cat $path0/02.mapping/$j/miRNA/$j.unAligned.fq $path0/02.mapping/$j/miRNA/$j.miRNA.reverseMap.fastq > $path0/02.mapping/$j/miRNA/$j.miRNA.unmapped.fastq

# map to piRNA
mkdir $path0/02.mapping/$j/piRNA
echo "start map to piRNA" >> $path0/log/$j.piRNA.log 2>> $path0/log/$j.piRNA.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/piRNA/$j.unAligned.fq -x $path0/src/piRNA $path0/02.mapping/$j/miRNA/$j.miRNA.unmapped.fastq -S $path0/02.mapping/$j/piRNA/$j.piRNA.sam >> $path0/log/$j.piRNA.log 2>> $path0/log/$j.piRNA.err
samtools view -b -f 16 $path0/02.mapping/$j/piRNA/$j.piRNA.sam > $path0/02.mapping/$j/piRNA/$j.piRNA.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/piRNA/$j.piRNA.sam > $path0/02.mapping/$j/piRNA/$j.piRNA.clean.bam
rsem-tbam2gbam $path0/src/piRNA $path0/02.mapping/$j/piRNA/$j.piRNA.clean.bam $path0/02.mapping/$j/piRNA/$j.piRNA.rsem.clean.bam > $path0/log/$j.piRNA.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/piRNA/$j.piRNA.reverseMap.bam > $path0/02.mapping/$j/piRNA/$j.piRNA.reverseMap.fastq
cat $path0/02.mapping/$j/piRNA/$j.unAligned.fq $path0/02.mapping/$j/piRNA/$j.piRNA.reverseMap.fastq > $path0/02.mapping/$j/piRNA/$j.piRNA.unmapped.fastq

# map to Y_RNA
mkdir $path0/02.mapping/$j/Y_RNA_exon
echo "start map to Y_RNA_exon" >> $path0/log/$j.Y_RNA_exon.log 2>> $path0/log/$j.Y_RNA_exon.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/Y_RNA_exon/$j.unAligned.fq -x $path0/src/Y_RNA_exon $path0/02.mapping/$j/piRNA/$j.piRNA.unmapped.fastq -S $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.sam >> $path0/log/$j.Y_RNA_exon.log 2>> $path0/log/$j.Y_RNA_exon.err
samtools view -b -f 16 $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.sam > $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.sam > $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.clean.bam
rsem-tbam2gbam $path0/src/Y_RNA_exon $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.clean.bam $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.rsem.clean.bam > $path0/log/$j.Y_RNA_exon.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.reverseMap.bam > $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.reverseMap.fastq
cat $path0/02.mapping/$j/Y_RNA_exon/$j.unAligned.fq $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.reverseMap.fastq > $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.unmapped.fastq

# map to snRNA
mkdir $path0/02.mapping/$j/snRNA
echo "start map to snRNA" >> $path0/log/$j.snRNA.log 2>> $path0/log/$j.snRNA.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/snRNA/$j.unAligned.fq -x $path0/src/snRNA $path0/02.mapping/$j/Y_RNA_exon/$j.Y_RNA_exon.unmapped.fastq -S $path0/02.mapping/$j/snRNA/$j.snRNA.sam >> $path0/log/$j.snRNA.log 2>> $path0/log/$j.snRNA.err
samtools view -b -f 16 $path0/02.mapping/$j/snRNA/$j.snRNA.sam > $path0/02.mapping/$j/snRNA/$j.snRNA.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/snRNA/$j.snRNA.sam > $path0/02.mapping/$j/snRNA/$j.snRNA.clean.bam
rsem-tbam2gbam $path0/src/snRNA $path0/02.mapping/$j/snRNA/$j.snRNA.clean.bam $path0/02.mapping/$j/snRNA/$j.snRNA.rsem.clean.bam > $path0/log/$j.snRNA.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/snRNA/$j.snRNA.reverseMap.bam > $path0/02.mapping/$j/snRNA/$j.snRNA.reverseMap.fastq
cat $path0/02.mapping/$j/snRNA/$j.unAligned.fq $path0/02.mapping/$j/snRNA/$j.snRNA.reverseMap.fastq > $path0/02.mapping/$j/snRNA/$j.snRNA.unmapped.fastq

# map to srpRNA
mkdir $path0/02.mapping/$j/srpRNA_exon
echo "start map to srpRNA_exon" >> $path0/log/$j.srpRNA_exon.log 2>> $path0/log/$j.srpRNA_exon.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/srpRNA_exon/$j.unAligned.fq -x $path0/src/srpRNA_exon $path0/02.mapping/$j/snRNA/$j.snRNA.unmapped.fastq -S $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.sam >> $path0/log/$j.srpRNA_exon.log 2>> $path0/log/$j.srpRNA_exon.err
samtools view -b -f 16 $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.sam > $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.sam > $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.clean.bam
rsem-tbam2gbam $path0/src/srpRNA_exon $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.clean.bam $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.rsem.clean.bam > $path0/log/$j.srpRNA_exon.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.reverseMap.bam > $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.reverseMap.fastq
cat $path0/02.mapping/$j/srpRNA_exon/$j.unAligned.fq $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.reverseMap.fastq > $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.unmapped.fastq

# map to tRNA
mkdir $path0/02.mapping/$j/tRNA
echo "start map to tRNA" >> $path0/log/$j.tRNA.log 2>> $path0/log/$j.tRNA.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/tRNA/$j.unAligned.fq -x $path0/src/tRNA $path0/02.mapping/$j/srpRNA_exon/$j.srpRNA_exon.unmapped.fastq -S $path0/02.mapping/$j/tRNA/$j.tRNA.sam >> $path0/log/$j.tRNA.log 2>> $path0/log/$j.tRNA.err
samtools view -b -f 16 $path0/02.mapping/$j/tRNA/$j.tRNA.sam > $path0/02.mapping/$j/tRNA/$j.tRNA.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/tRNA/$j.tRNA.sam > $path0/02.mapping/$j/tRNA/$j.tRNA.clean.bam
rsem-tbam2gbam $path0/src/tRNA $path0/02.mapping/$j/tRNA/$j.tRNA.clean.bam $path0/02.mapping/$j/tRNA/$j.tRNA.rsem.clean.bam > $path0/log/$j.tRNA.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/tRNA/$j.tRNA.reverseMap.bam > $path0/02.mapping/$j/tRNA/$j.tRNA.reverseMap.fastq
cat $path0/02.mapping/$j/tRNA/$j.unAligned.fq $path0/02.mapping/$j/tRNA/$j.tRNA.reverseMap.fastq > $path0/02.mapping/$j/tRNA/$j.tRNA.unmapped.fastq

# map to lncRNA
mkdir $path0/02.mapping/$j/lncRNA_exon
echo "start map to lncRNA_exon" >> $path0/log/$j.lncRNA_exon.log 2>> $path0/log/$j.lncRNA_exon.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/lncRNA_exon/$j.unAligned.fq -x $path0/src/lncRNA_exon $path0/02.mapping/$j/tRNA/$j.tRNA.unmapped.fastq -S $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.sam >> $path0/log/$j.lncRNA_exon.log 2>> $path0/log/$j.lncRNA_exon.err
samtools view -b -f 16 $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.sam > $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.sam > $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.clean.bam
rsem-tbam2gbam $path0/src/lncRNA_exon $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.clean.bam $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.rsem.clean.bam > $path0/log/$j.lncRNA_exon.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.reverseMap.bam > $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.reverseMap.fastq
cat $path0/02.mapping/$j/lncRNA_exon/$j.unAligned.fq $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.reverseMap.fastq > $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.unmapped.fastq

# map to mRNA
mkdir $path0/02.mapping/$j/mRNA_exon
echo "start map to mRNA_exon" >> $path0/log/$j.mRNA_exon.log 2>> $path0/log/$j.mRNA_exon.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/mRNA_exon/$j.unAligned.fq -x $path0/src/mRNA_exon $path0/02.mapping/$j/lncRNA_exon/$j.lncRNA_exon.unmapped.fastq -S $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.sam >> $path0/log/$j.mRNA_exon.log 2>> $path0/log/$j.mRNA_exon.err
samtools view -b -f 16 $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.sam > $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.reverseMap.bam
samtools view -b -F 16 $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.sam > $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.clean.bam
rsem-tbam2gbam $path0/src/mRNA_exon $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.clean.bam $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.rsem.clean.bam > $path0/log/$j.mRNA_exon.tbam2gbam.log
samtools fastq $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.reverseMap.bam > $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.reverseMap.fastq
cat $path0/02.mapping/$j/mRNA_exon/$j.unAligned.fq $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.reverseMap.fastq > $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.unmapped.fastq

# map to genome hg38other region
mkdir $path0/02.mapping/$j/hg38other
echo "start map to hg38other" >> $path0/log/$j.hg38other.log 2>> $path0/log/$j.hg38other.err
bowtie2 -p 4 --sensitive-local --un $path0/02.mapping/$j/hg38other/$j.unAligned.fq -x /BioII/lulab_b/shared/genomes/human_hg38/index/bowtie2_hg38_index/GRCh38_p10 $path0/02.mapping/$j/mRNA_exon/$j.mRNA_exon.unmapped.fastq -S $path0/02.mapping/$j/hg38other/$j.hg38other.sam >> $path0/log/$j.hg38other.log 2>> $path0/log/$j.hg38other.err
samtools view -S -b $path0/02.mapping/$j/hg38other/$j.hg38other.sam > $path0/02.mapping/$j/hg38other/$j.hg38other.clean.bam
ln -s $path0/02.mapping/$j/hg38other/$j.unAligned.fq $path0/02.mapping/$j/hg38other/$j.hg38other.unmapped.fastq
done;


## 05 ##
## summary the results and statistics
mkdir $path0/stat
## number of mapped readsfor each RNA types
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
# total library size
libSizeN=`cat $path0/log/$j.cutadapt.log | grep 'Total reads processed' | awk 'BEGIN{FS=OFS=":"}{print $2}' | sed 's/^ *//g' | sed 's/,//g'`
echo -e "$j\tpreprocess\tlibSizeN\t$libSizeN" >> $path0/stat/$j.readsN.stat.tsv
# too short
tooShortN=`cat $path0/log/$j.cutadapt.log | grep 'Reads that were too short' | awk 'BEGIN{FS=OFS=":"}{print $2}' | sed 's/^ *//g' | sed 's/,//g' | cut -d ' ' -f 1`
echo -e "$j\tpreprocess\ttooShortN\t$tooShortN" >> $path0/stat/$j.readsN.stat.tsv
# clean reads
cleanN=`cat $path0/log/$j.cutadapt.log | grep 'Reads written' | awk 'BEGIN{FS=OFS=":"}{print $2}' | sed 's/^ *//g' | sed 's/,//g' | cut -d ' ' -f 1`
echo -e "$j\tpreprocess\tcleanN\t$cleanN" >> $path0/stat/$j.readsN.stat.tsv
# rRNA mapped reads
rRNA_N=`samtools flagstat $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.clean.bam | awk 'NR==5{print $1}'`
echo -e "$j\tpreprocess\trRNA_N\t$rRNA_N" >> $path0/stat/$j.readsN.stat.tsv
# kept reads
keptN=`grep '@' $path0/02.mapping/$j/no_rRNA/$j.rRNA_exon.unmapped.fastq | wc -l`
echo -e "$j\tpreprocess\tkeptN\t$keptN" >> $path0/stat/$j.readsN.stat.tsv
# map to hg38
hg38_N=`samtools flagstat $path0/02.mapping/$j/hg38/$j.hg38.clean.bam | awk 'NR==5{print $1}'`
echo -e "$j\tmap2hg38\thg38\t$hg38_N" >> $path0/stat/$j.readsN.stat.tsv
# map to different RNA types
for k in miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon;
do echo $k;
RNAtypes_N=`samtools flagstat $path0/02.mapping/$j/$k/$j.$k.rsem.clean.bam | awk 'NR==5{print $1}'`
echo -e "$j\tsequentialMap\t$k\t$RNAtypes_N" >> $path0/stat/$j.readsN.stat.tsv
done;
# map to hg38 other region
hg38other_N=`samtools flagstat $path0/02.mapping/$j/hg38other/$j.hg38other.clean.bam | awk 'NR==5{print $1}'`
echo -e "$j\tmap2hg38other\thg38other\t$hg38other_N" >> $path0/stat/$j.readsN.stat.tsv
# non-human
nonHuman_N=`grep '@' $path0/02.mapping/$j/hg38other/$j.hg38other.unmapped.fastq | wc -l`
echo -e "$j\tmap2hg38other\tnonHuman_N\t$nonHuman_N" >> $path0/stat/$j.readsN.stat.tsv
done;


## length of mapped reads for each RNA types
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
for k in miRNA piRNA Y_RNA_exon snRNA srpRNA_exon tRNA lncRNA_exon mRNA_exon;
do echo $k;
samtools view $path0/02.mapping/$j/$k/$j.$k.rsem.clean.bam | awk 'BEGIN{FS=OFS="\t"}{print length($10)}' | sort -n | uniq -c | awk 'BEGIN{FS=" "; OFS="\t"}{print $2,$1}' | sort -nk1,1 | sed -e "s/^/$j\t$k\t/g" >> $path0/stat/$j.lengthN.stat.tsv
done;
sed -i -e "1i sample\ttype\tlen\tnum" $path0/stat/$j.lengthN.stat.tsv
done;

# plot the histogram
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
Rscript $path0/bin/plot_readsLens.R $path0/stat/$j.lengthN.stat.tsv $path0/stat/$j.lengthN.stat.pdf
done;


#---------------------------------------
# try collapse reads after adapter removement
# also consider the min read length.
path0=/home/younglee/projects/hcc_example
mkdir $path0/01.preprocess
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for l in 16 18 20 36;
do echo $l;
cutadapt -u -100 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m $l --trim-n --too-short-output=$path0/01.preprocess/$j.L${l}.tooShort.fastq -o $path0/01.preprocess/$j.cutadapt.L${l}.fastq $path0/00.rawdata/$j.fastq >> $path0/log/$j.cutadapt.L${l}.log 2>> $path0/log/$j.cutadapt.L${l}.err;
fastx_collapser -i $path0/01.preprocess/$j.cutadapt.L${l}.fastq -o $path0/01.preprocess/$j.cutadapt.L${l}.collapsed.fastq
done;
done;

# stat number
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for l in 16 18 20 36;
do echo $l;
readsN=`grep '@' $path0/01.preprocess/$j.cutadapt.L${l}.fastq | wc -l`;
collapsedReadsN=`grep '>' $path0/01.preprocess/$j.cutadapt.L${l}.collapsed.fastq | wc -l`;
echo -e "$j\t$l\t$readsN\t$collapsedReadsN" >> $path0/stat/readsL.collapse.num.sta.txt
done;
done;

# map to miRNA transcriptome
# cheak the mapping ratio of exactly and >1 mapped reads
for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
for l in 16 18 20 36;
do echo $l;
# map to miRNA
mkdir $path0/02.mapping/$j/miRNA
echo "start map to miRNA" >> $path0/log/$j.miRNA.log 2>> $path0/log/$j.L${l}.miRNA.err
bowtie2 -p 4 --sensitive-local --no-unal --un $path0/02.mapping/$j/miRNA/$j.L${l}.unAligned.fq -x $path0/src/miRNA $path0/01.preprocess/$j.cutadapt.L${l}.fastq -S $path0/02.mapping/$j/miRNA/$j.L${l}.miRNA.sam >> $path0/log/$j.L${l}.miRNA.log 2>> $path0/log/$j.L${l}.miRNA.err
done;
done;

for i in `ls $path0/00.rawdata | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
for l in 16 18 20 36;
do echo $l;
exactlyN=`cat $path0/log/$j.L${l}.miRNA.err | grep 'exactly' | awk 'BEGIN{FS=OFS=" "}{print $1}'`
multiN=`cat $path0/log/$j.L${l}.miRNA.err | grep '>1' | awk 'BEGIN{FS=OFS=" "}{print $1}'`
RNAtypes_N=`echo $(($exactlyN + $multiN))`
echo -e "$j\tsequentialMap\t$l\t$exactlyN\t$multiN\t$RNAtypes_N" >> $path0/stat/readsN.readsL.sta.tsv
done;
done;


