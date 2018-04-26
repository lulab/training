#!/bin/bash

PATH1=/BioII/lulab_b/shared/projects/exRNA/hcc_examples
PATH2=/home/zhangyanjun/projects/exRNA/hcc_examples
PATH3=/BioII/lulab_b/shared/genomes/human_hg38
PATH4=/home/zhangyanjun/genome/human_hg38

#Link fastq files to home directory
cd $PATH2
mkdir 00.rawdata
cd 00.rawdata
ln -s $PATH1/fastq/* ./

#create sample_name.txt to record names of fastq files
ls *fastq | while read id; do echo ${id%.*} >> sample_name; done

#create mapping_order.txt to record the order of mapping
echo -e "00.rRNA\n01.miRNA\n02.piRNA\n03.Y_RNA\n04.snRNA\n05.srpRNA\n06.tRNA\n07.lncRNA\n08.mRNA\n09.other_human_genome" > ../mapping_order

#build index
#link genomic annotation file to home directory
cd 
mkdir -p genome/human_hg38/index/bowtie2_index
mkdir $PATH4/gtf
mkdir $PATH4/sequence
ln -s $PATH3/gtf/* $PATH4/gtf/
ln -s $PATH3/sequence/GRCh38.p12.genome.fa $PATH4/sequence/
ln -s $PATH3/index/bowtie2_index/bowtie2_hg38_index/ $PATH4/index/bowtie2_index/

#get RNA sequence(.fa) from annotation file(.gtf)
cd $PATH4/index/bowtie2_index/
for i in `cat $PATH2/mapping_order`;do mkdir ${i%.*}.bowtie2_${i#*.}_index;done
cd ./00.bowtie2_rRNA_index
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/rRNA_exon.gff -fo rRNA.fa >getfasta.log 2>getfasta.err.log
cd ../01.bowtie2_miRNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/miRNA.gff -fo miRNA.fa >getfasta.log 2>getfasta.err.log
cd ../02.bowtie2_piRNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/piR_hg38_exon.gtf -fo piRNA.fa >getfasta.log 2>getfasta.err.log
cd ../03.bowtie2_Y_RNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/Y_RNA_exon.gff -fo Y_RNA.fa >getfasta.log 2>getfasta.err.log
cd ../04.bowtie2_snRNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/snRNA.gff -fo snRNA.fa >getfasta.log 2>getfasta.err.log
cd ../05.bowtie2_srpRNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/srpRNA_exon.gff -fo srpRNA.fa >getfasta.log 2>getfasta.err.log
cd ../06.bowtie2_tRNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/tRNA.gff -fo tRNA.fa >getfasta.log 2>getfasta.err.log
cd ../07.bowtie2_lncRNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/lncRNA_exon.gff -fo lncRNA.fa >getfasta.log 2>getfasta.err.log
cd ../08.bowtie2_mRNA_index/
bedtools getfasta -fi ../../../sequence/GRCh38.p12.genome.fa -bed ../../../gtf/mRNA_exon.gff -fo mRNA.fa >getfasta.log 2>getfasta.err.log
cd ../09.bowtie2_other_human_genome_index/
ln -s /BioII/lulab_b/shared/genomes/human_hg38/index/bowtie2_hg38_index/* ./

#build bowtie2 index
for i in `head -8 $PATH2/mapping_order`;do cd $PATH4/index/bowtie2_index/${i%.*}.bowtie2_${i#*.}_index; bowtie2-build ${i#*.}.fa ${i#*.} >${i#*.}.log 2>${i#*.}.err.log;done

#fastqc
cd $PATH2
mkdir 01.fastqc
cd 01.fastqc
for i in `cat $PATH2/00.rawdata/sample_name`;do fastqc ../00.rawdata/${i}.fastq -o ./ >./$i.fastqc.log 2>./$i.fastqc.err.log; done

#multiqc
#multiqc *zip >multiqc.log 2>multiqc.err.log

#remove adapters and trim reads
cd ../
mkdir 02.cutadapt
cd 02.cutadapt
for i in `cat ../00.rawdata/sample_name`;do mkdir $i;done
for i in `cat ../00.rawdata/sample_name`;do cutadapt -u -100 -q 30,30 --trim-n -m 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --too-short-output=./$i/$i.trimmed.cutAdapt3.TooShort.fastq -o ./$i/$i.trimmed.cutAdapt3.fastq ../00.rawdata/$i.fastq > ./$i/$i.cutAdapt3.log 2> ./$i/$i.cutAdapt3.err.log;done

#2nd fastQC
cd ../
mkdir 03.fastqc
for i in `cat ./00.rawdata/sample_name`; do fastqc ./02.cutadapt/$i/$i.trimmed.cutAdapt3.fastq -o ./03.fastqc > ./03.fastqc/fastqc2.log 2> ./03.fastqc/fastqc2.err.log;done

#multiqc
#multiqc *zip >multiqc.log 2>multiqc.err.log

#remove ribosomal RNA mapped reads
mkdir 04.rRNA
cd 04.rRNA
for i in `cat ../00.rawdata/sample_name`; do mkdir $i; done
for i in `cat ../00.rawdata/sample_name`; do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/00.bowtie2_rRNA_index/rRNA -S ./$i/$i.rRNA.sam --un ./$i/$i.no_rRNA.fastq ../02.cutadapt/$i/$i.trimmed.cutAdapt3.fastq > ./$i/$i.rRNA.log 2> ./$i/$i.rRNA.err.log;done

#sam to bam
#for i in `cat ../00.rawdata/sample_name.txt`;do samtools view -bS ./$i/Si.rRNA.sam > $i.rRNA.bam

#mapping reads to human different RNA types.
cd ../
mkdir 05.mapping
cd 05.mapping
mv ../mapping_order ./
for i in `cat ../sample_name`;do mkdir $i;done

for i in `cat ../sample_name`;do mkdir $i/00.hg38;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/09.bowtie2_other_human_genome_index/GRCh38_p10 -S ./$i/00.hg38/$i.hg38.sam --un ./$i/00.hg38/$i.hg38.unAligned.fq ../04.rRNA/$i/$i.no_rRNA.fastq > ../log/$i.hg38.log 2> ../log/$i.hg38.err.log;done

for i in `cat ../sample_name`;do mkdir $i/01.miRNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/01.bowtie2_miRNA_index/miRNA -S ./$i/01.miRNA/$i.miRNA.sam --un ./$i/01.miRNA/$i.miRNA.unAligned.fq ../04.rRNA/$i/$i.no_rRNA.fastq > ../log/$i.miRNA.log 2> ../log/$i.miRNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/02.piRNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/02.bowtie2_piRNA_index/piRNA -S ./$i/02.piRNA/$i.piRNA.sam --un ./$i/02.piRNA/$i.piRNA.unAligned.fq ./$i/01.miRNA/$i.miRNA.unAligned.fq > ../log/$i.piRNA.log 2> ../log/$i.piRNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/03.Y_RNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/03.bowtie2_Y_RNA_index/Y_RNA -S ./$i/03.Y_RNA/$i.Y_RNA.sam --un ./$i/03.Y_RNA/$i.Y_RNA.unAligned.fq ./$i/02.piRNA/$i.piRNA.unAligned.fq > ../log/$i.Y_RNA.log 2> ../log/$i.Y_RNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/04.snRNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/04.bowtie2_snRNA_index/snRNA -S ./$i/04.snRNA/$i.snRNA.sam --un ./$i/04.snRNA/$i.snRNA.unAligned.fq ./$i/03.Y_RNA/$i.Y_RNA.unAligned.fq > ../log/$i.snRNA.log 2> ../log/$i.snRNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/05.srpRNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/05.bowtie2_srpRNA_index/srpRNA -S ./$i/05.srpRNA/$i.srpRNA.sam --un ./$i/05.srpRNA/$i.srpRNA.unAligned.fq ./$i/04.snRNA/$i.snRNA.unAligned.fq > ../log/$i.srpRNA.log 2> ../log/$i.srpRNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/06.tRNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/06.bowtie2_tRNA_index/tRNA -S ./$i/06.tRNA/$i.tRNA.sam --un ./$i/06.tRNA/$i.tRNA.unAligned.fq ./$i/05.srpRNA/$i.srpRNA.unAligned.fq > ../log/$i.tRNA.log 2> ../log/$i.tRNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/07.lncRNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/07.bowtie2_lncRNA_index/lncRNA -S ./$i/07.lncRNA/$i.lncRNA.sam --un ./$i/07.lncRNA/$i.lncRNA.unAligned.fq ./$i/06.tRNA/$i.tRNA.unAligned.fq > ../log/$i.lncRNA.log 2> ../log/$i.lncRNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/08.mRNA;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index/08.bowtie2_mRNA_index/mRNA -S ./$i/08.mRNA/$i.mRNA.sam --un ./$i/08.mRNA/$i.mRNA.unAligned.fq ./$i/07.lncRNA/$i.lncRNA.unAligned.fq > ../log/$i.mRNA.log 2> ../log/$i.mRNA.err.log;done

for i in `cat ../sample_name`;do mkdir $i/09.other_human_genome;done
for i in `cat ../sample_name`;do bowtie2 --sensitive-local -x $PATH4/index/bowtie2_index//09.bowtie2_other_human_genome_index/GRCh38_p10 -S ./$i/09.other_human_genome/$i.other_human_genome.sam --un ./$i/09.other_human_genome/$i.other_human_genome.unAligned.fq ./$i/08.mRNA/$i.mRNA.unAligned.fq > ../log/$i.other_human_genome.log 2> ../$i.other_human_genome.err.log;done
