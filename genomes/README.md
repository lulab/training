# Genome hg38

Sequence/Annotation/index files

provider: GENCODE

Date: May 4, 2018

Maintainer: Yang Li

## Sequence
### download from gencode v27
#### located in folder: sequence/

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz
```

#### build genome index

```
samtools faidx GRCh38.p10.genome.fa
cut -f1,2 GRCh38.p10.genome.fa.fai > hg38.chrom.sizes
```

## Index 
#### build by bowtie2 and STAR
#### located in folder: index/

using bowtie2

```
mkdir bowtie2_hg38_index
cd bowtie2_hg38_index
bowtie2-build ../sequence/GRCh38.p10.genome.fa GRCh38.p10
```

using STAR

```
mkdir STAR_hg38_index
STAR --runMode genomeGenerate --runThreadN 15 --genomeDir STAR_hg38_index/ --genomeFastaFiles GRCh38.p10.genome.fa
```

## Annotaions
#### see README.md in folder: anno/[README.md](https://github.com/lulab/training/blob/master/genomes/anno/README.md)
