# Genome hg38

Sequence/Annotation/index files

provider: GENCODE

Date: May 4, 2018

Maintainer: Yang Li

## Sequence
### download from gencode v27
#### located in folder: sequence/

```
wget wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.p10.genome.fa.gz
```

#### parse genome

```
samtools faidx GRCh38.p10.genome.fa
faidx GRCh38.p10.genome.fa -i chromsizes > hg38.chrom.sizes
```

## Index 
#### build by bowtie2 and STAR
#### located in folder: index



## Annotaions
#### see README.md in folder: anno/
