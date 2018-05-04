# Annotation files for genome hg38
provider: GENCODE and other projects

format: gtf and gff

Date: May 4, 2018

Maintainer: Yang Li

local: cnode:/BioII/lulab_b/shared/genomes/human_hg38/anno


## statistics

| **RNA_type** | **gene_num** | **transcrips_num** | **source** | **file** |
| :------- |:-------|:------|:-----|:---------|
| rRNA  | 544 | 544 | Gencode27 | rRNA.gencode27.gtf / rRNA.gencode27.gff |
| miRNA | 1,881 | 1,881 | Gencode27 | miRNA.gencode27.gtf / miRNA.gencode27.gff |
| piRNA | 812,347 | 812,347 | piRBase | piRNA.piRBase.hg38.gtf / piRNA.piRBase.hg38.gff |
| snoRNA | 943 | 955 | Gencode27(misc_RNA) | snoRNA.gencode27.gtf / snoRNA.gencode27.gtf |
| snRNA | 1,900 | 1,900 | Gencode27 | snRNA.gencode27.gtf / snRNA.gencode27.gtf |
| srpRNA | 680 | 682 | Gencode27(misc_RNA) | srpRNA.gencode27.gtf / srpRNA.gencode27.gff |
| tRNA | 649 | 649 | Gencode27(predicted tRNA) | tRNA.gencode27.gtf / tRNA.gencode27.gff |
| lncRNA | 15,778 | 27,908 | Gencode27(lincRNA) | lncRNA.gencode27.gtf / lncRNA.gencode27.gff |
| lncRNA | 96,308 | 172,216 | NONCODEv5 | lncRNA.NONCODEv5.hg38.gtf / lncRNA.NONCODEv5.hg38.gff |
| lncRNA | 90,624 | 377,402 | mitranscritome | lncRNA.mitranscriptome.v2.hg38.gtf / lncRNA.mitranscriptome.v2.hg38.gff |
| lncRNA | 136464 | 541,901 | Gencode27+NONCODEv5+ MiTranscriptome+NC2017 | merged_lncRNA.combined.gtf / merged_lncRNA.combined.gff |
| mRNA | 19836 | 80,930 | Gencode27(protein_coding) | mRNA.gencode27.gtf / mRNA.gencode27.gff |
| allGenes | 58,288 | 200,401 | Gencode27 | gencode.v27.annotation.gtf / gencode.v27.annotation.gff |

## pre-process annotaion
### download gencode v27 annotations

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.tRNAs.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.tRNAs.gff3.gz
```

#### parse annotations

```
gtf=gencode.v27.annotation.gtf
for i in rRNA snRNA snoRNA srpRNA miRNA vaultRNA Y_RNA; do
grep "gene_type \"$i\"" $gtf > $i.gencode27.gtf
gffread -E $i.gencode27.gtf -o- > $i.gencode27.gff
done;

grep "gene_type \"protein_coding\"" $gtf | grep -v Selenocysteine > mRNA.gencode27.gtf
grep "RN7SL" $gtf > srpRNA.gencode27.gtf
grep 'Y_RNA' $gtf > Y_RNA.gencode27.gtf

for i in mRNA srpRNA Y_RNA; do
gffread -E $i.gencode27.gtf -o- > $i.gencode27.gff
done;

mv gencode.v27.tRNAs.gtf tRNA.gencode27.gtf
mv gencode.v27.tRNAs.gff3 tRNA.gencode27.gff

mv gencode.v27.long_noncoding_RNAs.gtf lncRNA.gencode27.gtf
mv gencode.v27.long_noncoding_RNAs.gff3 lncRNA.gencode27.gff

mv gencode.v27.annotation.gff3 gencode.v27.annotation.gff
```

### download piRNA from piRBase

### downlaod lncRNA from NONCODE(http://www.noncode.org/), mitranscriptome(http://mitranscriptome.org/), and Yang Yang's NC paper

```
wget http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz
wget http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz
wget https://media.nature.com/original/nature-assets/ncomms/2017/170213/ncomms14421/extref/ncomms14421-s3.txt
```

#### parse and convert
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# liftOver hg19 to hg38
gtfToGenePred -genePredExt mitranscriptome.v2.gtf mitranscriptome.v2.gp
liftOver -genePred mitranscriptome.v2.gp hg19ToHg38.over.chain mitranscriptome.v2.hg38.gp unmapped.gtf
genePredToGtf -utr -source=mitranscriptome file mitranscriptome.v2.hg38.gp lncRNA.mitranscriptome.v2.hg38.gtf
gffread -E lncRNA.mitranscriptome.v2.hg38.gtf -o- > lncRNA.mitranscriptome.v2.hg38.gff

gtfToGenePred -genePredExt ncomms14421-s3.txt lncRNA.lulab_ncomms14421.gp
liftOver -genePred lncRNA.lulab_ncomms14421.gp hg19ToHg38.over.chain lncRNA.lulab_ncomms14421.hg38.gp unmapped.gtf
genePredToGtf -utr -source=lulab_ncomms14421 file lncRNA.lulab_ncomms14421.hg38.gp lncRNA.lulab_ncomms14421.hg38.gtf
gffread -E lncRNA.lulab_ncomms14421.hg38.gtf -o- > lncRNA.lulab_ncomms14421.hg38.gff

rm -rf *.gp unmapped.gtf

mv NONCODEv5_human_hg38_lncRNA.gtf lncRNA.NONCODEv5.hg38.gtf
```

### merge lncRNA annotations

```
gffcompare -o merged_lncRNA -s ../sequence/GRCh38.p12.genome.fa lncRNA.NONCODEv5.hg38.gtf  lncRNA.mitranscriptome.v2.hg38.gtf  lncRNA.gencode27.gtf lncRNA.lulab_ncomms14421.hg38.gtf
gffread -E merged_lncRNA.combined.gtf -o- > merged_lncRNA.combined.gff
```

### cut lncRNA into bins
#### convert combined gtf to combined bed
```
gffread --bed merged_lncRNA.combined.gtf -o merged_lncRNA.combined.bed
```
#### grep the longest transciprt for each lncRNA gene and calculate its length
```
awk 'NR==FNR {split($10,m,"\"");split($12,n,"\"");D[m[2]]=n[2]} \
NR>FNR {split($11,ex,",");for (k in ex) len[$4]+=ex[k]; b[D[$4]] = len[$4] > a[D[$4]] ? $4 : b[D[$4]]; a[D[$4]] = len[$4] > a[D[$4]] ? len[$4] : a[D[$4]] } \
END { OFS = "\t"; for(x in a) print x, a[x],b[x]}' merged_lncRNA.combined.gtf merged_lncRNA.combined.bed >merged_lncRNA.combined.longest.tran
```
#### grep the exons of the longest transciprts
```
awk '{print $3}' merged_lncRNA.combined.longest.tran| grep -Ff - merged_lncRNA.combined.gtf|awk '$3=="exon"' >merged_lncRNA.combined.longest.exon.gtf
```
#### cut full length lncRNA transcripts into bins: bin_size_30, step_size_15
#### (1)deal with transcirpts longer than 30bp
```
awk 'BEGIN{ OFS="\t" } ($5-$4)>30 {split($10,a,"\""); \
for(i=0;i<=($5-$4-30)/15;i++) {print $1,$4+i*15, $4+i*15+30,a[2]"__"$4+i*15"__"$4+i*15+30,$6,$7} \
print $1, $4+(i+1)*15,$5,a[2]"__"$4+(i+1)*15"__"$5,$6,$7}' \
merged_lncRNA.combined.longest.exon.gtf |awk '$3>$2' >merged_lncRNA.combined.longest.exon.bin30.bed
```

#### (2)deal with transcirpts no longer than 30bp
```
awk 'BEGIN{ OFS="\t" } ($5-$4)<=30 {split($10,a,"\"");print $1,$4,$5,a[2]"__"$4"__"$5,$6,$7}' \
merged_lncRNA.combined.longest.exon.gtf >shorterthan30.bed
```
#### (3)concatenate two bed files generated abovely
```
cat merged_lncRNA.combined.longest.exon.bin30.bed shorterthan30.bed|sort -k1,1 -k2,2n -k3,3 > merged_lncRNA.combined.longest.exon.bin30.all.bed
```
