# Genome and Annotations

Path: `/BioII/lulab_b/shared/genomes/hg38/`

## Annotation summary table

| Type | Number of genes  | Source |
| ---- | ---- | ------ |
| miRNA | 1881 | GENCODE V27 (precursor) |
| piRNA | 23431 | piRNABank | 
| lncRNA | 15778 | GENCODE V27 and mitranscriptome |
| rRNA | 544 | GENCODE V27 |
| mRNA | 19836 | GENCODE V27 |
| snoRNA | 943 | GENCODE V27 |
| snRNA | 1900 | GENCODE V27 |
| srpRNA | 680 | GENCODE V27 |
| tRNA | 649 | GENCODE V27 |
| tucpRNA | 3734 | GENCODE V27 |
| Y_RNA | 756 | GENCODE V27 |
| circRNA | 140527 | circBase |
| repeats | - | UCSC Genome Browser (rmsk) |
| promoter | - | ChromHMM tracks from 9 cell lines from UCSC Genome Browser |
| enhancer | - | ChromHMM tracks from 9 cell lines from UCSC Genome Browser |


## Genome and annotation files

| File | Description |
| ---- | ------ |
| `fasta/genome.fa` | genome sequence |
| `fasta/circbase.junction.fa` | junction sequence in circBase |
| `gtf_by_biotype/${rna_type}.gtf` | separate GTF files for each RNA type |
| `gtf/gencode.gtf` | GENCODE GTF file |
| `gtf/mitranscriptome.gtf` | Mitranscriptome GTF file |
| `gtf/long_RNA.gtf` | GTF file of Long RNA (GENCODE + Mitranscriptome - miRNA) |
| `gtf/piRNABank.gtf` | piRNA GTF file from piRNABank |
| `gtf/gencode_tRNA.gtf` | GTF file of tRNA from GENCODE |
| `transcript_table/all.txt` | Table of transcript information (gene_id, transcript_id) |
| `rsem_index/bowtie2/${rna_type}` | RSEM index files for each RNA type |
| `rsem_index/bowtie2/${rna_type}.transcripts.fa` | Sequence for each RNA type |
| `gtf_longest_transcript/${rna_type}.gtf` | GTF files for the longest isoforms from GENCODE and Mitranscriptome |
| `bed/*.bed` | Transcript in BED12 format extracted from GTF files in `gtf/*.gtf |
| `index/bowtie2/circRNA` | Bowtie2 index for cirRNA |
| `long_index/star/` | STAR index including splicing junctions of long RNA |


## Generate the genome and annotation files

### Create genome directory
```bash
[ -d "genome/hg38/source" ] || mkdir -p "genome/hg38/source"
```

### Download chain files for CrossMap
```bash
wget -O genome/hg38/source/hg18ToHg38.over.chain.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
wget -O genome/hg38/source/NCBI36_to_GRCh38.chain.gz https://sourceforge.net/projects/crossmap/files/Ensembl_chain_files/homo_sapiens%28human%29/NCBI36_to_GRCh38.chain.gz
```
### Genome assembly (UCSC hg38)
```bash
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d -c genome/hg38/source/hg38.fa.gz > genome/hg38/fasta/genome.fa
samtools faidx genome/hg38/fasta/genome.fa
```


### ENCODE annotations
```bash
wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
#wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gff3.gz
wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gtf.gz
#wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gff3.gz
wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.tRNAs.gtf.gz
#wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.tRNAs.gff3.gz
zcat genome/hg38/source/gencode.v27.annotation.gtf.gz > genome/hg38/gtf/gencode.gtf
zcat genome/hg38/source/gencode.v27.long_noncoding_RNAs.gtf.gz > genome/hg38/gtf/gencode_lncRNA.gtf
zcat genome/hg38/source/gencode.v27.tRNAs.gtf.gz \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' > genome/hg38/gtf/gencode_tRNA.gtf
# Chain file for converting hg19 to hg38
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

### Mitranscriptome
```bash
wget -P genome/hg38/source http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz
tar -C genome/hg38/source --strip-components=1 -zxf genome/hg38/source/mitranscriptome.gtf.tar.gz mitranscriptome.gtf/mitranscriptome.v2.gtf.gz
# convert from hg19 to hg38
zcat genome/hg38/source/mitranscriptome.v2.gtf.gz \
    | CrossMap.py gff genome/hg38/source/hg19ToHg38.over.chain.gz /dev/stdin genome/hg38/source/mitranscriptome.v2.hg38.gtf
# remove invalid transcripts
bin/preprocess.py fix_gtf -i genome/hg38/source/mitranscriptome.v2.hg38.gtf -o genome/hg38/gtf/mitranscriptome.gtf
```

Extract lncRNA and TUCP RNA to separate GTF files:
```bash
grep 'tcat "lncrna"' genome/hg38/gtf/mitranscriptome.gtf > genome/hg38/gtf/mitranscriptome_lncRNA.gtf
# add exon feature
grep 'tcat "tucp"' genome/hg38/gtf/mitranscriptome.gtf \
    | awk 'BEGIN{OFS="\t";FS="\t"}{print;print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' > genome/hg38/gtf/mitranscriptome_tucp.gtf
cp genome/hg38/gtf/mitranscriptome_tucp.gtf genome/hg38/gtf_by_biotype/tucpRNA.gtf
```

### NONCODE
```bash
wget -P genome/hg38/source http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz
zcat genome/hg38/source/NONCODEv5_human_hg38_lncRNA.gtf.gz \
    | awk 'BEGIN{FS="\t";OFS="\t"}$7 != "." {print $1,"NONCODE",$3,$4,$5,$6,$7,$8,$9}' > genome/hg38/gtf/noncode.gtf
```

### lncRNAs identified in HCC (Nature communications 2017)
```bash
wget -P genome/hg38/source https://media.nature.com/original/nature-assets/ncomms/2017/170213/ncomms14421/extref/ncomms14421-s3.txt
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,"ncomms2017",$3,$4,$5,$6,$7,$8,$9}' genome/hg38/source/ncomms14421-s3.txt > genome/hg38/source/ncomms2017.gtf
CrossMap.py gff genome/hg38/source/hg19ToHg38.over.chain.gz genome/hg38/source/ncomms2017.gtf genome/hg38/source/ncomms2017.hg38.gtf
ln genome/hg38/source/ncomms2017.hg38.gtf genome/hg38/gtf/ncomms2017.gtf
```


### Merge lncRNA (GENCODE and Mitranscriptome)
```bash
cat genome/hg38/gtf/gencode_lncRNA.gtf \
    genome/hg38/gtf/mitranscriptome_lncRNA.gtf \
    > genome/hg38/gtf/merged_lncRNA.gtf
cp genome/hg38/gtf/merged_lncRNA.gtf genome/hg38/gtf_by_biotype/lncRNA.gtf
```

### piRBase (v1.0)
```bash
wget -O genome/hg38/source/piRBase-hsa-v1.0.bed.gz http://www.regulatoryrna.org/database/piRNA/download/archive/v1.0/bed/piR_hg19_v1.0.bed.gz
zcat genome/hg38/source/piRBase-hsa-v1.0.bed.gz \
    | CrossMap.py bed genome/hg38/source/hg19ToHg38.over.chain.gz /dev/stdin genome/hg38/source/piRBase-hsa-v1.0.hg38.bed
bedToGenePred genome/hg38/source/piRBase-hsa-v1.0.hg38.bed genome/hg38/source/piRBase-hsa-v1.0.hg38.genePred
genePredToGtf -source=piRBase file genome/hg38/source/piRBase-hsa-v1.0.hg38.genePred genome/hg38/source/piRBase-hsa-v1.0.hg38.gtf
ln genome/hg38/source/piRBase-hsa-v1.0.hg38.gtf genome/hg38/gtf/piRBase.gtf
```

### piRBase (v2.0)
```bash
wget -O genome/hg38/source/piRBase-hsa-v2.0.bed.gz http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/bed/hsa.bed.gz
zcat genome/hg38/source/piRBase-hsa-v2.0.bed.gz | bedtools sort > source/piRBase-hsa-v2.0.bed
bedToGenePred source/piRBase-hsa-v2.0.bed source/piRBase-hsa-v2.0.genePred
genePredToGtf -source=piRBase file source/piRBase-hsa-v2.0.genePred source/piRBase-hsa-v2.0.gtf
```

### Long RNA (GENCODE + Mitranscriptome - miRNA)
```bash
cat genome/hg38/gtf/gencode.gtf \
    genome/hg38/gtf/mitranscriptome_lncRNA.gtf \
    | grep -v 'gene_type "miRNA' \
    > genome/hg38/gtf/long_RNA.gtf
```

### piRNABank (NCBI36)
```bash
wget -O genome/hg38/source/ http://pirnabank.ibab.ac.in/downloads/all/human_all.zip
unzip genome/hg38/source/human_all.zip -d genome/hg38/source/
mv genome/hg38/source/human_pir.txt genome/hg38/source/piRNABank.human.txt
# Extract genomic coordinates from piRNABank
awk 'BEGIN{OFS="\t"}
    /^>/{na=split(substr($0,2),a,"|");split(a[na],b,":"); 
    if(b[5]=="Plus"){s="+"} else{s="-"}
    if(a[1]!=name){print "chr" b[2],b[3]-1,b[4],a[1],0,s}
    name=a[1]}' genome/hg38/source/piRNABank.human.txt \
    | bedtools sort > genome/hg38/source/piRNABank.human.bed
awk 'BEGIN{OFS="\t"}
    {if($0 ~ /^>/) {split(substr($0,2),a,"|"); 
        if((a[1] != name)&&(length(seq) > 0)){print ">" name;gsub(/U/,"T",seq);print seq} name=a[1]}
    else{seq=$0}}' genome/hg38/source/piRNABank.human.txt > genome/hg38/source/piRNABank.human.fa
bedToGenePred genome/hg38/source/piRNABank.human.bed genome/hg38/source/piRNABank.human.genePred
genePredToGtf -source=piRNABank file genome/hg38/source/piRNABank.human.genePred stdout \
    | awk '$3=="exon"' > genome/hg38/source/piRNABank.human.gtf

CrossMap.py gff genome/hg38/source/hg18ToHg38.over.chain.gz genome/hg38/source/piRNABank.human.gtf \
    genome/hg38/source/piRNABank.human.hg38.gtf
cp genome/hg38/source/piRNABank.human.hg38.gtf genome/hg38/gtf/piRNABank.gtf
cp genome/hg38/gtf/piRNABank.gtf genome/hg38/gtf_by_biotype/piRNA.gtf
gffread --bed -o genome/hg38/source/piRNABank.human.hg38.bed genome/hg38/source/piRNABank.human.hg38.gtf
bedtools getfasta -s -name -fi genome/hg38/fasta/genome.fa -bed genome/hg38/source/piRNABank.human.hg38.bed -split \
    > genome/hg38/source/piRNABank.human.hg38.fa
```

### miRBase
```bash
wget -O genome/hg38/source/miRBase.hsa.gff3 ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
```

### Intron
```bash
bin/preprocess.py extract_gene -i genome/hg38/gtf/long_RNA.gtf | bedtools sort > genome/hg38/bed/long_RNA.gene.bed
awk 'BEGIN{OFS="\t";FS="\t"} !/^#/{match($9,/gene_id "([^"]+)"/,a);print $1,$4-1,$5,a[1],0,$7}' genome/hg38/gtf/long_RNA.gtf \
    | bedtools sort > genome/hg38/bed/long_RNA.exon.bed
bedtools subtract -sorted -s -a genome/hg38/bed/long_RNA.gene.bed -b genome/hg38/bed/long_RNA.exon.bed \
    | bedtools sort > genome/hg38/bed/long_RNA.intron.bed
```

### Promoter/enhancer from ChromHMM (hg19)
```bash
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmH1hescHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHepg2HMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHmecHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHsmmHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHuvecHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhekHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed.g`

# hg19 => hg38
tracks="wgEncodeBroadHmmGm12878HMM wgEncodeBroadHmmH1hescHMM wgEncodeBroadHmmHepg2HMM
    wgEncodeBroadHmmHmecHMM wgEncodeBroadHmmHsmmHMM wgEncodeBroadHmmHuvecHMM
    wgEncodeBroadHmmK562HMM wgEncodeBroadHmmNhekHMM wgEncodeBroadHmmNhlfHMM"
for track in $tracks;do
    CrossMap.py bed genome/hg38/source/hg18ToHg38.over.chain.gz <(zcat genome/hg38/source/${track}.bed.gz) genome/hg38/source/${track}.hg38.bed
    awk 'BEGIN{OFS="\t";FS="\t"}($4=="1_Active_Promoter")||($4=="2_Weak_Promoter")||($4=="3_Poised_Promoter"){print $1,$2,$3,$4,$5,$6}' \
        genome/hg38/source/${track}.hg38.bed | bedtools sort > genome/hg38/bed/promoter.${track}.bed
    awk 'BEGIN{OFS="\t";FS="\t"}($4=="4_Strong_Enhancer")||($4=="5_Strong_Enhancer")||($4=="6_Weak_Enhancer")||($4=="7_Weak_Enhancer"){print $1,$2,$3,$4,$5,$6}' \
        genome/hg38/source/${track}.hg38.bed | bedtools sort > genome/hg38/bed/enhancer.${track}.bed
done
# merge promoters and enhancers from 9 cell lines
cat $(for track in $tracks;do echo genome/hg38/bed/promoter.${track}.bed;done) \
    | bedtools sort | bedtools merge -d 1 \
    | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,"promoter",0,"."}' > genome/hg38/bed/promoter.merged.bed
cat $(for track in $tracks;do echo genome/hg38/bed/enhancer.${track}.bed;done) \
    | bedtools sort | bedtools merge -d 1 \
    | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,"enhancer",0,"."}' > genome/hg38/bed/enhancer.merged.bed
```

### Repeats
UCSC GenomeBrowser -> Tools -> Table Browser

* assembly: GRCh38/hg38
* group: repeats
* track: RepeatMasker
* table: rmsk

Dowload to: genome/hg38/source/rmsk.bed.gz
```bash
gunzip -c genome/hg38/source/rmsk.bed.gz > genome/hg38/bed/rmsk.bed
```

### circRNA database (circBase)
```bash
wget -O genome/hg38/source/circbase.hg19.fa.gz http://www.circbase.org/download/human_hg19_circRNAs_putative_spliced_sequence.fa.gz
zcat genome/hg38/source/circbase.hg19.fa.gz | bin/preprocess.py extract_circrna_junction -s 50 -o genome/hg38/fasta/circbase.junction.fa
samtools faidx genome/hg38/fasta/circbase.junction.fa
```


