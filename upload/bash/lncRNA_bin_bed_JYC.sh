#! bin/bash

```
This is a script written for selecting lncRNAs from annotation file and converting the coordinates into 30bp bins

Input file: 
gencode.v27.annotation.gtf
biotype.list, which includes all the RNA types you think belonging to lncRNA category
```
#grep lncRNAs from the whole v27.annotation.gtf
awk 'NR==FNR {a[$1]=$1} NR>FNR {for (x=1;x<=NF;x++) if ($x~"^transcript_type"&&($(x+1) in a)) print $0}' biotype.list gencode.v27.annotation.gtf \
 >v27.lncRNA.gtf

#grep the longest transciprts for each lncRNA gene and calculate its length
awk '$3=="transcript" {b[$10]=($5-$4)>a[$10]?$12:b[$10]; a[$10]= ($5-$4)>a[$10]?($5-$4):a[$10]} END { OFS = "\t"; for(x in a) print x, a[x],b[x]}' \
v27.lncRNA.gtf >v27.lncRNA.longest

#grep the exons of the longest transciprts
awk '{print $3}' v27.lncRNA.longest| grep -Ff - v27.lncRNA.gtf |awk '$3=="exon"' >v27.lncRNA.longest.exon.gtf

#convert gtf file to bed file
gffread v27.lncRNA.longest.exon.gtf >v27.lncRNA.longest.exon.bed

#replace gene ID in bed file with gene name;
#extend upstream and downstream coordinates by 200bp seperately
awk 'BEGIN {OFS="\t"} NR==FNR {split($12,a,"\"");G[a[2]]=$16} NR>FNR&&G[$4] {print $1,$2-200,$3+200,G[$4],".",$6}' \
v27.lncRNA.longest.exon.gtf v27.lncRNA.longest.exon.bed >v27.lncRNA.longest.exon.ext.bed

#cut full length lncRNA into bins: bin_size_30, step_size_15
# (1)deal with transcirpts longer than30bp
awk 'BEGIN{ OFS="\t" } ($3-$2)>30 {split($4,a,"\""); \
for(i=0;i<=($3-$2-30)/15;i++) {print $1,$2+i*15, $2+i*15+30,a[2]"__"$2+i*15"__"$2+i*15+30,$6,$7} \
print $1, $2+(i+1)*15,$3,a[2]"__"$2+(i+1)*15"__"$3,$6,$7}' \
v27.lncRNA.longest.exon.ext.bed |awk '$3>$2' >v27.lncRNA.longest.exon.ext.bin30.bed

# (2)deal with transcirpts no longer than 30bp
awk 'BEGIN{ OFS="\t" } ($3-$2)<=30 {split($4,a,"\"");print $1,$2,$3,a[2]"__"$2"__"$3,$6,$7}' \
v27.lncRNA.longest.exon.ext.bed >v27.shoterthan30.bed

# (3)concatenate two bed files
cat v27.lncRNA.longest.exon.ext.bin30.bed|sort -k1,1 -k3,3 -k2,2n >v27.lncRNA.longest.exon.bin30.all.bed








