#!/bin/bash

PATH0=/Share/home/wangsiqi/projects/01.HCC_biomarker/01.exRNA-seq/01.LuLab_HCC_exRNA

for i in `cat ./file_name.new`
do
        awk '(NR%4==2){print}' $PATH0/01.Pre_cutadapt/${i}/${i}.trimmed.cutAdapt3.fastq | awk '{print length($NF)}' | sort -n | uniq -c >clean_reads_length/${i}.clean_reads_length;

        samtools view -F 4 $PATH0/03.ordered_mapping_bowtie2/${i}/00.hg38/${i}.hg38.sam | awk '{print $10}' | awk '{print length($NF)}' | sort -n | uniq -c > 00.hg38_reads_length/${i}.hg38_reads_length;

        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/01.miRNA/${i}.miRNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 01.miRNA_reads_length/${i}.miRNA_reads_length;
        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/02.piRNA/${i}.piRNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 02.piRNA_reads_length/${i}.piRNA_reads_length;
        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/03.Y_RNA/${i}.Y_RNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 03.Y_RNA_reads_length/${i}.Y_RNA_reads_length;
        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/04.snRNA/${i}.snRNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 04.snRNA_reads_length/${i}.snRNA_reads_length;
        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/05.srpRNA/${i}.srpRNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 05.srpRNA_reads_length/${i}.srpRNA_reads_length;
        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/06.tRNA/${i}.tRNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 06.tRNA_reads_length/${i}.tRNA_reads_length;
        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/07.other_lncRNA/${i}.other_lncRNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 07.other_lncRNA_reads_length/${i}.other_lncRNA_reads_length;
        awk '(NR%4==2){print}' $PATH0/03.ordered_mapping_bowtie2/${i}/08.mRNA/${i}.mRNA.mapped.fastq |awk '{print length($NF)}' | sort -n | uniq -c > 08.mRNA_reads_length/${i}.mRNA_reads_length;
done

        paste  00.hg38_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >00.hg38_reads_length/hg38_reads.length;

        paste  01.miRNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >01.miRNA_reads_length/miRNA_reads.length;

        paste  02.piRNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >02.piRNA_reads_length/piRNA_reads.length;

        paste  03.Y_RNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >03.Y_RNA_reads_length/Y_RNA_reads.length;

        paste  04.snRNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >04.snRNA_reads_length/snRNA_reads.length;

        paste  05.srpRNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >05.srpRNA_reads_length/srpRNA_reads.length;

        paste  06.tRNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >06.tRNA_reads_length/tRNA_reads.length;

        paste  07.other_lncRNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >07.other_lncRNA_reads_length/other_lncRNA_reads.length;

        paste  08.mRNA_reads_length/* | awk '{print $2,$1,$3,$5,$7,$9,$11,$13,$15,$17,$19,$21,$23,$25,$27,$29,$31,$33,$35,$37,$39,$41,$43,$45,$47,$49,$51,$53,$55,$57,$59,$61,$63,$65,$67,$69,$71,$73,$75,$77,$79,$81,$83,$85,$87,$89,$91,$93,$95,$97,$99,$101,$103,$105,$107,$109,$111,$113,$115,$117,$119,$121}' >08.mRNA_reads_length/mRNA_reads.length;
