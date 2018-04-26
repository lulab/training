#THis is a bash script to summarize RNA-seq mapping ratio step by step

#!/bin/bash

PATH0=/home/zhangyanjun/projects/exRNA/hcc_examples

cd $PATH0
mkdir ./data_summary
touch ./data_summary/data_summary.header
echo -e "Sample\tRaw_reads\tClean_reads\tClean_reads%\trRNA\trRNA%\tKept_reads\thg38\thg38%\tmiRNA\tmiRNA%\tpiRNA\tpiRNA%\tY_RNA\tY_RNA%\tsnRNA\tsnRNA%\tsrpRNA\tsrpRNA%\ttRNA\ttRNA%\tlncRNA\tlncRNA%\tmRNA\tmRNA%\tOther_human_region\tOther_human_region%\thg38_ordered\thg38_ordered%\tNonhuman\tNonhuman%" > ./data_summary/data_summary.header

for i in `cat ./sample_name`
do
    cat ./02.cutadapt/$i/$i.cutAdapt3.log | grep 'Total reads processed' | awk 'BEGIN{FS=OFS=":"}{print $2}' | sed 's/^[ \t]*//g' | sed 's/,//g' >> ./data_summary/01.Raw_reads;
    cat ./02.cutadapt/$i/$i.cutAdapt3.log | grep 'Reads that were too short' | awk 'BEGIN{FS=OFS=":"}{print $2}' | sed 's/^[ \t]*//g' | sed 's/,//g' | sed 's/ //g' >> ./data_summary/02.reads_shorter_than_15nt;
    cat ./02.cutadapt/$i/$i.cutAdapt3.log | grep 'Reads written' | awk 'BEGIN{FS=OFS=":"}{print $2}' | sed 's/^[ \t]*//g' | sed 's/,//g' | sed 's/ //g' | sed 's/(.*)//g' >> ./data_summary/03.Clean_reads_n;
    cat ./02.cutadapt/$i/$i.cutAdapt3.log | grep 'Reads written' | awk 'BEGIN{FS=OFS=":"}{print $2}' | sed 's/^[ \t]*//g' | sed 's/,//g' | sed 's/ //g' | sed 's/^.*(://g' | sed 's/).*$//g' >> ./data_summary/03.Clean_reads_ratio;
    paste ./data_summary/03.Clean_reads_n ./data_summary/03.Clean_reads_ratio > ./data_summary/03.Clean_reads;
    wc -l ./04.rRNA/$i/$i.no_rRNA.fastq |awk 'reads=(($1/4)){printf("%d\n",reads)}' >> ./data_summary/05.Kept_reads;
    
    tail -6  ./04.rRNA/$i/$i.rRNA.err.log | grep 'exactly' | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/04.rRNA_reads_1;
    tail -6  ./04.rRNA/$i/$i.rRNA.err.log | grep '>1' |  awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/04.rRNA_reads_n;
    paste -d "\t" ./data_summary/04.rRNA_reads_1 ./data_summary/04.rRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/04.rRNA_reads;
    paste ./data_summary/04.rRNA_reads ./data_summary/03.Clean_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/04.rRNA_reads - > ./data_summary/04.rRNA;
    
    tail -6 ./log/$i.hg38.err.log | grep 'exactly' | head -1 | awk 'BEGIN{F S=OFS=" "}{print $1}' >> ./data_summary/06.hg38_reads_1;
    tail -6 ./log/$i.hg38.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/06.hg38_reads_n;
    paste -d "\t" ./data_summary/06.hg38_reads_1 ./data_summary/06.hg38_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/06.hg38_reads;
    paste ./data_summary/06.hg38_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/06.hg38_reads - > ./data_summary/06.hg38;

    tail -6 ./log/$i.miRNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/11.miRNA_reads_1;
    tail -6 ./log/$i.miRNA.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/11.miRNA_reads_n;
    paste -d "\t" ./data_summary/11.miRNA_reads_1 ./data_summary/11.miRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/11.miRNA_reads;
    paste ./data_summary/11.miRNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/11.miRNA_reads - > ./data_summary/11.miRNA;
    
    tail -6 ./log/$i.piRNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/12.piRNA_reads_1;
    tail -6 ./log/$i.piRNA.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/12.piRNA_reads_n;   
    paste -d "\t" ./data_summary/12.piRNA_reads_1 ./data_summary/12.piRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/12.piRNA_reads;
    paste ./data_summary/12.piRNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/12.piRNA_reads - > ./data_summary/12.piRNA;

    tail -6 ./log/$i.Y_RNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/13.Y_RNA_reads_1; 
    tail -6 ./log/$i.Y_RNA.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/13.Y_RNA_reads_n;
    paste -d "\t" ./data_summary/13.Y_RNA_reads_1 ./data_summary/13.Y_RNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/13.Y_RNA_reads;
    paste ./data_summary/13.Y_RNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/13.Y_RNA_reads - > ./data_summary/13.Y_RNA;
    
    tail -6 ./log/$i.snRNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/14.snRNA_reads_1;
    tail -6 ./log/$i.snRNA.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/14.snRNA_reads_n;
    paste -d "\t" ./data_summary/14.snRNA_reads_1 ./data_summary/14.snRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/14.snRNA_reads;
    paste ./data_summary/14.snRNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/14.snRNA_reads - > ./data_summary/14.snRNA;
                                                                                                                                
    tail -6 ./log/$i.srpRNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/15.srpRNA_reads_1;             
    tail -6 ./log/$i.srpRNA.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/15.srpRNA_reads_n;
    paste -d "\t" ./data_summary/15.srpRNA_reads_1 ./data_summary/15.srpRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/15.srpRNA_reads;
    paste ./data_summary/15.srpRNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/15.srpRNA_reads - > ./data_summary/15.srpRNA;
                                                                                                                                                
    tail -6 ./log/$i.tRNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/16.tRNA_reads_1;
    tail -6 ./log/$i.tRNA.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/16.tRNA_reads_n;     
    paste -d "\t" ./data_summary/16.tRNA_reads_1 ./data_summary/16.tRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/16.tRNA_reads;
    paste ./data_summary/16.tRNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/16.tRNA_reads - > ./data_summary/16.tRNA;                               

    tail -6 ./log/$i.lncRNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/17.other_lncRNA_reads_1;    
    tail -6 ./log/$i.lncRNA.err.log | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/17.other_lncRNA_reads_n;            
    paste -d "\t" ./data_summary/17.other_lncRNA_reads_1 ./data_summary/17.other_lncRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/17.other_lncRNA_reads;
    paste ./data_summary/17.other_lncRNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/17.other_lncRNA_reads - > ./data_summary/17.other_lncRNA;                                                                                                                       
    tail -6 ./log/$i.mRNA.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/18.mRNA_reads_1;
    tail -6 ./log/$i.mRNA.err.log  | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/18.mRNA_reads_n;     
    paste -d "\t" ./data_summary/18.mRNA_reads_1 ./data_summary/18.mRNA_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/18.mRNA_reads;
    paste ./data_summary/18.mRNA_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/18.mRNA_reads - > ./data_summary/18.mRNA;

    tail -6 ./log/$i.other_human_genome.err.log | grep 'exactly' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/19.other_human_genome_reads_1;
    tail -6 ./log/$i.other_human_genome.err.log  | grep '>1' | head -1 | awk 'BEGIN{FS=OFS=" "}{print $1}' >> ./data_summary/19.other_human_genome_reads_n;
    paste -d "\t" ./data_summary/19.other_human_genome_reads_1 ./data_summary/19.other_human_genome_reads_n | awk 'total=(($1+$2)){print total}' > ./data_summary/19.other_human_genome_reads;
    paste ./data_summary/19.other_human_genome_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/19.other_human_genome_reads - > ./data_summary/19.other_human_genome;

    paste -d "\t" ./data_summary/11.miRNA_reads ./data_summary/12.piRNA_reads ./data_summary/13.Y_RNA_reads ./data_summary/14.snRNA_reads ./data_summary/15.srpRNA_reads ./data_summary/16.tRNA_reads ./data_summary/17.other_lncRNA_reads ./data_summary/18.mRNA_reads ./data_summary/19.other_human_genome_reads | awk 'total=(($1+$2+$3+$4+$5+$6+$7+$8+$9)){print total}' > ./data_summary/20.hg38_ordered_reads;
    paste -d "\t" ./data_summary/20.hg38_ordered_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/20.hg38_ordered_reads - > ./data_summary/20.hg38_ordered;

    paste ./data_summary/05.Kept_reads ./data_summary/20.hg38_ordered_reads | awk 'nonhuman=(($1-$2)){print nonhuman}' > ./data_summary/21.nonhuman_reads;
    paste ./data_summary/21.nonhuman_reads ./data_summary/05.Kept_reads | awk 'ratio=(($1/$2*100)){printf("%.2f\n",ratio)}' | sed 's/$/%/g' | paste -d "\t" ./data_summary/21.nonhuman_reads - > ./data_summary/21.nonhuman_genome;

paste -d "\t" ./sample_name ./data_summary/01.Raw_reads ./data_summary/03.Clean_reads ./data_summary/04.rRNA ./data_summary/05.Kept_reads ./data_summary/06.hg38 ./data_summary/11.miRNA ./data_summary/12.piRNA ./data_summary/13.Y_RNA ./data_summary/14.snRNA ./data_summary/15.srpRNA ./data_summary/16.tRNA ./data_summary/17.other_lncRNA ./data_summary/18.mRNA ./data_summary/19.other_human_genome ./data_summary/20.hg38_ordered ./data_summary/21.nonhuman_genome | cat ./data_summary/data_summary.header -  > ./data_summary/00.data_summary              
done
