## length of mapped reads for each RNA types
cd /home/zhangyanjun/projects/exRNA/hcc_examples
mkdir ./stat
for i in `cat ./sample_name`
do
for j in `tail -9 ./05.mapping/mapping_order | head -8`
do
grep -v '@' ./05.mapping/$i/$j/$i.${j#*.}.sam | awk 'BEGIN{FS=OFS="\t"}($2!=4){print length($10)}' | sort -n | uniq -c | awk 'BEGIN{FS=" "; OFS="\t" }{print $2,$1}' | sort -nk1,1 | sed -e "s/^/$i\t${j#*.}\t/g" >> ./stat/$i.lengthN.stat.tsv
done
sed -i -e "1i sample\ttype\tlen\tnum" ./stat/$i.lengthN.stat.tsv
done
