# plot the histogram
cd /home/zhangyanjun/projects/exRNA/hcc_examples
for i in `cat ./sample_name`
do
Rscript ./bin/plot_readsLens.R ./stat/$i.lengthN.stat.tsv ./stat/$i.lengthN.    stat.pdf
done
