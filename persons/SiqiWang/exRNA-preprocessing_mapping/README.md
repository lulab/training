
# Snakemake package for preprocesing and mapping of exRNA-seq
please click: https://github.com/SiqiWang94/Preprocessing_and_mapping_of_exRNA-seq_data

> Shared by Siqi Wang

> Email: zzuwsq@163.com

|File|Description|Function|
|:-------------|:---------------------------|:--------|
|Snakefile|Snakemake workflows|Python scripts extended by declarative code to define rules. Rules describe how to create output files from input files.| 
|config.json| Standard configuration file |To make your workflows more flexible|
|cluster.json| Cluster configuration file  |To specify cluster submission parameters outside the Snakefile|
|jobscript.pbs| Auxilliary file| To submit your Snakemake commandto the cluster vis the `qsub` or `bsub`| 
|run.sh| Bash file| Record the snakemake command to run snakefile on your system|
|log|Log files direction| To store the log files|

All scripts have been tested on Lu Lab Server(IBME)
