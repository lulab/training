#!/bin/bash
#This script compares MD5 of files with known MD5 strings and shows the names of truncated files.
#Last edited by ZLMZ, 2018.3.29

echo -e "Now try to find the truncated file.\n"
for file in file*
do
	md5_new=$(md5sum $file | cut -d " " -f 1)
	md5_old=$(grep "$file" md5sum.txt | cut -d " " -f 4)
	[ $md5_new != $md5_old ] && echo -e "The file $file is truncated.\n"
done
echo "Job completed."



