#!/bin/bash

#load R module before running script
#module load R/3.4.1
all_chr=all_chr.converted #this is just the header for the file
converted_prefix=CONVERTED.chr
R_dir=/u/home/p/parkyj/nobackup-eeskin2/scripts/R/concordance
java_dir=/u/home/p/parkyj/nobackup-eeskin2/scripts/java

#make a copy so we don't overwrite existing file
cp $all_chr $all_chr.length.txt
for x in `seq 1 22`; do cat $converted_prefix$x >> $all_chr.length.txt; done
Rscript $R_dir/convert_old_format.R > /dev/null 2>&1
Rscript $R_dir/convert_new_format.R > /dev/null 2>&1
java -Xmx1G -cp $java_dir RecodeBothOldNewSTR old_str.txt new_str.txt old_str.data new_str.data old_str.allele new_str.allele > /dev/null 2>&1
Rscript $R_dir/compare_STR.R > /dev/null 2>&1