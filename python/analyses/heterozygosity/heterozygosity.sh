#get the allele frequencies at each loci
#input name of output het.[output]
#change lines 1, 3
file_suffix="bi.processed"
for x in `cat chrlist`; \
do vcftools --vcf $x.bi.processed.recode.vcf --freq --out $x; \
cat $x.frq >> $file_suffix; done

grep -v CHROM $file_suffix > $file_suffix.frq
python heterozygosity.py $file_suffix.frq > het.$file_suffix

rm -f *.frq $file_suffix