#get the allele frequencies at each loci
#input name of output het.[output]
file_suffix="het"
dir=/u/home/p/parkyj/nobackup-eeskin2/bipolar_hipSTR/orig/filtered/multiallelic
vcf=multiallelic.recode
for x in `seq 1 22`; \
do \
vcftools --gzvcf $dir/chr$x.$vcf.vcf.gz --freq --out $x; \
grep -v CHROM $x.frq > $x.tmp.frq; \
python heterozygosity.py $x.tmp.frq > chr$x.$file_suffix; done

rm -f *.frq