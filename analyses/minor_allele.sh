#get maf
dir=/u/home/p/parkyj/nobackup-eeskin2/bipolar_hipSTR/orig/filtered/multiallelic
vcf=multiallelic.recode

for x in `seq 1 22`; 
do \
vcftools --gzvcf $dir/chr$x.$vcf.vcf.gz --freq --out $x; \
python minor_allele.py $x.frq 1 > chr$x.maf; done

rm -f *.frq