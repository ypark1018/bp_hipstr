dir=/u/home/p/parkyj/nobackup-eeskin2/bipolar_hipSTR/orig/filtered/multiallelic

for x in `seq 1 22`; \
do vcftools --gzvcf $dir/chr$x.multiallelic.recode.vcf.gz --counts --out chr$x; \
cat chr$x.frq.count >> count.alleles; done