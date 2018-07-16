for x in `cat chrlist`; \
do vcftools --vcf $x.include.recode.vcf --counts --out $x; \
python paste.py $x.frq.count >> count.alleles; done