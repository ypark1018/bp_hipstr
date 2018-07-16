rm -f *.freq
for x in `cat chrlist`; do \
    vcftools --vcf ${x}.bi.recode.vcf --freq --out $x; \
    #get minor allele frequencies
    python minor_allele.py $x.frq > $x.freq; done
rm -f *.frq
cat *.freq > minor_alleles.bi.freq