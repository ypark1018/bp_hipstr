#make sure vcftools module is loaded
#module load vcftools

for x in `seq 1 22`; do \
vcf-concat chr$x*.vcf.gz | gzip -c > chr$x.hipSTR.vcf.gz; \
done