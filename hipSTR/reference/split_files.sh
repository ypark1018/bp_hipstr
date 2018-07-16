for x in `seq 1 22`; do \
mkdir chr$x; \
grep chr$x hg19.hipstr_reference.bed > chr$x/chr$x.hg19.hipstr_reference.bed; \
split -l 10000 -a 2 -d chr$x/chr$x.hg19.hipstr_reference.bed chr$x/chr$x.hg19.hipstr.reference.; \
done