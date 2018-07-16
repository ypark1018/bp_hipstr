for x in `seq 1 22`; do \
for y in $(seq 0 $(($(ls chr$x | wc -l)-2))); do \
echo $x $y >> hg19.hipstr.reference.list; \
done; done