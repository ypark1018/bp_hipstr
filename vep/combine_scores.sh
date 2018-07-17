dir=/u/home/p/parkyj/nobackup-eeskin2/bipolar_hipSTR/orig/filtered/multiallelic/vep

head -n 1 $dir/variant_scores1_1/chr1.allSTR.sum > tmp.head

for x in `seq 1 9`; do \
#create headers
cp -f tmp.head allSTR$x.sum; \
cp -f tmp.head del$x.sum; \
cp -f tmp.head allSTR$x.rv.sum; \
cp -f tmp.head del$x.rv.sum; \
cp -f tmp.head allSTR$x.cm.sum; \
cp -f tmp.head del$x.cm.sum; done

#concat
for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores1_1/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR1.sum; \
tail -n 1 $dir/variant_scores1_1/chr$x.del.sum | sed "s/$/\n/" >> del1.sum; \
tail -n 1 $dir/variant_scores1_1/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR1.rv.sum; \
tail -n 1 $dir/variant_scores1_1/chr$x.rv.del.sum | sed "s/$/\n/" >> del1.rv.sum; \
tail -n 1 $dir/variant_scores1_1/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR1.cm.sum; \
tail -n 1 $dir/variant_scores1_1/chr$x.cm.del.sum | sed "s/$/\n/" >> del1.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores1_2/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR2.sum; \
tail -n 1 $dir/variant_scores1_2/chr$x.del.sum | sed "s/$/\n/" >> del2.sum; \
tail -n 1 $dir/variant_scores1_2/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR2.rv.sum; \
tail -n 1 $dir/variant_scores1_2/chr$x.rv.del.sum | sed "s/$/\n/" >> del2.rv.sum; \
tail -n 1 $dir/variant_scores1_2/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR2.cm.sum; \
tail -n 1 $dir/variant_scores1_2/chr$x.cm.del.sum | sed "s/$/\n/" >> del2.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores1_3/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR3.sum; \
tail -n 1 $dir/variant_scores1_3/chr$x.del.sum | sed "s/$/\n/" >> del3.sum; \
tail -n 1 $dir/variant_scores1_3/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR3.rv.sum; \
tail -n 1 $dir/variant_scores1_3/chr$x.rv.del.sum | sed "s/$/\n/" >> del3.rv.sum; \
tail -n 1 $dir/variant_scores1_3/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR3.cm.sum; \
tail -n 1 $dir/variant_scores1_3/chr$x.cm.del.sum | sed "s/$/\n/" >> del3.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores2_1/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR4.sum; \
tail -n 1 $dir/variant_scores2_1/chr$x.del.sum | sed "s/$/\n/" >> del4.sum; \
tail -n 1 $dir/variant_scores2_1/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR4.rv.sum; \
tail -n 1 $dir/variant_scores2_1/chr$x.rv.del.sum | sed "s/$/\n/" >> del4.rv.sum; \
tail -n 1 $dir/variant_scores2_1/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR4.cm.sum; \
tail -n 1 $dir/variant_scores2_1/chr$x.cm.del.sum | sed "s/$/\n/" >> del4.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores2_2/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR5.sum; \
tail -n 1 $dir/variant_scores2_2/chr$x.del.sum | sed "s/$/\n/" >> del5.sum; \
tail -n 1 $dir/variant_scores2_2/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR5.rv.sum; \
tail -n 1 $dir/variant_scores2_2/chr$x.rv.del.sum | sed "s/$/\n/" >> del5.rv.sum; \
tail -n 1 $dir/variant_scores2_2/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR5.cm.sum; \
tail -n 1 $dir/variant_scores2_2/chr$x.cm.del.sum | sed "s/$/\n/" >> del5.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores2_3/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR6.sum; \
tail -n 1 $dir/variant_scores2_3/chr$x.del.sum | sed "s/$/\n/" >> del6.sum; \
tail -n 1 $dir/variant_scores2_3/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR6.rv.sum; \
tail -n 1 $dir/variant_scores2_3/chr$x.rv.del.sum | sed "s/$/\n/" >> del6.rv.sum; \
tail -n 1 $dir/variant_scores2_3/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR6.cm.sum; \
tail -n 1 $dir/variant_scores2_3/chr$x.cm.del.sum | sed "s/$/\n/" >> del6.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores3_1/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR7.sum; \
tail -n 1 $dir/variant_scores3_1/chr$x.del.sum | sed "s/$/\n/" >> del7.sum; \
tail -n 1 $dir/variant_scores3_1/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR7.rv.sum; \
tail -n 1 $dir/variant_scores3_1/chr$x.rv.del.sum | sed "s/$/\n/" >> del7.rv.sum; \
tail -n 1 $dir/variant_scores3_1/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR7.cm.sum; \
tail -n 1 $dir/variant_scores3_1/chr$x.cm.del.sum | sed "s/$/\n/" >> del7.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores3_2/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR8.sum; \
tail -n 1 $dir/variant_scores3_2/chr$x.del.sum | sed "s/$/\n/" >> del8.sum; \
tail -n 1 $dir/variant_scores3_2/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR8.rv.sum; \
tail -n 1 $dir/variant_scores3_2/chr$x.rv.del.sum | sed "s/$/\n/" >> del8.rv.sum; \
tail -n 1 $dir/variant_scores3_2/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR8.cm.sum; \
tail -n 1 $dir/variant_scores3_2/chr$x.cm.del.sum | sed "s/$/\n/" >> del8.cm.sum; done

for x in `seq 1 22`; do \
tail -n 1 $dir/variant_scores3_3/chr$x.allSTR.sum | sed "s/$/\n/" >> allSTR9.sum; \
tail -n 1 $dir/variant_scores3_3/chr$x.del.sum | sed "s/$/\n/" >> del9.sum; \
tail -n 1 $dir/variant_scores3_3/chr$x.rv.allSTR.sum | sed "s/$/\n/" >> allSTR9.rv.sum; \
tail -n 1 $dir/variant_scores3_3/chr$x.rv.del.sum | sed "s/$/\n/" >> del9.rv.sum; \
tail -n 1 $dir/variant_scores3_3/chr$x.cm.allSTR.sum | sed "s/$/\n/" >> allSTR9.cm.sum; \
tail -n 1 $dir/variant_scores3_3/chr$x.cm.del.sum | sed "s/$/\n/" >> del9.cm.sum; done

rm -f tmp.head