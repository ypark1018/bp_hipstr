for x in `seq 1 22`; do \
cat chr$x.*individual.summary.txt | grep "TRUE DEVONO MUTATION" | cut -d' ' -f7 | paste -sd+ - | bc >> DM.snv; \
cat chr$x.*individual.summary.txt | grep "AMBIGUOUS DEVONO MUTATION" | cut -d' ' -f7 | paste -sd+ - | bc >> ADM.snv; \
cat chr$x.*individual.summary.txt | grep "MENDEL ERRORS IN DOUBLE TRIOS" | cut -d' ' -f9 | paste -sd+ - | bc >> ME.snv; \
done
