This document is the guide to calling STR using hipSTR and perform preliminary analyses which include the following:
1. hipSTR
2. QC for hipSTR output
3. Isolate polymorphic loci
4. Identify mendelian errors & correct to missing
5. Minor allele frequency (MAF)
6. Heterozygosity
7. De novo mutation analysis
8. Burden Analysis

**for all python scripts list of possible input is displayed when run with no input

Contact Youngjun Park at parkyj2350@gmail.com for any questions.


#              1. hipSTR               #

JOB FILE: ./hipSTR/hipSTR.job
    	  ./hipSTR/concat.sh
REFERENCE: https://github.com/tfwillems/HipSTR

INPUT:
--bams			comma separated list of BAM files generated from BWA-MEM, sorted and indexed using samtools
--fasta			FASTA file containing sequence for each chromosome in BED file
--regions		BED file containing coordinates for each STR region of interst. The BED file in script is divided for parallelization
--str-vcf		name of the output file
--chrom			chromosome number; can be omited if you want to call the entire genome at once <- would not recommend

NOTE:
The BED file was divided into 226 different regions so 226 jobs total will be submitted
If SNP vcf is accessible, we can use it to get a more accurate phasing information (more detailed info on github above)

After hipSTR.job is finished, combine all vcfs for each chromosome using catcat.sh


#                2. QC                 #

JOB FILE: ./hipSTR/filter.job
REFERENCE: https://github.com/tfwillems/HipSTR

INPUT:
filter.job
--dir			path to vcf files
--vcf			name of vcf files
--out_dir		path to output
--filter		filtering python script provided by hipSTR: filter_vcf.py
--outliers		File containing specific loci to be removed from the output vcf. List one sample per line.

filter_vcf.py
--vcf			vcf file
--min-call-qual		minimum Q score for each call (notated under Q in INFO in vcf file)
--max-call-flak-indel	total number of reads containing an indel in the regions flanking the STR
--max-call-stutter	total number of reads at a locus with what HipSTR thinks is a stutter artifact
--min-loc-calls		minimum number of calls required for each loci (call rate)


#     3. Isolate Polymorphic Loci      #

JOB FILE: ./hipSTR/multiallelic.job
    	  ./hipSTR/biallelic.job
INPUT:
--dir			path to vcf files
--vcf			name of vcf
--out_dir		path to output

NOTE:
Make sure vcftools is loaded. multiallelic.job filters for loci containing 2 or more unique alleles; biallelic.job filters for loci containing 2 unique alleles.


#       4. Identify ME/Correct         #

JOBFILE: ./mendel_errors/get.me.job
	 ./mendel_errors/set.me.missing.job
	 ./mendel_errors/me.sh

JAVA SCRIPTS: ConvertVCFToLinkage
     	      SetLinkageFileMEMissingPedCheck2
	      ConvertPedToGIGILongFormat
	      MergePedCheckMEStatFile
	      ConvertGIGIOutputVCF

INPUT:
get.me.job & set.me.missing.job
--java			path to java program
--javacp		path to java scripts
--vcftools		path to vcftools software
--pedcheck		path to pedcheck software

--origdir		path to family info data
--pedinfo		pedfile containing family structure
--pedstruct		pedfile containing family structure + sample ID etc
--idmapfile		^
--familyfile		file containing list of different families

--dir			path to vcf files
--name			name of vcf file is in this format chr$chr.$name.vcf.gz
--outdir		path to output

--seqfam		another ped info file (ask JaeHoon)

NOTE:
Use get.me.job to acquire mendelian errors and me.sh to compile the results into one file. The file you want to look at for the ME is *.level.1.pedcheck.err.gz.
Compiled output will be in file called mendel.errors

Use set.me.missing.job to set the genotypes with mendelian errors to missing and filter out the loci with missing rate > 10%.


#      5. Minor Allele Frequency       #

JOB FILE: ./analyses/minor_allele.sh

PYTHON SCRIPT: minor_allele.py

NOTE:
self explanatory - acquire MAF for each loci. For multiallelic variants MAF is the frequency of second most prevalent allele in a variant


#          6. Heterozygosity           #

JOB FILE: ./analyses/heterozygosity.sh

PYTHON SCRIPT: heterozygosity.py

INPUT:
--file_suffix		output file will be in this format chr$chr.$file_suffix
--dir			path to vcf
--vcf			vcf file

NOTE:
self explanatory - acquire heterozygosity for each variant


#     7. De Novo Mutation Analysis     #

JOB FILE: ./denovo/denovovcf.job
    	  ./denovo/denovo.job
	  ./denovo/denovo.SNV.job
	  ./denovo/add_DM.sh

PYTHON SCRIPT: ./python/denovo/denovo.2.py

INPUT:
denovovcf.job
--f_dir			path to filtered files
--output_dir		path to output

--vcf_file		path to VCF after QC
--keep_file		path to VCF after correcting genotypes with mendelian errors to missing

denovo.job
--f_dir			path to directory of vcf file
--vcf_file		name of vcf file
--output_dir		path to output
--family		file containing conversion from family ID to sample ID on each line: [FAMILY ID] [SAMPLE ID]
--ped			pedfile with family info
--denovo		path to python script

NOTE:
denovovcf.job filters out the VCF after QC with the list of loci present after correcting genotypes with mendelian errors to missing (output of set.me.missing.job). After filtering, only loci with little mendelian inconsistencies should be left.


#          8. Burden Analysis          #

JOB FILE: ./vep/vep.str.job
    	  ./vep/vep_parser.job
	  ./vep/map.variants.v2.job
	  ./vep/rv_overlap.job
	  ./vep/vep_filter.job
	  ./vep/variant_score.job
	  ./vep/combine_scores.sh

PYTHON SCRIPTS: ./python/vep_parser.py
       		./python/map.variants.v2.py
		./python/minor_allele.py
		./python/rv_overlap.py
		./python/variant_score.1.py
		./python/variant_score.2.py
		./python/variant_score.3.py

R SCRIPT: ./R/burden_test.R

INPUT:
variant_score.*.py
python variant_score.*.py [input vcf] [impact level: MODIFIER=0, LOW=1, MODERATE=2, HIGH=3] [output prefix]

JAVA SCRIPTS: ParseVEPVCFFile

software: VEP

NOTE:
- Use vep.str.job to annotate the variants using VEP
- vep_parser.job parses the vcf file generated from vep.str.job to output in easily readable format
- map_variants.v2.job maps vcf to known gene locations provided by user
- rv_overlap.job outputs the list of rare variants from 1KGP and hipSTR. Specific thresholds can be set for each data set can be set in the job file. It will look for loci present in 1KGP first and look for loci in hipSTR that are not present int 1KGP data
- vep_filter.job filters out the vep vcf based on list of overlapping genes provided and rare variants provided by user
- variant_score.job performs burden analyses based on following categories respectively: presence of deleterious genes, sum of lengths of deleterious alleles, sum of differences in length from reference allele (variance)
- combine_scores.sh is a helper script to compile all chromosomes into one file -> each line is sum score for each chromosome
-burden_test.R takes the files created from combine_scores.sh and runs a linear mixed model; requires the following files:
1. ./vep/burden_test/indiv_pheno_info.txt
2. ./vep/burden_test/PedFileBPDB.txt
3. ./vep/burden_test/admixing.csv