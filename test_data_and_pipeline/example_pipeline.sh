##Example pipeline for example data
##22 February 2015

##identify unique high coverage SNVs
echo "identifying SNVs"
python ../identify_high_cov_SNVs.py variants 10

##Determine coverage of bam files for at those SNV loci############################
##
## ls bams/* | xargs -I % -P 6 -n 1 bash get_cov.sh % high_cov_unique_SNVs.bed 10 &
##
###################################################################################

##Annotate SNVs at unique loci#####################################################
##
## perl table_annovar.pl ../high_cov_unique_SNVs.anno_input humandb/ -buildver hg19 \
## -out unique_snvs -remove -protocol refGene -operation g -nastring . -csvout
##
###################################################################################

echo "Creating SNVM"
##Create a SNVM with every high confidence SNV
python ../create_SNVM.py coverages/ freqs/ unique_snvs.hg19_multianno.csv


##Filter SNVM of SNVs that do not have enough samples with adequate coverage
## 11 is the first group size
## Removing SNVs that have more than 2 Na's in either each group
echo "Filtering SNVM"
python ../remove_NAs.py SNVM.unfiltered.csv 11 2 2

echo "Creating Training Sets (exonic)"
##Create training sets for exonic SNVs, at different sizes, for every combination of
##leave 2 out.
python ../uni_rank_and_set_creation.py SNVM.filtered.csv 11 10,25,50,75,100 exonic exonic_every_combo every_combo

echo "Creating Training Sets (UTRs)"
##Create training sets for UTR3 and UTR5 SNVs, at different sizes, for 50 random
## combinations of leave 2 out.
python ../uni_rank_and_set_creation.py SNVM.filtered.csv 11 10,25,50,75,100 UTR3,UTR5 UTRS_rand_subset_50 rand_subset 50
