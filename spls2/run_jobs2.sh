#Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv none all_parallel 10 1500 50 
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv nonsyn nonsyn_perm 488 1000 .801
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv intronic intronic_perm 18 1000 .805
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv UTR3 UTR3_perm 36 1000 .730
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv ncRNA ncRNA_perm 13 1000 .689


