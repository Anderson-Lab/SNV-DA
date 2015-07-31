Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv nonsyn nonsyn_parallel 10 1000 40 
Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv exonic exonic_parallel 10 1500 50 
Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv intergenic intergenic_parallel 10 500 25 
Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv intronic intronic_parallel 10 250 10 
Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv UTR3 UTR3_parallel 10 800 40 
Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv exonic,UTR3 exonic_UTR3_parallel 10 1500 50 
Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv ncRNA ncRNA_parallel 10 160 5 
Rscript parallel_find_opt_sPLS.R SNVM.unsup_filtered SNVM.unsup_loci.csv all all_parallel 10 1500 50 
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv nonsyn nonsyn_perm 530 1000 .8343
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv intergenic intergenic_perm 456 1000 .8434
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv intronic intronic_perm 20 1000 .8014
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv UTR3 UTR3_perm 42 1000 .7598
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv ncRNA ncRNA_perm 12 1000 .7043
Rscript parallel_permutation_sPLS.R SNVM.unsup_filtered.csv SNVM.unsup_loci.csv none all_perm 32 1000 .669
cd ~/pipelines/alignment/
python align_wrapper.py pvr_fastqs.txt ~/Raw/NGS/PvR/bams/








