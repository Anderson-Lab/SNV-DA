Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv nonsyn nonsyn_parallel 10 1000 40 
Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv exonic exonic_parallel 10 1500 50 
Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv intergenic intergenic_parallel 10 500 25 
Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv intronic intronic_parallel 10 250 10 
Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv UTR3 UTR3_parallel 10 800 40 
Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv exonic,UTR3 exonic_UTR3_parallel 10 1500 50 
Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv ncRNA ncRNA_parallel 10 160 5 
Rscript parallel_find_opt_sPLS.R SNVM.less_filtered.csv SNVM_loci_less.csv all all_parallel 10 1500 50 
