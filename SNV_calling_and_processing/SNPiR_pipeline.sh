#!/bzin/bash

#export PERL5LIB=$PERL5LIB:/software/SNPIR_pipeline/SNPiR/
#export PATH=$PATH:/software/bwa-0.7.12/
#export PERL5OPT=$PERL5OPT:/software/SNPIR_pipeline/SNPiR/

FILE1="$1"

sample_name="$2"
numT="$3"
mkdir inter_files
mkdir inter_files/$sample_name.dir
mkdir pipe_logs
mkdir variants
annot_files="$4"


#echo  "Trimming Files" >> pipe_logs/$sample_name.log

#java -jar /software/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $numT \
#-phred33 $FILE1 inter_files/$sample_name.dir/$sample_name.trimmed \
#ILLUMINACLIP:/software/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 \
#LEADING:3 TRAILING:3 \
#SLIDINGWINDOW:4:20 MINLEN:36

#make sure that you have already created a SNPiR index using bwa index
echo "Aligning trimmed forward file to SNPiR.fa" >> pipe_logs/$sample_name.log
bwa mem -t $numT $annot_files/SNPiR.fa $FILE1  \
> inter_files/$sample_name.dir/$sample_name.forward.sam


#Aligning reads to normal reference sequence to get accurate coverage estiamtes
nohup bwa mem -t $numT ~/Raw/NGS/genome.builds/hg19/hg19_ref.fa $FILE1  > inter_files/$sample_name.dir/$sample_name.normal.sam &

#Convert Coordinates
echo "Converting coordinates of aligned file to SNPiR format" >> pipe_logs/$sample_name.log
java -Xmx8g -classpath /software/SNPIR_pipeline/SNPiR/ convertCoordinates < inter_files/$sample_name.dir/$sample_name.forward.sam > inter_files/$sample_name.dir/$sample_name.converted.sam

#Convert sam files to sorted bam files
echo "Converting sam to bam" >> pipe_logs/$sample_name.log
samtools view -bS -F 4 inter_files/$sample_name.dir/$sample_name.converted.sam | samtools sort - inter_files/$sample_name.dir/$sample_name.converted.sorted

if [ -s inter_files/$sample_name.dir/$sample_name.converted.sorted.bam ];
then 
        rm inter_files/$sample_name.dir/$sample_name.converted.sam
fi
 
#Make sure that there is a reference dictionary created by picard before using picard tools
echo "Removing duplicates with Picard" >> pipe_logs/$sample_name.log
java -Xmx8g -jar /software/SNPIR_pipeline/software/picard-tools-1.110/MarkDuplicates.jar VERBOSITY=ERROR MAX_RECORDS_IN_RAM=30000000 INPUT=inter_files/$sample_name.dir/$sample_name.converted.sorted.bam \
OUTPUT=inter_files/$sample_name.dir/$sample_name.rmdup.bam METRICS_FILE=inter_files/$sample_name.dir/$sample_name.picard_info.txt \
REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

if [ -s inter_files/$sample_name.dir/$sample_name.rmdup.bam ];
then 
        rm inter_files/$sample_name.dir/$sample_name.converted.sorted.bam
fi

#Change options as necessary
echo "Add read group info (required for GATK)" >> pipe_logs/$sample_name.log
java -Xmx8g -jar /software/SNPIR_pipeline/software/picard-tools-1.110/AddOrReplaceReadGroups.jar VERBOSITY=ERROR INPUT=inter_files/$sample_name.dir/$sample_name.rmdup.bam OUTPUT=inter_files/$sample_name.dir/$sample_name.filtered_w_RG.bam RGLB=$sample_name RGPL=Illumina RGPU=Lane1 RGSM=$sample_name VALIDATION_STRINGENCY=SILENT

if [ -s inter_files/$sample_name.dir/$sample_name.filtered_w_RG.bam ];
then 
        rm inter_files/$sample_name.dir/$sample_name.rmdup.bam
fi


echo "reordering bam" >> pipe_logs/$sample_name.log
java -Xmx8g -jar /software/SNPIR_pipeline/software/picard-tools-1.110/ReorderSam.jar VERBOSITY=ERROR INPUT=inter_files/$sample_name.dir/$sample_name.filtered_w_RG.bam output=inter_files/$sample_name.dir/$sample_name.filtered.bam REFERENCE=$annot_files/SNPiR.fa VALIDATION_STRINGENCY=SILENT


if [ -s inter_files/$sample_name.dir/$sample_name.filtered.bam ];
then 
        rm inter_files/$sample_name.dir/$sample_name.filtered_w_RG.bam
fi

echo "create index for bam" >> pipe_logs/$sample_name.log
samtools index inter_files/$sample_name.dir/$sample_name.filtered.bam

touch  inter_files/$sample_name.dir/$sample_name.intervals

echo "indel realigner" >> pipe_logs/$sample_name.log
java -Xmx8g -jar /software/GenomeAnalysisTK.jar -T IndelRealigner -R $annot_files/hg19_ref.fa  -I inter_files/$sample_name.dir/$sample_name.filtered.bam -targetIntervals inter_files/$sample_name.dir/$sample_name.intervals -o inter_files/$sample_name.dir/$sample_name.indel.aligned.bam -U ALLOW_N_CIGAR_READS --maxReadsForRealignment 10000 

echo "BaseRecalibrator(creates grp file)" >> pipe_logs/$sample_name.log
java -Xmx8g -jar /software/GenomeAnalysisTK.jar  -T BaseRecalibrator -nct $numT -I inter_files/$sample_name.dir/$sample_name.indel.aligned.bam  -R $annot_files/hg19_ref.fa  -U ALLOW_N_CIGAR_READS   -knownSites $annot_files/dbSNP_hg19_20150605.vcf -o inter_files/$sample_name.dir/$sample_name.grp

echo "Print Reads (creates recalibrated bam from grp file)" >> pipe_logs/$sample_name.log
java -Xmx8g -jar /software/GenomeAnalysisTK.jar -T PrintReads -nct $numT -R $annot_files/hg19_ref.fa -I inter_files/$sample_name.dir/$sample_name.indel.aligned.bam -U ALLOW_N_CIGAR_READS -BQSR inter_files/$sample_name.dir/$sample_name.grp -o inter_files/$sample_name.dir/$sample_name.final.bam 


if [ -s inter_files/$sample_name.dir/$sample_name.final.bam ];
then 
        #rm inter_files/$sample_name.dir/$sample_name.filtered.bam
        rm inter_files/$sample_name.dir/$sample_name.indel.aligned.bam
fi

echo "UnifiedGenotyper" >> pipe_logs/$sample_name.log
java -Xmx8g -jar /software/GenomeAnalysisTK.jar -T UnifiedGenotyper  -nct $numT -R $annot_files/hg19_ref.fa -I inter_files/$sample_name.dir/$sample_name.final.bam --dbsnp $annot_files/dbSNP_hg19_20150605.vcf -o inter_files/$sample_name.dir/$sample_name.snps.raw.vcf -stand_call_conf 0 -stand_emit_conf 0  -U ALLOW_N_CIGAR_READS --output_mode EMIT_VARIANTS_ONLY -glm SNP 

echo "Convert vcf to SNPiR format (using SNPiR tool)" >> pipe_logs/$sample_name.log
sh /software/SNPIR_pipeline/SNPiR/convertVCF.sh inter_files/$sample_name.dir/$sample_name.snps.raw.vcf inter_files/$sample_name.dir/$sample_name.snpir.txt 20

python filter_low_cov_single.py 10 inter_files/$sample_name.dir/$sample_name.snpir.txt inter_files/$sample_name.dir/$sample_name.snpir.hicov.txt

echo "filter mismatches at read end" >> pipe_logs/$sample_name.log
perl /software/SNPIR_pipeline/SNPiR/filter_mismatch_first6bp.pl -infile inter_files/$sample_name.dir/$sample_name.snpir.hicov.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.trimmed.txt -bamfile inter_files/$sample_name.dir/$sample_name.final.bam

python filter_low_cov_single.py 10 inter_files/$sample_name.dir/$sample_name.snpir.trimmed.txt inter_files/$sample_name.dir/$sample_name.snpir.hicov.trimmed.txt

echo "Remove variants that lie in repetitive regions" >> pipe_logs/$sample_name.log
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' inter_files/$sample_name.dir/$sample_name.snpir.hicov.trimmed.txt | intersectBed -a stdin -b $annot_files/RepeatMasker.bed -v | cut -f1,3-7 > inter_files/$sample_name.dir/$sample_name.snpir.rmsk.txt

echo "Filter variants in intronic regions (on this step)" >> pipe_logs/$sample_name.log
perl /software/SNPIR_pipeline/SNPiR/filter_intron_near_splicejuncts.pl -infile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.txt -genefile $annot_files/all_gene_annot.bed

echo "filter variants in homopolymers" >> pipe_logs/$sample_name.log
perl /software/SNPIR_pipeline/SNPiR/filter_homopolymer_nucleotides.pl -infile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.txt -refgenome $annot_files/hg19_ref.fa
#
#
echo "filter variants that were caused by mismapped reads" >> pipe_logs/$sample_name.log
perl /software/SNPIR_pipeline/SNPiR/BLAT_candidates.pl -infile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.rmblat.txt -bamfile inter_files/$sample_name.dir/$sample_name.final.bam -refgenome $annot_files/hg19_ref.fa
#
#
echo "remove known RNA editing sites" >> pipe_logs/$sample_name.log
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.rmblat.txt | intersectBed -a stdin -b $annot_files/rna_editing_sites.bed -v > variants/$sample_name.variants.bed 

#echo "It finished running! Happy hunting!" >> pipe_logs/$sample_name.log
