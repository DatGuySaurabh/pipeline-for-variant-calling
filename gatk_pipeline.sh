#!/bin/bash;
echo -e "\n\n\n";
echo -e "\n";
echo "#######################################";
echo "#variables and intialization#";
echo "#######################################";
echo -e "\n";


#Root Directories

#pipeline_root="";
gatk_root="/home/ubuntu/Downloads/gatk-4.2.6.1/";
picard_root="/home/ubuntu/Downloads/";
read1="${1}";
seq1=$(basename $read1)

read2="${2}";
seq2=$(basename $read2)


referance="${3}";
ref=$(basename $referance)

known_variants="${4}";
kv=$(basename $known_variants)

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check of read 1#";
echo "#######################################";
echo -e "\n";
#fastqc $seq1;
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#FASTA quality check of read 2#";
echo "#######################################";
echo -e "\n";
#fastqc $seq2;
date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#indexing reference genome#";
echo "#######################################";
echo -e "\n";
samtools faidx $ref;


bwa index $ref;

#create sequence dictionary
java -jar "${picard_root}"picard.jar CreateSequenceDictionary R=$ref;

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#indexing known variants#";
echo "#######################################";
echo -e "\n";

tabix -fp vcf $kv;


date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#mapping/alignment of reads with reference genome#";
echo "#######################################";
echo -e "\n";

bwa mem -M -t 4 $ref $seq1 $seq2 > "$seq1"_"$seq2"_mapped.sam;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#convert sam to bam#";
echo "#######################################";
echo -e "\n";

samtools view -b "$seq1"_"$seq2"_mapped.sam -o "$seq1"_"$seq2"_mapped.bam;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#indexing bam file#";
echo "#######################################";
echo -e "\n";

samtools index "$seq1"_"$seq2"_mapped.bam;
 
date +"%D %T":"##Process Completed Time##";

 
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#sorting bam file#";
echo "#######################################";
echo -e "\n";

samtools sort "$seq1"_"$seq2"_mapped.bam -o "$seq1"_"$seq2"_mapped_sorted.bam;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#remove duplicates#";
echo "#######################################";
echo -e "\n";

samtools rmdup "$seq1"_"$seq2"_mapped_sorted.bam "$seq1"_"$seq2"_mapped_sorted_duprem.bam;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#add read groups to bam file#";
echo "#######################################";
echo -e "\n";

java -jar "${picard_root}"picard.jar AddOrReplaceReadGroups I="$seq1"_"$seq2"_mapped_sorted_duprem.bam O="$seq1"_"$seq2"_mapped_sorted_duprem_RG.bam RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=sample;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#indexing final bam file#";
echo "#######################################";
echo -e "\n";

samtools index "$seq1"_"$seq2"_mapped_sorted_duprem_RG.bam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Displaying stats of BAM file#";
echo "#######################################";
echo -e "\n";

samtools stat "$seq1"_"$seq2"_mapped_sorted_duprem_RG.bam > "$seq1"_"$seq2"_stats.txt;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Base quality score Recalibration 1 (BQSR)#";
echo "#######################################";
echo -e "\n";

#Creating base recal table

"${gatk_root}"./gatk BaseRecalibrator -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG.bam -R $ref --known-sites $kv -O "$seq1"_"$seq2"_base_recal_1.table;

date +"%D %T":"##Process Completed Time##";


#Apply BQSR
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#apply BQSR#";
echo "#######################################";
echo -e "\n";

"${gatk_root}"./gatk ApplyBQSR -R $ref -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG.bam --bqsr-recal-file "$seq1"_"$seq2"_base_recal_1.table -O "$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Base quality score Recalibration 2 (BQSR)#";
echo "#######################################";
echo -e "\n";

#Creating base recal table 2 for covarities

"${gatk_root}"./gatk BaseRecalibrator -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam -R $ref --known-sites $kv -O "$seq1"_"$seq2"_base_recal_2.table;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Analyzing Covariaties#";
echo "#######################################";
echo -e "\n";

#Analyzing Covarities

"${gatk_root}"./gatk AnalyzeCovariates -before "$seq1"_"$seq2"_base_recal_1.table -after "$seq1"_"$seq2"_base_recal_1.table -plots "$seq1"_"$seq2"_recalibration_plot.pdf ;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Collecting Allignment Summary Metrics#";
echo "#######################################";
echo -e "\n";
#Collecting Alignment metrics
java -jar "${picard_root}"picard.jar CollectAlignmentSummaryMetrics R="$ref" I="$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam O="$seq1"_"$seq2"_allignment_metrics.txt;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Collect Insert size metrics#";
echo "#######################################";
echo -e "\n";
#Collecting Insert Metrics
java -jar "${picard_root}"picard.jar CollectInsertSizeMetrics I="$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam O="$seq1"_"$seq2"_insert_metrics.txt HISTOGRAM_FILE="$seq1"_"$seq2"_insert_metrics_HISTOGRAM.pdf;

date +"%D %T":"##Process Completed Time##";



echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Variant calling using Haplotypecaller#";
echo "#######################################";
echo -e "\n";

#haplotypecaller variant calling
"${gatk_root}"./gatk HaplotypeCaller  -R $ref -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam  -O "$seq1"_"$seq2"_raw_variants.vcf  --sample-name sample ;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Hard Filtration of Variants#";
echo "#######################################";
echo -e "\n";

#Vaariant hard filtering
"${gatk_root}"./gatk VariantFiltration  -R $ref -V "$seq1"_"$seq2"_raw_variants.vcf  -O "$seq1"_"$seq2"_Filtered_variants.vcf  -filter-name "QD_2" -filter "QD < 2.00"  -filter-name "DP_10" -filter "DP < 10.00" -filter-name "MQ_40" -filter "MQ < 40.00" -filter-name "verylow_quality" -filter "QUAL < 30.00"  -filter-name "low_quality" -filter "QUAL >30 && QUAL<50";

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Excluding Filtered Variables#";
echo "#######################################";
echo -e "\n";

#Excluding Filtered Variables
"${gatk_root}"./gatk SelectVariants  --exclude-filtered  -V "$seq1"_"$seq2"_Filtered_variants.vcf  -O "$seq1"_"$seq2"_final_variants.vcf ;

date +"%D %T":"##Process Completed Time##";


echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo -e "\n";
echo "#######################################";
echo "#Moving results to desktop#";
echo "#######################################";
echo -e "\n";

#creating result directory and moving bam , stats , recalibration_plot , histogram and final_vcf file .

mkdir /home/ubuntu/Desktop/"$seq1"_"$seq2"_variant_discovery_results

mv "$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam /home/ubuntu/Desktop/"$seq1"_"$seq2"_variant_discovery_results
mv "$seq1"_"$seq2"_stats.txt /home/ubuntu/Desktop/"$seq1"_"$seq2"_variant_discovery_results
mv "$seq1"_"$seq2"_recalibration_plot.pdf /home/ubuntu/Desktop/"$seq1"_"$seq2"_variant_discovery_results
mv "$seq1"_"$seq2"_insert_metrics_HISTOGRAM.pdf /home/ubuntu/Desktop/"$seq1"_"$seq2"_variant_discovery_results
mv "$seq1"_"$seq2"_final_variants.vcf /home/ubuntu/Desktop/"$seq1"_"$seq2"_variant_discovery_results
echo done 

exit 0;


















