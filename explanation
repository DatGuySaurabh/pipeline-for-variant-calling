 Here is an explanation of each command in the script:

#!/bin/bash is called a shebang line. It specifies the interpreter to be used for the script that follows. In this case, bash is the interpreter.

gatk_root="/home/ubuntu/Downloads/gatk-4.2.6.1/"; is setting the value of the gatk_root variable to be the path to the gatk program.

picard_root="/home/ubuntu/Downloads/"; is setting the value of the picard_root variable to be the path to the picard program.

read1="${1}"; seq1=$(basename $read1) is setting the value of the read1 variable to be the first argument passed to the script, and then setting the value of seq1 to be the base name of read1.

read2="${2}"; seq2=$(basename $read2) is setting the value of the read2 variable to be the second argument passed to the script, and then setting the value of seq2 to be the base name of read2.

referance="${3}"; ref=$(basename $referance) is setting the value of the referance variable to be the third argument passed to the script, and then setting the value of ref to be the base name of referance.

known_variants="${4}"; kv=$(basename $known_variants) is setting the value of the known_variants variable to be the fourth argument passed to the script, and then setting the value of kv to be the base name of known_variants.

#fastqc $seq1; is a commented-out command to run the fastqc program on seq1.

#fastqc $seq2; is a commented-out command to run the fastqc program on seq2.

samtools faidx $ref; is a command to index the reference genome using the samtools faidx command.

bwa index $ref; is a command to index the reference genome using the bwa index command.

java -jar "${picard_root}"picard.jar CreateSequenceDictionary R=$ref; is a command to use java to run the picard program and create a sequence dictionary for the reference genome.

tabix -fp vcf $kv; is a command to index the known variants file using the tabix command.

bwa mem -M -t 4 $ref $seq1 $seq2 > "$seq1"_"$seq2"_mapped.sam; is a command to align seq1 and seq2 with the reference genome using the bwa mem command, and output the resulting alignment in SAM
samtools view -b "$seq1"_"$seq2"_mapped.sam -o "$seq1"_"$seq2"_mapped.bam; is a command to convert the alignment file from SAM to BAM format using the samtools view command.

samtools index "$seq1"_"$seq2"_mapped.bam; is a command to index the BAM file using the samtools index command.

samtools sort "$seq1"_"$seq2"_mapped.bam -o "$seq1"_"$seq2"_sorted.bam; is a command to sort the BAM file by coordinate using the samtools sort command.

java -jar "${picard_root}"picard.jar MarkDuplicates I="$seq1"_"$seq2"_sorted.bam O="$seq1"_"$seq2"_marked_duplicates.bam M="$seq1"_"$seq2"_marked_dup_metrics.txt; is a command to mark duplicates in the sorted BAM file using java and the picard program.

samtools index "$seq1"_"$seq2"_marked_duplicates.bam; is a command to create an index file for the marked duplicates BAM file using the samtools index command.

java -jar "${gatk_root}"gatk-package-4.2.6.1-local.jar BaseRecalibrator -R $ref -I "$seq1"_"$seq2"_marked_duplicates.bam -known-sites $kv -O "$seq1"_"$seq2"_recal_data.table; is a command to perform base quality score recalibration on the marked duplicates BAM file using java and the gatk program.

java -jar "${gatk_root}"gatk-package-4.2.6.1-local.jar ApplyBQSR -R $ref -I "$seq1"_"$seq2"_marked_duplicates.bam --bqsr-recal-file "$seq1"_"$seq2"_recal_data.table -O "$seq1"_"$seq2"_recal.bam; is a command to apply the base quality score recalibration to the marked duplicates BAM file using java and the gatk program.
"${gatk_root}"./gatk BaseRecalibrator -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG.bam -R $ref --known-sites $kv -O "$seq1"_"$seq2"_base_recal_1.table; is a command to perform base quality score recalibration on the BAM file using the gatk BaseRecalibrator program.

"${gatk_root}"./gatk ApplyBQSR -R $ref -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG.bam --bqsr-recal-file "$seq1"_"$seq2"_base_recal_1.table -O "$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam; is a command to apply the base quality score recalibration to the BAM file using the gatk ApplyBQSR program.

"${gatk_root}"./gatk BaseRecalibrator -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam -R $ref --known-sites $kv -O "$seq1"_"$seq2"_base_recal_2.table; is a command to perform a second round of base
"${gatk_root}"./gatk AnalyzeCovariates -before "$seq1"_"$seq2"_base_recal_1.table -after "$seq1"_"$seq2"_base_recal_1.table -plots "$seq1"_"$seq2"_recalibration_plot.pdf ; is a command to analyze the covarities of the base quality score recalibration using the gatk AnalyzeCovariates program.

java -jar "${picard_root}"picard.jar CollectAlignmentSummaryMetrics R="$ref" I="$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam O="$seq1"_"$seq2"_allignment_metrics.txt; is a command to collect alignment summary metrics for the BAM file using java and the picard program.

java -jar "${picard_root}"picard.jar CollectInsertSizeMetrics I="$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam O="$seq1"_"$seq2"_insert_metrics.txt H="$seq1"_"$seq2"_insert_histogram.pdf; is a command to collect insert size metrics for the BAM file using java and the picard program.

"${gatk_root}"./gatk HaplotypeCaller -R $ref -I "$seq1"_"$seq2"_mapped_sorted_duprem_RG_BQSR.bam -O "$seq1"_"$seq2"_raw_variants.vcf; is a command to generate a VCF file of raw variants using the gatk HaplotypeCaller program.

"${gatk_root}"./gatk VariantFiltration -R $ref -V "$seq1"_"$seq2"_raw_variants.vcf -O "$seq1"_"$seq2"_filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "my_snp_filter"; is a command to filter the raw variants VCF file using the gatk VariantFiltration program.
