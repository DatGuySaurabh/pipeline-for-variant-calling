# pipeline-for-variant-calling
Variant Calling for NGS data 
Here is an explanation of the script:

This script appears to be a pipeline for processing and aligning DNA sequence data. It takes four input files as arguments: two FASTQ files containing paired-end DNA sequence reads (read1 and read2), a reference genome in FASTA format (referance), and a VCF file containing known variants (known_variants). The script performs the following steps:

Quality check the input sequences. This step is commented out in the script, but it likely involves running the fastqc program on seq1 and seq2 to assess the quality of the reads.

Index the reference genome and create a sequence dictionary. This step involves running the samtools faidx and bwa index commands on the reference genome, and then using java and the picard program to create a sequence dictionary for the reference genome.

Index the known variants. This step involves running the tabix program on the VCF file to index it.

Align the input sequences with the reference genome. This step involves running the bwa mem command to align seq1 and seq2 with the reference genome and output the resulting alignment in SAM format.

Convert the alignment file from SAM to BAM format. This step involves running samtools view to convert the SAM file to BAM format.

Index the BAM file. This step involves running samtools index on the BAM file.

Sort the BAM file. This step involves running samtools sort on the BAM file to sort it by coordinate.

Mark duplicates in the sorted BAM file. This step involves using java and the picard program to mark duplicates in the BAM file.

Create an index file for the marked duplicates BAM file. This step involves running samtools index on the marked duplicates BAM file.

Perform base quality score recalibration on the marked duplicates BAM file. This step involves using java and the gatk program to recalibrate the base quality scores in the BAM file, using the known variants file as a reference.

The script includes several date statements which print the current date and time to the console. It also includes several echo statements which print messages to the console, such as the start and end times for each step in the pipeline and the names of the steps being performed.
