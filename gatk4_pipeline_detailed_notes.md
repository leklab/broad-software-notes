## GATK4+ pipeline notes
---
- By: Prech Brian Uapinyoying 
- Created: 08/10/2018
- Last updated: 08/22/2018

### Requirements and resources

1. Notes: 
    - <span style="color:red">**STILL INCOMPLETE**</span>- problems with output from haplotypecaller step
    - Most of the early steps were performed locally in the official GATK4 docker container that uses Ubuntu, so any install commands will be for a Debian distribution (see part [Part XIII](#part-xiii-installing-and-running-gatk40-in-docker-to-run-locally-so-that-you-can-try-out-commands-without-starting-any-jobs-on-the-cluster))

2. Requirements:
    - Whole exome, whole genome or RNA-seq data
    - Fastq files or pre-aligned BAMs
    - Paired-end sequencing data
    - GATK/4.0.6.0
    - BWA aligner
    - UCSC hg19 from broad reference_bundle

3. Resources
    - [Documentation for first half of pipeline in GATK4](https://gatkforums.broadinstitute.org/gatk/discussion/7899/reference-implementation-pairedendsinglesamplewf-pipeline)
    - [Detailed tutorial of each step (from GATK3)](https://software.broadinstitute.org/gatk/documentation/topic?name=tutorials&v=3)
    - [What are gVCFs compared to regular VCFs?](https://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf)

### Index
1. [Part I. Convert input data into an unmapped BAM (uBAM) as input to the pipeline.](#part-i-convert-input-data-into-an-unmapped-bam-ubam-as-input-to-the-pipeline)
2. [Part II. Mark Illumina adapters using `MarkIlluminaAdapters`](#part-ii-mark-illumina-adapters-markilluminaadapters)
3. [Part III. Convert uBAM to fastq prior to BWA and discount adapter sequences using `SamToFastq`](#part-iii-convert-ubam-to-fastq-prior-to-bwa-and-discount-adapter-sequences-using-samtofastq)
4. [Part IV. Alignment to reference genome using `BWA mem`](#part-iv-alignment-to-reference-genome-using-bwa-mem)
5. [Part V. Quality measurements: `CollectQualityYieldMetrics` & `CollectAlignmentSummaryMetrics`](#part-v-quality-measurements-collectqualityyieldmetrics-collectalignmentsummarymetrics)
6. [Part VI. Mark duplicates](#part-vi-mark-duplicates)
7. [Part VII. Sort the bam file by coordinate and fix specific tags using `SortSam` and `SetNmAndUqTags`](#part-vii-sort-the-bam-file-by-coordinate-and-fix-specific-tags-using-sortsam-and-setnmanduqtags)
8. [Part VIII. Recalibrate base quality scores to produce analysis ready reads](#part-viii-recalibrate-base-quality-scores-to-produce-analysis-ready-reads)
9. [Part VIV. Call variants per sample using `HaplotypeCaller` in gVCF mode.](#part-viv-call-variants-per-sample-using-haplotypecaller-in-gvcf-mode)
10. [Part X. Process a small call set (of gVCFs) with `GenotypeGVCFs`](#part-x-process-a-small-call-set-of-gvcfs-with-genotypegvcfs)
11. [Part XI. Perform quality controls on the final BAM file to get metrics and see if there is contamination](#part-xi-perform-quality-controls-on-the-final-bam-file-to-get-metrics-and-see-if-there-is-contamination)
12. [Part XII. Other useful tools](#part-xii-other-useful-tools)
13. [Part XIII. Installing and running GATK4.0 in docker to run locally so that you can try out commands without starting any jobs on the cluster](#part-xiii-installing-and-running-gatk40-in-docker-to-run-locally-so-that-you-can-try-out-commands-without-starting-any-jobs-on-the-cluster)

### Part I. Convert input data into an unmapped BAM (uBAM) as input to the pipeline.

[Full instructions for this part from Broad](https://gatkforums.broadinstitute.org/gatk/discussion/6484/how-to-generate-an-unmapped-bam-from-fastq-or-aligned-bam)

1. Go download example data sets from the tutorial to try these out on your docker image.
    - [tutorial_6484_FastqToSam.tar.gz](https://drive.google.com/open?id=0BzI1CyccGsZiUXVNT3hsNldvUFk)
    - [tutorial_6464_RevertSam.tar.gz](https://drive.google.com/open?id=0BzI1CyccGsZiMWZacmVWWnV2VFE)
    - The aligned files were mapped to `human_g1k_v37_decoy.fasta.gz`, I couldn't find that, but you can use any reference like hg19

2. If your raw data are FASTQ files, use PicardTools `FastqToSam`, otherwise skip to step 3.
    - Picard's `FastqToSam` transforms a FASTQ file to an unmapped BAM (uBAM), requires two read group fields and makes optional specification of other read group fields. In the command below we note which fields are required for GATK Best Practices Workflows. All other read group fields are optional.
```bash
# example through PicardTools (for older versions of GATK):
# With GATK4.0, we can use the wrapper script. Don't forget to specify
# the java arguments.  If you have multiple args, wrap them in quotes like so:
# gatk --java-options "-Xmx4G -XX:+PrintGCDetails" [program arguments]
# You also need to change the argument syntax to "kabob" style e.g. --OUTPUT or -O

gatk --java-options "-Xmx8G" FastqToSam \
    --FASTQ /gatk/my_data/test_data/tutorial_6484_FastqToSam/6484_snippet_1.fastq \
    --FASTQ2 /gatk/my_data/test_data/tutorial_6484_FastqToSam/6484_snippet_2.fastq \
    --OUTPUT /gatk/my_data/results/fastqtosam/6484_snippet_fastqtosam.bam \ # Consider naming <sample>.unmapped.bam
    --READ_GROUP_NAME H0164.2 \
    --SAMPLE_NAME NA12878 \
    --LIBRARY_NAME Solexa-272222 \
    --PLATFORM_UNIT H0164ALXX140820.2 \
    --PLATFORM illumina \
    --SEQUENCING_CENTER BI \
    --STRIP_UNPAIRED_MATE_NUMBER=true \
    --RUN_DATE 2014-08-20T00:00:00-0400
```
    - on the cluster Slurm run
```bash
#! /bin/bash

gatk --java-options "-Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m" FastqToSam \
    --FASTQ /gatk/my_data/test_data/tutorial_6484_FastqToSam/6484_snippet_1.fastq \
    --FASTQ2 /gatk/my_data/test_data/tutorial_6484_FastqToSam/6484_snippet_2.fastq \
    --OUTPUT /gatk/my_data/results/fastqtosam/6484_snippet_fastqtosam.bam \ # Consider naming <sample>.unmapped.bam
    --READ_GROUP_NAME H0164.2 \
    --SAMPLE_NAME NA12878 \
    --LIBRARY_NAME Solexa-272222 \
    --PLATFORM_UNIT H0164ALXX140820.2 \
    --PLATFORM illumina \
    --SEQUENCING_CENTER BI \
    --STRIP_UNPAIRED_MATE_NUMBER true \
    --RUN_DATE 2014-08-20T00:00:00-0400
```

3. If your raw data are aligned BAM files from another source use PicardTool's `RevertSam`.
    - This process takes more work, because we need to clear specific tags and sanitize/discard problematic records that that cause problems for certain downstream tools e.g. MarkIlluminaAdapters
    - What attributes should you clear from the bam file? Take a look at some of the tags in your bam file after the main commands below for explainations
```bash
# Using GATK4 wrapper tool for docker example
gatk --java-options "-Xmx8G" RevertSam \
    -I /gatk/my_data/test_data/tutorial_6484_RevertSam/6484_snippet.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_revertsam.bam \ # Consider naming <sample>.unmapped.bam
    --SANITIZE true \
    --MAX_DISCARD_FRACTION 0.005 \
    --ATTRIBUTE_TO_CLEAR XT \
    --ATTRIBUTE_TO_CLEAR XN \
    --ATTRIBUTE_TO_CLEAR AS \
    --ATTRIBUTE_TO_CLEAR OC \
    --ATTRIBUTE_TO_CLEAR OP \
    --SORT_ORDER queryname \
    --RESTORE_ORIGINAL_QUALITIES true \
    --REMOVE_DUPLICATE_INFORMATION true \
    --REMOVE_ALIGNMENT_INFORMATION true
    # --TMP_DIR /<path>/tmp

# Ruddle example:
gatk --java-options "-Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m" RevertSam \
    -I /home/pbu2/scratch60/UC394/UC394.1.bam \
    -O /home/pbu2/scratch60/results/2018-08-03/s1_revertsam/UC394.1.unmapped.bam \
    --SANITIZE true \
    --MAX_DISCARD_FRACTION 0.005 \
    --ATTRIBUTE_TO_CLEAR XT \
    --ATTRIBUTE_TO_CLEAR XN \
    --ATTRIBUTE_TO_CLEAR AS \
    --ATTRIBUTE_TO_CLEAR OC \
    --ATTRIBUTE_TO_CLEAR OP \
    --SORT_ORDER queryname \
    --RESTORE_ORIGINAL_QUALITIES true \
    --REMOVE_DUPLICATE_INFORMATION true \
    --REMOVE_ALIGNMENT_INFORMATION true
    # --TMP_DIR /<path>/tmp
/home/pbu2/scratch60/results/2018-08-03/s1_revertsam

# Monkol's params for his script are very similar except:
# ATTRIBUTE_TO_CLEAR = X0, X1, XA, XC, XG, XM, XN, XO, XT, AM, AS, BQ, CC, CP,
# E2, H0, H1, H2, HI, IH, MF, NH, OC, OP, PQ, R2, S2, SM, SQ, U2, XQ
# MAX_DISCARD_FRACTION = 0.01
```
    - To check an see what tags are in your aligned BAM file, you can run
```bash
samtools view 6484_snippet.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}'

# SAM Tags found in this particular file (There are more):
# ---------------------------------------
# MC - * CIGAR string for mate/next segment
# MD - * String for mismatching positions
# PG - * Program
# RG - Read group
# NM - * Edit distance to the reference
# MQ - * Mapping quality of the mate/next segment
# OQ - Original base quality
# UQ - * Phred likelihood of the segment, conditional on the mapping being correct
# AS - * Genome assembly identifier.
# SA - * Other canonical alignments in a chimeric alignment
# XT - tags starting with X?, Y?, Z? are reserved for end users
# OC - Original CIGAR
# OP - Original mapping position
# XN - tags starting with X?, Y?, Z? are reserved for end users

# * Standard tags cleared by default by revert SAM
# Non-standard tags need to be removed manually using the `ATTRIBUTE_TO_CLEAR` flag
```
    - Read the SAM/BAM file specifications for more details on tags: 
        + [Samtools specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
        + [Documentation on SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf)

    - `SANITIZE` flag - Downstream tools will have problems with paired reads with missing mates, duplicated records, and records with mismatches in length of bases and qualities. Any paired reads file subset for a genomic interval requires sanitizing to remove reads with lost mates that align outside of the interval.

    - In this command, we've set `MAX_DISCARD_FRACTION` to a more strict threshold of 0.005 instead of the default 0.01. Whether or not this fraction is reached, the tool informs you of the number and fraction of reads it discards. This parameter asks the tool to additionally inform you of the discarded fraction via an exception as it finishes processing.

    - Optional, but it's recommended that you create a subset of the aligned BAM file and test out `RevertSam` first so it won't take so long. To do this, we can use the `PrintReads` tool. This requires the the reference file too. [See this discussion thread](https://gatkforums.broadinstitute.org/gatk/discussion/6517).

### Part II. Mark Illumina adapters `MarkIlluminaAdapters`

```bash
# Picard example:
java -Xmx8G -jar /path/picard.jar MarkIlluminaAdapters \
    I=6483_snippet_revertsam.bam \
    O=6483_snippet_markilluminaadapters.bam \
    M=6483_snippet_markilluminaadapters_metrics.txt \ # naming required
    TMP_DIR=/path/shlee #optional to process large files

# Gatk4 example:
# /gatk/my_data/results/revertsam
gatk --java-options "-Xmx8G" MarkIlluminaAdapters \
    -I /gatk/my_data/results/revertsam/6484_snippet_revertsam.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    -M /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters_metrics.txt \
    --ADAPTERS PAIRED_END
    # --TMP_DIR <path>/tmp
```

### Part III. Convert uBAM to fastq prior to BWA and discount adapter sequences using `SamToFastq`

```bash
java -Xmx8G -jar /path/picard.jar SamToFastq \
    I=6483_snippet_markilluminaadapters.bam \
    FASTQ=6483_snippet_samtofastq_interleaved.fq \
    CLIPPING_ATTRIBUTE=XT \ # take those marked illumina adapters
    CLIPPING_ACTION=2 \ # and change quality scores to 2?
    # Interleave = F and R reads are arranged one below the other in a single fastq
    INTERLEAVE=true \ # BWA can accept interleave fastq with -p option
    NON_PF=true \ # Include reads that do not pass QC
    TMP_DIR=/path/shlee #optional to process large files  

# Gatk4
gatk --java-options "-Xmx8G" SamToFastq \
    -I /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    --FASTQ /gatk/my_data/results/revertsam/6484_snippet_samtofastq_interleave.fastq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INTERLEAVE true \
    --NON_PF true
    # --TMP_DIR <path>/tmp
```

### Part IV. Alignment to reference genome using `BWA mem`

1. If you are using the GATK in the docker container, it will not have several peices of software installed. Docker images are linux based, so there is apt-get. You can install basic tools like vim and nano command line text editors.
```bash
apt-get update
apt-get install vim
apt-get install nano
```

2. BWA aliner is not pre-installed in the docker container and its not in the apt-get repository. You will need to install it yourself.
```bash
mkdir /usr/local/tools
mkdir /usr/local/tools/bwa
cd /usr/local/tools/bwa
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
tar -xvjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make # compiles the source code
./bwa # test to see if it runs
cd /usr/local/bin 
ln -s /usr/local/tools/bwa/bwa-0.7.17/bwa bwa # create a soft-link so we can run bwa from any directory
bwa # see if it works when in a different directory
```

3. Before we can align the sequence we need the reference genome and index it using BWA.  
    - For this example/tutorial data, I will use hg19/gch37 reference + decoy sequences to speed things up
    - Read more about [decoy genomes](http://www.cureffi.org/2013/02/01/the-decoy-genome/) and why it speeds up alignments
```bash
cd /gatk/my_data/references/hg19_decoy
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
bwa index hs37d5.fa # took me 1.5 hrs on a core i7, 3.1 ghz, 10 gb ram allocation
samtools faidx hs37d5.fa # create the hs37d5.fa.fai
```
    - However, if you need the real genome only, you can grab a copy of the BWA indexed hg19 right off of Ruddle. Copy it to your `my_data` directory that is mounted with your docker container.
    - Or broad's [hg38 build](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/)
```bash
cd /Users/uapinyoyingpb/Notebook/Archive/2018/molecular_diagnostics/gatk_pipeline
scp -r pbu2@ruddle.hpc.yale.edu:/home/bioinfo/genomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/version_0.6.0/ .

# You also want to index the genome.fa file (fasta), you get a genome.fa.fai file
cd BWAIndex/version_0.6.0/
samtools faidx genome.fa
```
    - You can also download the iGenomes pre-indexed genomes from the [Illumina website](https://support.illumina.com/sequencing/sequencing_software/igenome.html)
    - [Notes on alternate index files](https://software.broadinstitute.org/gatk/documentation/article?id=8017#2)

4. Time to do the actual alignment with `BWA mem`
```bash
bwa mem -M -t 7 -p /gatk/my_data/references/hg19_decoy/hs37d5.fa /gatk/my_data/results/revertsam/6484_snippet_samtofastq_interleave.fastq > 6484_snippet_bwa_mem.sam

# Real time: 67.013 sec; CPU: 152.360 sec
```
    - `-M` to flag shorter split hits as secondary. This is optional for Picard compatibility as MarkDuplicates can directly process BWA's alignment, whether or not the alignment marks secondary hits. However, if we want MergeBamAlignment to reassign proper pair alignments, to generate data comparable to that produced by the Broad Genomics Platform, then we must mark secondary alignments.
    - `-p` to indicate the given file contains interleaved paired reads.
    - `-t` followed by a number for the number of processor threads to use concurrently. Here we use seven threads which is one less than the total threads available on my Mac laptap. Check your server or system's total number of threads with the following command `getconf _NPROCESSORS_ONLN`.  If your computer has hyper-threading you would have 2x the number of physical cores. (for me this is 4 physical cores and 8 cores)

5. Restore altered data and apply & adjust meta information using `MergeBamAlignment`
    - [First create a fasta dict file](https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference)
```bash
java -jar CreateSequenceDictionary.jar R=Homo_sapiens_assembly18.fasta O=Homo_sapiens_assembly18.dict

# GATK 4.0
gatk --java-options "-Xmx8G" CreateSequenceDictionary -R /gatk/my_data/references/hg19_decoy/hs37d5.fa -O /gatk/my_data/references/hg19_decoy/hs37d5.dict
```
    - Now run the actual `MergeBamALignment` script
```bash
java -Xmx16G -jar /path/picard.jar MergeBamAlignment \
    R=/path/Homo_sapiens_assembly19.fasta \
    UNMAPPED_BAM=6383_snippet_revertsam.bam \
    ALIGNED_BAM=6483_snippet_bwa_mem.sam \ #accepts either SAM or BAM
    O=6483_snippet_mergebamalignment.bam \
    CREATE_INDEX=true \ #standard Picard option for coordinate-sorted outputs
    ADD_MATE_CIGAR=true \ #default; adds MC tag
    CLIP_ADAPTERS=false \ #changed from default
    CLIP_OVERLAPPING_READS=true \ #default; soft-clips ends so mates do not extend past each other
    INCLUDE_SECONDARY_ALIGNMENTS=true \ #default
    MAX_INSERTIONS_OR_DELETIONS=-1 \ #changed to allow any number of insertions or deletions
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \ #changed from default BestMapq
    ATTRIBUTES_TO_RETAIN=XS \ #specify multiple times to retain tags starting with X, Y, or Z 
    TMP_DIR=/path/shlee #optional to process large files

# GATK4.0 version + monkol's options
gatk --java-options "-Xmx8G" MergeBamAlignment \
    -R /gatk/my_data/references/hg19_decoy/hs37d5.fa \
    --UNMAPPED_BAM /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    --ALIGNED_BAM /gatk/my_data/results/revertsam/6484_snippet_bwa_mem.sam \
    -O /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment.bam \
    --CREATE_INDEX true \
    --ADD_MATE_CIGAR true \
    --CLIP_ADAPTERS false \
    --CLIP_OVERLAPPING_READS true \
    --INCLUDE_SECONDARY_ALIGNMENTS true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ATTRIBUTES_TO_RETAIN XS \
    --ATTRIBUTES_TO_RETAIN X0 \ # options from monkol scripts
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --PAIRED_RUN true \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    #--MAX_RECORDS_IN_RAM=2000000 \ # maybe for clusters with high ram alloc
    --SORT_ORDER unsorted \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true  
    #--TMP_DIR=/<path>/tmp # optional to process large files
```
Easy to paste GATK4 all in one line:
```bash
gatk --java-options "-Xmx8G" MergeBamAlignment -R /gatk/my_data/references/hg19_decoy/hs37d5.fa --UNMAPPED_BAM /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam --ALIGNED_BAM /gatk/my_data/results/revertsam/6484_snippet_bwa_mem.sam -O /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment.bam --CREATE_INDEX true --ADD_MATE_CIGAR true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --ATTRIBUTES_TO_RETAIN X0 --VALIDATION_STRINGENCY SILENT --EXPECTED_ORIENTATIONS FR --PAIRED_RUN true --IS_BISULFITE_SEQUENCE false --ALIGNED_READS_ONLY false --SORT_ORDER unsorted --UNMAPPED_READ_STRATEGY COPY_TO_TAG --ALIGNER_PROPER_PAIR_FLAGS true --UNMAP_CONTAMINANT_READS true  
```
    - Lots of interesting details about what this actually does! https://software.broadinstitute.org/gatk/documentation/article?id=6483#step3C

### Part V. Quality measurements: `CollectQualityYieldMetrics` & `CollectAlignmentSummaryMetrics`

1. `CollectQualityYieldMetrics` 
    - [Overview of the tool](https://broadinstitute.github.io/picard/command-line-overview.html#CollectQualityYieldMetrics)
    - [Metric definitions](https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectQualityYieldMetrics.QualityYieldMetrics)
```bash
java -jar picard.jar CollectQualityYieldMetrics \
       I=input.bam \
       O=quality_yield_metrics.txt \

# GATK4
gatk --java-options "-Xmx8G" CollectQualityYieldMetrics \
    -I /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters_qualityyieldmetrics.txt
```

2. `CollectAlignmentSummaryMetrics`
    - Produces a summary of alignment metrics from a SAM or BAM file. This tool takes a SAM/BAM file input and produces metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters. Note that these quality filters are specific to Illumina data
```bash
java -jar picard.jar CollectAlignmentSummaryMetrics \
          R=reference_sequence.fasta \
          I=input.bam \
          O=output.txt

# GATK4 
gatk --java-options "-Xmx8G" CollectAlignmentSummaryMetrics \
    -R /gatk/my_data/references/hg19_decoy/hs37d5.fa \
    -I /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment_alignmentsummarymetrics.txt 
```

### Part VI. Mark duplicates

1. Most of the time, we will use MarkDuplicates, but there is also MarkDuplicatesWithMateCigar for special applications
    - [Background on MarkDuplicates vs MarkDuplicatesWithMateCigar](https://software.broadinstitute.org/gatk/documentation/article?id=6747)
```bash
java -Xmx32G -jar picard.jar MarkDuplicates \
    INPUT=6747_snippet.bam \ #specify multiple times to merge 
    OUTPUT=6747_snippet_markduplicates.bam \
    METRICS_FILE=6747_snippet_markduplicates_metrics.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \ #changed from default of 100
    CREATE_INDEX=true \ #optional
    TMP_DIR=/tmp

# GATK4 version:
gatk --java-options "-Xmx8G" MarkDuplicates \
    --INPUT /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment.bam \
    --OUTPUT /gatk/my_data/results/revertsam/6484_snippet_markduplicates.bam \
    --METRICS_FILE /gatk/my_data/results/revertsam/6484_snippet_markduplicates_metrics.txt \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER queryname \
    --CREATE_MD5_FILE true \
    --VALIDATION_STRINGENCY SILENT \
    --CREATE_INDEX true
    #--TMP_DIR=/tmp
```

### Part VII. Sort the bam file by coordinate and fix specific tags using `SortSam` and `SetNmAndUqTags`

1. `SortSam`
    - [Documentation on the SortSam tool](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.0/picard_sam_SortSam.php)
```bash
java -jar picard.jar SortSam \
    I=input.bam \
    O=sorted.bam \
    SORT_ORDER=coordinate

gatk --java-options "-Xmx8G" SortSam \
    -I /gatk/my_data/results/revertsam/6484_snippet_markduplicates.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_markduplicates_sorted.bam \
    --SORT_ORDER coordinate
```

2. `SetNmMdAndUqTags`
    - [Documentation on the SetNmMdAndUqTags](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_SetNmAndUqTags.php)
    - This tool takes in a coordinate-sorted SAM or BAM and calculatesthe NM, MD, and UQ tags by comparing with the reference.
    - This may be needed when MergeBamAlignment was run with SORT_ORDER other than 'coordinate' and thuscould not fix these tags then. The input must be coordinate sorted in order to run. If specified,the MD and NM tags can be ignored and only the UQ tag be set.
```bash
java -jar picard.jar SetNmMdAndUqTags 
    R=reference_sequence.fasta \
    I=sorted.bam \
    O=fixed.bam 

gatk --java-options "-Xmx8G" SetNmMdAndUqTags \
    -R /gatk/my_data/references/hg19_decoy/hs37d5.fa \
    -I /gatk/my_data/results/revertsam/6484_snippet_markduplicates_sorted.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_sorted_fixed.bam
```

### Part VIII. Recalibrate base quality scores to produce analysis ready reads

[What is base quality recalibration?](https://software.broadinstitute.org/gatk/documentation/article?id=2801)

1. Analyze patterns of covariation in the sequence dataset. 
    - This creates a GATKReport file called recal_data.table containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data.
    - You need to provide the program with a set of known sites or it won't run. The known sites are used to build the covariation model and estimate empirical base qualities. 
    - For each sample we are going to recalibrate by
```bash
# Example:
java -jar GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R reference.fa \
    -I input_reads.bam \
    -L 20 \ # restricts analysis to chromosome 20 (needs to match ref file nomenclature)
    -knownSites dbsnp.vcf \
    -knownSites gold_indels.vcf \
    -o recal_data.table \
    --useOriginalQualities \
    -U \ # What does this depreciated command do ?
    --disable_auto_index_creation_and_locking_when_reading_rods # not avail in gatk4

# GATK4, I don't have the dbsnp files for the decoy genome
gatk --java-options "-Xmx8G" BaseRecalibrator \
    -R reference.fa \
    -I input_reads.bam \
    -known-sites dbsnp.vcf \ # dashes instead of camel case
    -known-sites gold_indels.vcf \
    -O chr1_recal_data.table \ # big O
    # --disable_auto_index_creation_and_locking_when_reading_rods \ # Nonger available in gatk4
    -L 1 \ # depends on chromosome
    --use-original-qualities # dashes instead of camel case
```
    - We want to break up the recalibration of the bam files into chunks of the genome. Monkol divided it in this Perl script by processing the larger chromosomes alone, and combining the smaller ones together. Basically, we will run the above command multiple times on the same BAM file but with different -L and -o for each batch of chromosomes. 
    - Ideally, each batch would be run on a different node on the cluster so they can be processed in parallel. We will bring them back together later using the gather commands
```perl
%chr_lookup = (
    "1" => "-L 1",
    "2" => "-L 2",
    "3" => "-L 3",
    "4" => "-L 4",
    "5" => "-L 5",
    "6" => "-L 6",
    "7" => "-L 7",
    "8" => "-L 8",
    "9" => "-L 9",
    "10" => "-L 10",
    "11" => "-L 11",
    "12_13" => "-L 12 -L 13",
    "14_15" => "-L 14 -L 15",
    "16_17" => "-L 16 -L 17",
    "18_21" => "-L 18 -L 19 -L 20 -L 21",
    "22_X"  => "-L 22 -L X",
    "Y_NC"  => "-L Y -L MT -L GL000207.1 -L GL000226.1 -L GL000229.1 \
    -L GL000231.1 -L GL000210.1 -L GL000239.1 -L GL000235.1 -L GL000201.1 \
    -L GL000247.1 -L GL000245.1 -L GL000197.1 -L GL000203.1 -L GL000246.1 \
    -L GL000249.1 -L GL000196.1 -L GL000248.1 -L GL000244.1 -L GL000238.1 \
    -L GL000202.1 -L GL000234.1 -L GL000232.1 -L GL000206.1 -L GL000240.1 \
    -L GL000236.1 -L GL000241.1 -L GL000243.1 -L GL000242.1 -L GL000230.1 \
    -L GL000237.1 -L GL000233.1 -L GL000204.1 -L GL000198.1 -L GL000208.1 \
    -L GL000191.1 -L GL000227.1 -L GL000228.1 -L GL000214.1 -L GL000221.1 \
    -L GL000209.1 -L GL000218.1 -L GL000220.1 -L GL000213.1 -L GL000211.1 \
    -L GL000199.1 -L GL000217.1 -L GL000216.1 -L GL000215.1 -L GL000205.1 \
    -L GL000219.1 -L GL000224.1 -L GL000223.1 -L GL000195.1 -L GL000212.1 \
    -L GL000222.1 -L GL000200.1 -L GL000193.1 -L GL000194.1 -L GL000225.1 \
    -L GL000192.1 -L NC_007605"
);
```

2. Do a second pass to analyze covariation remaining after recalibration
    - This creates another GATKReport file, which we will use in the next step to generate plots. Note the use of the `-BQSR` flag, which tells the GATK engine to perform on-the-fly recalibration based on the first recalibration data table.
```bash
java -jar GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R reference.fa \
    -I input_reads.bam \
    -L 20 \
    -knownSites dbsnp.vcf \
    -knownSites gold_indels.vcf \
    -BQSR recal_data.table \ # From first pass!
    -o post_recal_data.table 

# GATK4
gatk --java-options "-Xmx8G" BaseRecalibrator \
    -R reference.fa \
    -I input_reads.bam \
    -L 1 \
    -known-sites Homo_sapiens_assembly19.dbsnp138.vcf \
    -known-sites Homo_sapiens_assembly19.known_indels.vcf \
    -known-sites Homo_sapiens_assembly19.variantEvalGoldStandard.vcf \
    -BQSR chr1_recal_data.table \
    -o chr1_post_recal_data.table
```
    - make sure your feature files (vcfs) are indexed `gatk IndexFeatureFile`
```bash
gatk --java-options "-Xmx8G" IndexFeatureFile --feature-file Homo_sapiens_assembly19.known_indels.vcf
gatk --java-options "-Xmx8G" IndexFeatureFile --feature-file Homo_sapiens_assembly19.variantEvalGoldStandard.vcf
```

3. Generate before/after comparison plots
    - This generates recalibration_plots.pdf that show how the reported base qualities match up to the empirical qualities calculated by the BaseRecalibrator. Comparing the before and after plots allows you to check the effect of the base recalibration process before you actually apply the recalibration to your sequence data. For details on how to interpret the base recalibration plots, please see the online GATK documentation.
```bash
java -jar GenomeAnalysisTK.jar \
    -T AnalyzeCovariates \
    -R reference.fa \
    -L 20 \
    -before recal_data.table \ # from first step
    -after post_recal_data.table \ # from second step
    -plots recalibration_plots.pdf

# GATK4 
gatk --java-options "-Xmx8G" AnalyzeCovariates \
    -R reference.fa \
    -L 1
    -before chr1_recal_data.table \
    -after chr1_post_recal_data.table \
    -plots chr1_recalibration_plots.pdf
```

4. Gather all the base quality recalibration reports together using `GatherBqsrReports`.
    - We split the base recalibration by chromosomes or chromosome groups so now we need to bring them all back into a single file.  This combined file will be used to recalibrate each of the broken up bam files
```bash
java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.GatherBqsrReports \
          -I input.list \ # a list of the recal_data.table for each -L
          -O output.grp

#GATK
gatk --java-options "-Xmx8G" GatherBQSRReports \
    -I input.list \
    -O combined_recal_data.table
```

5. Apply the recalibration to your sequence data
    - This creates a file called `recal_reads.bam` containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down. However, you have the option to retain them by adding the flag `–emit_original_quals` to the PrintReads command, in which case the original qualities will also be written in the file, tagged `OQ`.
```bash
java -jar GenomeAnalysisTK.jar \
    -T PrintReads \
    -R reference.fa \
    -I input_reads.bam \
    -L 1 \
    -BQSR recal_data.table \
    -o recal_reads.bam 

# GATK4 no longer uses PrintReads for applying base recalibration
# Instead there is a dedicated ApplyBQSR function
gatk --java-options "-Xmx8G" ApplyBQSR \
    -R reference.fa \
    -I input_reads.bam \
    -L 1 \
    --bqsr-recal-file combined_recal_data.table \
    -O chr1_recal_reads.bam # switch to big `O`
```
    - Notice how this step uses a very simple tool, PrintReads, to apply the recalibration. What’s happening here is that we are loading in the original sequence data, having the GATK engine recalibrate the base qualities on-the-fly thanks to the `-BQSR` flag (as explained earlier), and just using PrintReads to write out the resulting data to the new file.


6. Gather all the base recalibrated BAM files together using `filesGatherBamFiles`
```bash
java -jar picard.jar GatherBamFiles \
      I=input1.bam \ # must be specified in order, 1 per line
      I=input2.bam \
      O=gathered_files.bam

# GATK4, Can also provide a list of bams in txt file, 1 per line
gatk --java-options "-Xmx8G" GatherBamFiles \
    -I chr1_recal_reads.bam \
    -I chr2_recal_reads.bam \
    -I chr3_recal_reads.bam \
    -I chr4_recal_reads.bam \ # and so on until chrY etc
    -O gathered_files.bam
```

### Part VIV. Call variants per sample using `HaplotypeCaller` in gVCF mode.

```
# Monkol's pipeline
batch_vcf.pl 
    calls --> make_gvcf.pl 
                calls --> 1. scatter_intervals.pl (IntervalListTools)
                calls --> 2. haplotype_caller.pl  (HaplotypeCaller) 
                calls --> 3. gather_and_index.pl  (GatherVcfs)
```

1. Break up the regions you plan to call variants on into lists of even chunks / intervals
    - [Documentation on IntervalListTools](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.1/picard_util_IntervalListTools.php)
```bash
gatk --java-options "-Xmx8G" IntervalListTools \
    --PADDING 0 \
    --SCATTER_COUNT 10 \ # make 10 interval files
    --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --UNIQUE true \
    --SORT true \
    --BREAK_BANDS_AT_MULTIPLES_OF 0 \
    --INPUT /ycga-gpfs/project/ysm/lek/ml2529/reference/exome_calling_regions.v1.interval_list \
    --OUTPUT <dir>/ # This needs to be a directory not a file
```
    - `--PADDING` - The amount to pad each end of the intervals by before other operations are undertaken. Negative numbers are allowed and indicate intervals should be shrunk. Resulting intervals < 0 bases long will be removed. Padding is applied to the interval lists (both INPUT and SECOND_INPUT, if proivided) before the ACTION is performed.
    - `SCATTER_COUNT` - The number of files into which to scatter the resulting list by locus; in some situations, fewer intervals may be emitted. Note - if > 1, the resultant scattered intervals will be sorted and uniqued. The sort will be inverted if the INVERT flag is set.
    - `--BREAK_BANDS_AT_MULTIPLES_OF` - If set to a positive value will create a new interval list with the original intervals broken up at integer multiples of this value. Set to 0 to NOT break up intervals
    - `--UNIQUE` - If true, merge overlapping and adjacent intervals to create a list of unique intervals. Implies SORT=true.
    - `--SORT` - If true, sort the resulting interval list by coordinate.
    - `--INPUT` - One or more interval lists. If multiple interval lists are provided the output is theresult of merging the inputs. Supported formats are interval_list and VCF.
    - `--OUTPUT` - The output interval list file to write (if SCATTER_COUNT == 1), or the directory into which to write the scattered interval sub-directories (if SCATTER_COUNT > 1).
    - `--SUBDIVISION_MODE` - Selects between various ways in which scattering of the interval-list can happen
        + `BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW` - A scatter approach that tries to balance the number of bases in each interval list by estimating the remaining interval lists sizes. This is computed from the total number of unique bases and the bases we have consumed. This means that the interval list with the most number of unique bases is at most the ideal split length larger than the smallest interval list (unique number of bases).
        + [Documentation on subdivision mode](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.1/picard_util_IntervalListTools.php#--SUBDIVISION_MODE)

2. Run HaplotypeCaller
    - [Overview of HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/article?id=11068)
    - [Detailed documentation](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)
    - You want to run HaplotypeCaller on each sample for the number of SCATTER_COUNT (10) you defined in the previous step.
```python
# Pseudocode
for <interval_file> in list_of_interval_files:
    run (haplotypeCaller.py -L <interval_file>)
```
    - Actual command:
```bash
# GATK4
gatk --java-options "-Xmx8G" HaplotypeCaller \
    -L interval_list \
    -I input.bam \
    -O output.raw.snps.indels.vcf \ # changed form 'o'
    -R /ycga-gpfs/project/ysm/lek/ml2529/reference/Homo_sapiens_assembly19.fasta \
    -ERC GVCF \
    --min-pruning 3 \ # changed from '--minPruning'
    --max-alternate-alleles 3 \ # changed from '--max_alternate_alleles'
    --max-num-haplotypes-in-population 200 \ # changed from '--maxNumHaplotypesInPopulation'
    -contamination 0.0 \
    --read-filter OverclippedReadFilter \ # changed from '--read_filter OverclippedRead'
    -bamout output_HC_local_realigned.bam
    
    # These args seem to be gone in the transition from gatk3.4 to gatk4
    # --disable_auto_index_creation_and_locking_when_reading_rods \
    # -variant_index_parameter 128000 \
    # -variant_index_type LINEAR \
    # -pairHMM VECTOR_LOGLESS_CACHING \
```
    - `-ERC` / `--emitRefConfidence` - Mode for emitting reference confidence scores
    - `--minPruning` - Minimum support to not prune paths in the graph
    - `--max_alternate_alleles` - Maximum number of alternate alleles to genotype
    - `--maxNumHaplotypesInPopulation` - Maximum number of haplotypes to consider for your population
    - ` -contamination` / `--contamination_fraction_to_filter` - Fraction of contamination to aggressively remove
    - `--read-filter` / `-RF` - Read filters to be applied before analysis. This argument may be specified 0 or more times.
    - `-bamout` / `--bam-output` - File to which assembled haplotypes should be written
    - `--contamination-fraction-per-sample-file` / `-contamination-file` - Tab-separated File containing fraction of contamination in sequencing data (per sample) to aggressively remove. Format should be "" (Contamination is double) per line; No header.


3. Now combine all the output vcfs for all of the samples together and generate index using `GatherVcfs`
    - Gathers multiple VCF files from a scatter operation into a single VCF file. Input files must be supplied in genomic order and must not have events at overlapping positions.
    - Simple little class that combines multiple VCFs that have exactly the same set of samples and totally discrete sets of loci.
```bash
gatk --java-options "-Xmx8G" GatherVcfs \
    --INPUT $out_dir/$sample/$sample.vcf.gz.list \
    --OUTPUT $out_dir/$sample.vcf.gz
    # --CREATE_INDEX true \ # no longer in gatk4

# what is this?
print "tabix $out_dir/$sample.vcf.gz\n";
`tabix $out_dir/$sample.vcf.gz`;
```
    - `--CREATE_INDEX` - Whether to create a BAM index when writing a coordinate-sorted BAM file

### Part X. Process a small call set (of gVCFs) with `GenotypeGVCFs`

1. GenotypeGVCFs merges gVCF records that were produced as part of the Best Practices workflow for variant discovery using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the HaplotypeCaller, or result from combining such gVCF files using CombineGVCFs. 
    - This tool performs the multi-sample joint aggregation step and merges the records together in a sophisticated manner: at each position of the input gVCFs, this tool will combine all spanning records, produce correct genotype likelihoods, re-genotype the newly merged record, and then re-annotate it.
    - Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods for every position) that GenotypeGVCFs requires for its operation.
    - If the gVCF files contain allele specific annotations, add -G Standard -G AS_Standard to the command line.
```bash
gatk --java-options "-Xmx8G" GenotypeGVCFs \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -R /ycga-gpfs/project/ysm/lek/ml2529/reference/Homo_sapiens_assembly19.fasta \
    -D /ycga-gpfs/project/ysm/lek/ml2529/reference/Homo_sapiens_assembly19.dbsnp138.vcf \
    -V sample1.g.vcf \
    -V sample2.g.vcf \
    -o output.vcf \
    -L /ycga-gpfs/project/ysm/lek/ml2529/reference/exome_calling_regions.v1.interval_list
```

### Part XI. Perform quality controls on the final BAM file to get metrics and see if there is contamination

```
# Monkol's Pipeline
bam_qc.pl 
    calls --> collect_hs_metrics.pl.       (CollectHsMetrics)
    calls --> contamination_check.pl       (CollectMultipleMetrics)
    calls --> collect_multiple_metrics.pl  (verifyBamID)
```

1. Collect QC metrics pertaining to hybrid-selection methods (exome-capture) using `CollectHsMetrics`
    - This tool takes a SAM/BAM file input and collects metrics that are specific for sequence datasets generated through hybrid-selection. Hybrid-selection (HS) is the most commonly used technique to capture exon-specific sequences for targeted sequencing experiments such as exome sequencing
    - his tool requires an aligned SAM or BAM file as well as bait and target interval files in Picard interval_list format. You should use the bait and interval files that correspond to the capture kit that was used to generate the capture libraries for sequencing, which can generally be obtained from the kit manufacturer. If the baits and target intervals are provided in BED format, you can convert them to the Picard interval_list format using Picard's BedToInterval tool
    - If a reference sequence is provided, this program will calculate both AT_DROPOUT and GC_DROPOUT metrics. Dropout metrics are an attempt to measure the reduced representation of reads, in regions that deviate from 50% G/C content. This reduction in the number of aligned reads is due to the increased numbers of errors associated with sequencing regions with excessive or deficient numbers of G/C bases, ultimately leading to poor mapping efficiencies and lowcoverage in the affected regions.
    - If you are interested in getting G/C content and mean sequence depth information for every target interval, use the PER_TARGET_COVERAGE option.
    - Note: Metrics labeled as percentages are actually expressed as fractions!
```bash
gatk --java-options "-Xmx8G" CollectHsMetrics \
    --REFERENCE_SEQUENCE /ycga-gpfs/project/ysm/lek/ml2529/reference/Homo_sapiens_assembly19.fasta \
    --INPUT $bam \
    --OUTPUT out_dir/$sample.hybrid_selection_metrics \
    --TARGET_INTERVALS /ycga-gpfs/project/ysm/lek/ml2529/reference/exome_calling_regions.v1.interval_list \
    --BAIT_INTERVALS /ycga-gpfs/project/ysm/lek/ml2529/reference/exome_calling_regions.v1.interval_list \
    --TMP_DIR /ycga-gpfs/scratch60/ysm/lek/ml2529/tmp
```
    - `--BAIT_INTERVALS` - An interval list file that contains the locations of the baits used.
    - `--TARGET_INTERVALS` - An interval list file that contains the locations of the targets.

2. Collect multiple classes of QC metrics using `CollectMultipleMetrics`
    - [Documentation](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_analysis_CollectMultipleMetrics.php)
    - This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time to cut down on the time spent reading in data from input files. Available modules include: 
        + CollectAlignmentSummaryMetrics
        + CollectInsertSizeMetrics
        + QualityScoreDistribution 
        + MeanQualityByCycle 
        + CollectBaseDistributionByCycle
        + CollectGcBiasMetrics
        + RnaSeqMetrics
        + CollectSequencingArtifactMetrics
        + CollectQualityYieldMetrics. 
    - The tool produces outputs of '.pdf' and '.txt' files for each module, except for the CollectAlignmentSummaryMetrics module, which outputs only a '.txt' file. Output files are named by specifying a base name (without any file extensions). (may require R libraries to gen pdfs?)
```bash
gatk --java-options "-Xmx8G" CollectMultipleMetrics \
    --TMP_DIR /ycga-gpfs/scratch60/ysm/lek/ml2529/tmp \
    --REFERENCE_SEQUENCE /ycga-gpfs/project/ysm/lek/ml2529/reference/Homo_sapiens_assembly19.fasta \
    --INPUT input.bam \
    --OUTPUT $out_dir/$sample. \
    --PROGRAM null \
    --PROGRAM MeanQualityByCycle \
    --PROGRAM QualityScoreDistribution \
    --PROGRAM CollectInsertSizeMetrics \
    --PROGRAM CollectAlignmentSummaryMetrics \
    --METRIC_ACCUMULATION_LEVEL null \
    --METRIC_ACCUMULATION_LEVEL ALL_READS \
    --ASSUME_SORTED true
```

3. Contamination check using `VerifyBamID` (a non-GATK/Picard tool)
    - [Documentation](https://genome.sph.umich.edu/wiki/VerifyBamID#Basic_Usage)
    - verifyBamID is a software that verifies whether the reads in particular file match previously known genotypes for an individual (or group of individuals), and checks whether the reads are contaminated as a mixture of two samples. 
    - It can detect sample contamination and swaps when external genotypes are available. When external genotypes are not available, verifyBamID still robustly detects sample swaps.
```bash
VerifyBamID --vcf /ycga-gpfs/project/ysm/lek/ml2529/reference/ExomeContam.vcf \
    --verbose \
    --ignoreRG \
    --bam sample.bam \ #works on cram too
    --out sample # prefix
```
    - Each option provides the following features:
        + `--vcf` : specify required VCF file
        + `--bam` : specify required BAM file (indexed with .bam.bai or .bai file)
        + `--free-mix` : (default) Estimate contamination using sequence-only method with Brent's single dimensional optimization.
        + `--chip-mix` : (default) Estimate contamination using sequence+array method with Brent's single dimensional optimization.
        + `--ignoreRG` : ignore the read grouup level comparison and compare samples only (recommended for an expedited run)
        + `--out` : output file prefix (required)
        + `--verbose` : print the progress of the method on the screeen

### Part XII. Other useful tools

1. Convert BAM to CRAM for storage (reduces file size in half!). 
    - Make sure you have samtools version 1.x, not the really old 0.1.x versions
    - [Download link for samtools, bcftools and htslib](http://www.htslib.org/download/)
```bash
samtools view -C -T /gatk/my_data/references/hg19_decoy/hs37d5.fa /gatk/my_data/results/revertsam/6484_snippet_sorted_fixed.bam | \
    tee 6484_snippet_sorted_fixed.cram | \
    md5sum > 6484_snippet_sorted_fixed.cram.md5 && \
    samtools index 6484_snippet_sorted_fixed.cram && \
    mv 6484_snippet_sorted_fixed.cram.crai 6484_snippet_sorted_fixed.crai
```

2. Validate BAM files to see if there are errors
```bash
gatk --java-options "-Xmx8G" ValidateSamFile \
    -I input.bam \
    --MODE SUMMARY
```

3. Creating a python virtualenv in conda to use python 3 and snakemake
```bash
module load GATK/4.0.6.0
module load Anaconda2
conda create -n py3 python=3
source activate py3
anaconda show bioconda/snakemake
conda install --channel https://conda.anaconda.org/bioconda snakemake
```

### Part XIII. Installing and running GATK4.0 in docker to run locally so that you can try out commands without starting any jobs on the cluster.

[Link to official article](https://software.broadinstitute.org/gatk/documentation/article?id=11090)

1. Sign-up for a docker account, download and install it for your OS

2. Make sure it works by running `docker --version`

3. Pull the latest GATK docker container
```bash
docker pull broadinstitute/gatk # defaults to latest

# you can also specify a version
# docker pull broadinstitute/gatk:4.0.5.1
```

4. Launch the GATK container
```bash
# docker image # shows you all the different containers you have on disk

docker run -it broadinstitute/gatk

# if it works you should see your terminal change to something like this:
# root@ea3a5218f494:/gatk#
```

5. Use the GATK wrapper script to issue commands.  GATK 4.0 has combined GATK and PicardTools into one neat package. The container has the gatk wrapper script all set up and ready to go, so you can now run any GATK or Picard command you want. Note that if you want to run a Picard command, you need to use the new syntax, which follows GATK conventions (-I instead of I= and so on). Let's use --list to list all tools available in this version.
```bash
./gatk --list

# the container also has gatk in the enviornmental path so you can also run
# the gatk wrapper from any directory
gatk --list # is just as valid
```

6. Setup a mounted drive to transfer files to and from within the GATK docker container.  To do this we must first shutdown the currently running container and launch a new one with a new command
```bash
exit # shutdown/detatch container and takes you back to a normal terminal prompt

# Check to see if any other instances are running and shutdown any you forgot about
docker ps # Check to see if any other instances are running
docker rm <container_id>
```

7. Increase the amount of RAM allocated to docker. I increased mine to 10gb (my computer has 16gb), so I can run gatk using 8gb of RAM. To do this, just go into the docker preferences on your task bar and the advanced tab. Change it to an amount you have available.

8. Now launch the GATK container with the directory where all your GATK related files are!
```bash
docker run -v /<path_to_data_on_my_mac>:/<path_to_mount_in_docker_container>/ -it broadinstitute/gatk

# My specific run command
# This command allows me to access `../gatk_pipeline` folder inside the docker container
# through the `my_data` folder i specified
docker run -v /Users/uapinyoyingpb/Notebook/Archive/2018/molecular_diagnostics/gatk_pipeline:/gatk/my_data -it broadinstitute/gatk
```
