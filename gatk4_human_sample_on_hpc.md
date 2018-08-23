## Running a real sample Real sample (UC394.1.bam) using theGATK4 pipeline 
- By: Prech Brian Uapinyoying 
- Created: 08/11/2018
- Last updated: 08/22/2018
- <span style="color:red">**STILL INCOMPLETE**</span>- problems with output from haplotypecaller step

1. Optional - Start tmux and interative mode and run the commands directly or send the full scritps to sbatch.
```bash
tmux new -s ruddle2_sesh

srun --pty --mem-per-cpu=6000 --mem=48000 --cpus-per-task=8 --nodes=1 -p interactive bash

#srun --pty --mem-per-cpu=6000 --mem=60000 --cpus-per-task=10 --nodes=1 -p interactive bash

# Get number of cores
echo $SLURM_CPUS_ON_NODE

# Load the modules
module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b
```

2. Relevant directories
```bash
# Reference
/<path>/<to>/reference/
	Homo_sapiens_assembly19.dict
	Homo_sapiens_assembly19.fasta
	Homo_sapiens_assembly19.fasta.amb
	Homo_sapiens_assembly19.fasta.ann
	Homo_sapiens_assembly19.fasta.bwt
	Homo_sapiens_assembly19.fasta.fai
	Homo_sapiens_assembly19.fasta.pac
	Homo_sapiens_assembly19.fasta.rpac
	Homo_sapiens_assembly19.fasta.sa

# Real sample
/<path>/<to>/data/UC394/UC394.1.bam

# Results folder
/<path>/<to>/results/2018-08-08/

# make all the folders to store our results
cd /<path>/<to>/results/2018-08-08/
mkdir s10_gatherBqsrReports s11_recalibratedBam s12_gatherBamFiles s13_genIntervalList s14_haplotypeCaller s15_gatherVcfs s16_gatherHcBams s1_revertsam s2_markIlluminaAdapters s3_samToFastq s4_bwaAlignment s5_mergeBamAlignments s6_markDuplicates s7_sortsam s8_setNmMdAndUqTags s9_baseRecalibrator bamqc tmp
```

3. s1_RevertSam
    - Check what tags to remove
```bash
# step1
# Check attributes first, might take a while if its a huge bam file
samtools view /<path>/<to>/data/UC394/UC394.1.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | awk '{ if(!x[$1]++) { print }}'

# results = X0 MD PG RG AM NM SM MQ OQ UQ OC OP XT XN 
```
    - Now run the command itself
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" RevertSam \
    -I /<path>/<to>/data/UC394/UC394.1.bam \
    -O /<path>/<to>/results/2018-08-08/s1_revertsam/UC394-1_revertsam.unmapped.bam \
    --SANITIZE true \
    --MAX_DISCARD_FRACTION 0.005 \
    --ATTRIBUTE_TO_CLEAR X0 \
    --ATTRIBUTE_TO_CLEAR MD \
    --ATTRIBUTE_TO_CLEAR PG \
    --ATTRIBUTE_TO_CLEAR RG \
    --ATTRIBUTE_TO_CLEAR AM \
    --ATTRIBUTE_TO_CLEAR NM \
    --ATTRIBUTE_TO_CLEAR SM \
    --ATTRIBUTE_TO_CLEAR MQ \
    --ATTRIBUTE_TO_CLEAR OQ \
    --ATTRIBUTE_TO_CLEAR UQ \
    --ATTRIBUTE_TO_CLEAR OC \
    --ATTRIBUTE_TO_CLEAR OP \
    --ATTRIBUTE_TO_CLEAR XT \
    --ATTRIBUTE_TO_CLEAR XN \
    --SORT_ORDER queryname \
    --RESTORE_ORIGINAL_QUALITIES true \
    --REMOVE_DUPLICATE_INFORMATION true \
    --REMOVE_ALIGNMENT_INFORMATION true \
    --TMP_DIR /<path>/<to>/results/2018-08-08/tmp
```

4. s2_markIlluminaAdapters
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" MarkIlluminaAdapters \
    -I /<path>/<to>/results/2018-08-08/s1_revertsam/UC394-1_revertsam.unmapped.bam \
    -O /<path>/<to>/results/2018-08-08/s2_markIlluminaAdapters/UC394-1_markilluminaadapters.unmapped.bam \
    -M /<path>/<to>/results/2018-08-08/s2_markIlluminaAdapters/UC394-1_markilluminaadapters_metrics.txt \
    --ADAPTERS PAIRED_END \
    --TMP_DIR /<path>/<to>/results/2018-08-08/tmp
```

5. s3_samToFastq
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" SamToFastq \
    -I /<path>/<to>/results/2018-08-08/s2_markIlluminaAdapters/UC394-1_markilluminaadapters.unmapped.bam \
    --FASTQ /<path>/<to>/results/2018-08-08/s3_samToFastq/UC394-1_samtofastq_interleave.fastq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INTERLEAVE true \
    --NON_PF true \
    --TMP_DIR /<path>/<to>/results/2018-08-08/tmp
```

6. s4_bwaAlignment - Alignment with bwa mem
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

bwa mem -M -t 7 \
	-p /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
	/<path>/<to>/results/2018-08-08/s3_samToFastq/UC394-1_samtofastq_interleave.fastq \
	> /<path>/<to>/results/2018-08-08/s4_bwaAlignment/UC394-1_bwa_mem.sam
```

7. s5_mergeBamAlignments - GATK4.0 version + monkol's options
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" MergeBamAlignment \
    -R /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
    --UNMAPPED_BAM /<path>/<to>/results/2018-08-08/s2_markIlluminaAdapters/UC394-1_markilluminaadapters.unmapped.bam \
    --ALIGNED_BAM /<path>/<to>/results/2018-08-08/s4_bwaAlignment/UC394-1_bwa_mem.sam \
    -O /<path>/<to>/results/2018-08-08/s5_mergeBamAlignments/UC394-1_mergebamalignment.bam \
    --CREATE_INDEX true \
    --ADD_MATE_CIGAR true \
    --CLIP_ADAPTERS false \
    --CLIP_OVERLAPPING_READS true \
    --INCLUDE_SECONDARY_ALIGNMENTS true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ATTRIBUTES_TO_RETAIN XS \
    --ATTRIBUTES_TO_RETAIN X0 \
    --VALIDATION_STRINGENCY SILENT \
    --EXPECTED_ORIENTATIONS FR \
    --PAIRED_RUN true \
    --IS_BISULFITE_SEQUENCE false \
    --ALIGNED_READS_ONLY false \
    --SORT_ORDER unsorted \
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true \
    --TMP_DIR /<path>/<to>/results/2018-08-08/tmp
```

8. s6_markDuplicates (skipping `CollectQualityYieldMetrics` & `CollectAlignmentSummaryMetrics`)
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" MarkDuplicates \
    --INPUT /<path>/<to>/results/2018-08-08/s5_mergeBamAlignments/UC394-1_mergebamalignment.bam \
    --OUTPUT /<path>/<to>/results/2018-08-08/s6_markDuplicates/UC394-1_markduplicates.bam \
    --METRICS_FILE /<path>/<to>/results/2018-08-08/s6_markDuplicates/UC394-1_markduplicates_metrics.txt \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER queryname \
    --CREATE_MD5_FILE true \
    --VALIDATION_STRINGENCY SILENT \
    --CREATE_INDEX true \
    --TMP_DIR /<path>/<to>/results/2018-08-08/tmp
```

9. s7_sortsam
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" SortSam \
    -I /<path>/<to>/results/2018-08-08/s6_markDuplicates/UC394-1_markduplicates.bam \
    -O /<path>/<to>/results/2018-08-08/s7_sortsam/UC394-1_markduplicates_sorted.bam \
    --SORT_ORDER coordinate
```

10. s8_setNmMdAndUqTags
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" SetNmMdAndUqTags \
    -R /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
    -I /<path>/<to>/results/2018-08-08/s7_sortsam/UC394-1_markduplicates_sorted.bam \
    -O /<path>/<to>/results/2018-08-08/s8_setNmMdAndUqTags/UC394-1_sorted_fixed.bam
```

11. Index the bam
```bash
samtools index /<path>/<to>/results/2018-08-08/s8_setNmMdAndUqTags/UC394-1_sorted_fixed.bam
```

##IMPORTANT NOTE: It looks like by starting from a Broad/GATK processed file and RevertSam, there looks to be no base recalibration done on the samples. 

12. s9_baseRecalibrator
    - The recal.table files have no useful information and so the next few steps after this one below did not work for me using this sample (e.g. `GatherBQSRReports` and therefore nothing to `ApplyBQSR`) So I will keep notes on these steps for testing later on samples from another institution. In the mean time I will skip to `IntervalListTools`, step 17 and use the `UC394-1_sorted_fixed.bam` from step 11 as input into `HaplotypeCaller`.
```bash
# make sure your feature files (vcfs) are indexed `gatk IndexFeatureFile`
# can't write to monkol's dir so made soft-links to his files here:
# /<path>/<to>/reference/*

gatk --java-options "-Xmx8G" IndexFeatureFile --feature-file Homo_sapiens_assembly19.dbsnp138.vcf
gatk --java-options "-Xmx8G" IndexFeatureFile --feature-file Homo_sapiens_assembly19.known_indels.vcf
gatk --java-options "-Xmx8G" IndexFeatureFile --feature-file Homo_sapiens_assembly19.variantEvalGoldStandard.vcf
```
    - Create shell script to generate job submission scripts for each chromosome(s) `br_script_gen.sh`
```bash
#! /bin/bash

# Split the chromosomes
CHR1="-L 1"
CHR2="-L 2"
CHR3="-L 3"
CHR4="-L 4"
CHR5="-L 5"
CHR6="-L 6"
CHR7="-L 7"
CHR8="-L 8"
CHR9="-L 9"
CHR10="-L 10"
CHR11="-L 11"
CHR12_13="-L 12 -L 13"
CHR14_15="-L 14 -L 15"
CHR16_17="-L 16 -L 17"
CHR18_21="-L 18 -L 19 -L 20 -L 21"
CHR22_X="-L 22 -L X"
CHRY_NC="-L Y -L MT -L GL000207.1 -L GL000226.1 -L GL000229.1 \
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

# Create an array of all of them

CHR_ARRAY=( "${CHR1}" "${CHR2}" "${CHR3}" "${CHR4}" "${CHR5}" "${CHR6}" 
    "${CHR7}" "${CHR8}" "${CHR9}" "${CHR10}" "${CHR11}" "${CHR12_13}" 
    "${CHR14_15}" "${CHR16_17}" "${CHR18_21}" "${CHR22_X}" "${CHRY_NC}" )

CHR_NAMES=( "CHR1" "CHR2" "CHR3" "CHR4" "CHR5" "CHR6" "CHR7" "CHR8" "CHR9" 
    "CHR10" "CHR11" "CHR12_13" "CHR14_15" "CHR16_17" "CHR18_21" "CHR22_X" 
    "CHRY_NC" )

JAVA_OPTS='-Xmx8G'

for ((i=0;i<${#CHR_ARRAY[@]};++i)); do
interval_num=$(($i+1))
printf """#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options '-Xmx8G' BaseRecalibrator \
-R /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
-I /<path>/<to>/results/2018-08-08/s8_setNmMdAndUqTags/UC394-1_sorted_fixed.bam \
-known-sites /<path>/<to>/reference/Homo_sapiens_assembly19.dbsnp138.vcf \
-known-sites /<path>/<to>/reference/Homo_sapiens_assembly19.known_indels.vcf \
-known-sites /<path>/<to>/reference/Homo_sapiens_assembly19.variantEvalGoldStandard.vcf \
-O /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_${CHR_NAMES[i]}_recal_data.table \
${CHR_ARRAY[i]} \
--use-original-qualities
""" > UC394-1_${CHR_NAMES[i]}_submit_script.sh
done
```
    - Now create a bash script to submit all the scripts `s9_submit.sh`
```bash
#! /bin/bash

sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR10_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR11_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR12_13_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR14_15_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR16_17_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR18_21_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR1_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR22_X_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR2_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR3_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR4_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR5_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR6_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR7_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR8_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR9_submit_script.sh
sbatch /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHRY_NC_submit_script.sh
```
    - This should generate a bunch of data table files.  Make sure they all are properly generated. You can check the submit script logs.
13. s10_gatherBqsrReports
```bash
# Generate the list
cd /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/
vim recal_data.list

# Paste below in
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR1_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR2_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR3_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR4_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR5_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR6_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR7_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR8_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR9_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR10_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR11_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR12_13_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR14_15_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR16_17_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR18_21_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHR22_X_recal_data.table
/<path>/<to>/results/2018-08-08/s9_baseRecalibrator/UC394-1_CHRY_NC_recal_data.table
```
    - Create a submit script for it running the command `s10_submit.sh`
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" GatherBQSRReports \
    -I /<path>/<to>/results/2018-08-08/s9_baseRecalibrator/recal_data.list \
    -O /<path>/<to>/results/2018-08-08/s10_gatherBqsrReports/combined_recal_data.table
```
    - Run the submit script at the command line `sbatch s10_submit.sh`

14. s11_recalibratedBam - `ApplyBQSR` which replaced PrintReads step in gatk3.4
```bash
#! /bin/bash

# Split the chromosomes
CHR1="-L 1"
CHR2="-L 2"
CHR3="-L 3"
CHR4="-L 4"
CHR5="-L 5"
CHR6="-L 6"
CHR7="-L 7"
CHR8="-L 8"
CHR9="-L 9"
CHR10="-L 10"
CHR11="-L 11"
CHR12_13="-L 12 -L 13"
CHR14_15="-L 14 -L 15"
CHR16_17="-L 16 -L 17"
CHR18_21="-L 18 -L 19 -L 20 -L 21"
CHR22_X="-L 22 -L X"
CHRY_NC="-L Y -L MT -L GL000207.1 -L GL000226.1 -L GL000229.1 \
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

# Create an array of all of them

CHR_ARRAY=( "${CHR1}" "${CHR2}" "${CHR3}" "${CHR4}" "${CHR5}" "${CHR6}" 
    "${CHR7}" "${CHR8}" "${CHR9}" "${CHR10}" "${CHR11}" "${CHR12_13}" 
    "${CHR14_15}" "${CHR16_17}" "${CHR18_21}" "${CHR22_X}" "${CHRY_NC}" )

CHR_NAMES=( "CHR1" "CHR2" "CHR3" "CHR4" "CHR5" "CHR6" "CHR7" "CHR8" "CHR9" 
    "CHR10" "CHR11" "CHR12_13" "CHR14_15" "CHR16_17" "CHR18_21" "CHR22_X" 
    "CHRY_NC" )

JAVA_OPTS='-Xmx8G'

for ((i=0;i<${#CHR_ARRAY[@]};++i)); do
interval_num=$(($i+1))
printf """#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options '-Xmx8G' ApplyBQSR \
    -R /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
    -I /<path>/<to>/results/2018-08-08/s8_setNmMdAndUqTags/UC394-1_sorted_fixed.bam \
    ${CHR_ARRAY[i]} \
    --bqsr-recal-file /<path>/<to>/results/2018-08-08/s10_gatherBqsrReports/combined_recal_data.table \
    -O /<path>/<to>/results/2018-08-08/s11_recalibratedBam/UC394-1_${CHR_NAMES[i]}_recal_reads.bam
""" > UC394-1_${CHR_NAMES[i]}_applybqsr_submit_script.sh
done
```
    - Create a script to submit them all through `sbatch`

16. s12_gatherBamFiles
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" GatherBamFiles \
    -I /<path>/<to>/results/2018-08-08/s11_recalibratedBam/recal_bam.list \
    -O /<path>/<to>/results/2018-08-08/s12_gatherBamFiles/UC394-1_gathered_recal.bam

# Index the bam file before running haplotype caller
samtools index /<path>/<to>/results/2018-08-08/s12_gatherBamFiles/UC394-1_gathered_recal.bam
```

17. s13_genIntervalList - generate 10 interval lists for calling exome samples using `IntervalListTools`
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" IntervalListTools \
    --PADDING 0 \
    --SCATTER_COUNT 10 \
    --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --UNIQUE true \
    --SORT true \
    --BREAK_BANDS_AT_MULTIPLES_OF 0 \
    --INPUT /<path>/<to>/reference/exome_calling_regions.v1.interval_list \
    --OUTPUT /<path>/<to>/results/2018-08-08/s13_genIntervalList/ # This needs to be a directory not a file
```

18. s14_haplotypeCaller - Run haplotype caller on all 10 scattered intervals
    - Generate scripts for running the sample over each interval
    - Will use the fixed tags bam from s8 folder for this test since BQSR did nothing
    - Use this in future: `# -I /<path>/<to>/results/2018-08-08/s12_gatherBamFiles/UC394-1_gathered_recal.bam \`
    - `hc_script_gen.sh`
```bash
#! /bin/bash

INTERVAL_ROOT='/<path>/<to>/results/2018-08-08/s13_genIntervalList'
OUTPUT_ROOT='/<path>/<to>/results/2018-08-08/s14_haplotypeCaller'

# Create an array of the 10 temp directories produced by the scatter step
TEMP_DIRS=( 'temp_0001_of_10' 'temp_0002_of_10' 'temp_0003_of_10' 
	'temp_0004_of_10' 'temp_0005_of_10' 'temp_0006_of_10' 'temp_0007_of_10' 
	'temp_0008_of_10' 'temp_0009_of_10' 'temp_0010_of_10' )

JAVA_OPTS='-Xmx8G'

for ((i=0;i<${#TEMP_DIRS[@]};++i)); do
interval_num=$(($i+1))
printf """#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options ${JAVA_OPTS} HaplotypeCaller \
-L ${INTERVAL_ROOT}/${TEMP_DIRS[i]}/scattered.interval_list \
-I /<path>/<to>/results/2018-08-08/s8_setNmMdAndUqTags/UC394-1_sorted_fixed.bam \
-O ${OUTPUT_ROOT}/UC394-1_interval_${interval_num}_output.raw.snps.indels.vcf \
-R /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
-ERC GVCF \
--min-pruning 3 \
--max-alternate-alleles 3 \
--max-num-haplotypes-in-population 200 \
-contamination 0.0 \
--read-filter OverclippedReadFilter \
-bamout UC394-1_interval_${interval_num}_local_realigned.bam
""" > UC394-1_interval_${interval_num}_hc_submit_script.sh
done
```
    - Generate script
```
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_10_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_1_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_2_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_3_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_4_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_5_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_6_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_7_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_8_hc_submit_script.sh
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_9_hc_submit_script.sh
```

19. s15_gatherVCFs
```bash
ls -1 /<path>/<to>/results/2018-08-08/s14_haplotypeCaller/*.vcf > /<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_vcf.list
### IMPORTANT make sure its in the correct order
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_1_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_2_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_3_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_4_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_5_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_6_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_7_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_8_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_9_output.raw.snps.indels.vcf
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_10_output.raw.snps.indels.vcf
```
    - running the script
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b
gatk --java-options "-Xmx8G" GatherVcfs \
    --INPUT /<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_vcf.list \
    --OUTPUT /<path>/<to>/results/2018-08-08/s15_gatherVcfs/UC394-1.vcf.gz
```

20. s16_gatherHcBams - gather the bams again, but this time ones locally realigned and output from haplotype caller
```bash
# /<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_bam.list
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_1_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_2_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_3_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_4_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_5_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_6_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_7_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_8_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_9_local_realigned.bam
/<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_interval_10_local_realigned.bam

gatk --java-options "-Xmx8G" GatherBamFiles \
    -I /<path>/<to>/results/2018-08-08/s14_haplotypeCaller/UC394-1_bam.list \
    -O /<path>/<to>/results/2018-08-08/s16_gatherHcBams/UC394-1_gathered_hc.bam

# Index the bam file before running haplotype caller
samtools index /<path>/<to>/results/2018-08-08/s16_gatherHcBams/UC394-1_gathered_hc.bam
```

21. bamqc - Had lots of problems, going to try QC metrics to troubleshoot 
    - `CollectHsMetrics`
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" CollectHsMetrics \
    --REFERENCE_SEQUENCE /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
    --INPUT /<path>/<to>/results/2018-08-08/s8_setNmMdAndUqTags/UC394-1_sorted_fixed.bam \
    --OUTPUT /<path>/<to>/results/2018-08-08/bamqc/UC394-1_sorted_fixed.hybrid_selection_metrics \
    --TARGET_INTERVALS /<path>/<to>/reference/exome_calling_regions.v1.interval_list \
    --BAIT_INTERVALS /<path>/<to>/reference/exome_calling_regions.v1.interval_list \
    --TMP_DIR /<path>/<to>/results/2018-08-08/tmp
```
    - CollectHsMetrics
```bash
#! /bin/bash

# Specify partition
#SBATCH -p general
 
# Specify number of nodes 
#SBATCH --nodes=1
 
# Specify number of cores
#SBATCH --cpus-per-task=8

# Total memory
#SBATCH --mem=48000 
#SBATCH --mem-per-cpu=6000
 
# 6 hour timelimit:
#SBATCH --time 0-06:00:00

module load GATK/4.0.6.0-Java-1.8.0_121
module load SAMtools/1.9-foss-2016b
module load BWA/0.7.17-foss-2016b

gatk --java-options "-Xmx8G" CollectMultipleMetrics \
    --REFERENCE_SEQUENCE /<path>/<to>/reference/Homo_sapiens_assembly19.fasta \
    --INPUT /<path>/<to>/results/2018-08-08/s8_setNmMdAndUqTags/UC394-1_sorted_fixed.bam \
    --OUTPUT /<path>/<to>/results/2018-08-08/bamqc/UC394-1_sorted_fixed \
    --PROGRAM null \
    --PROGRAM MeanQualityByCycle \
    --PROGRAM QualityScoreDistribution \
    --PROGRAM CollectInsertSizeMetrics \
    --PROGRAM CollectAlignmentSummaryMetrics \
    --METRIC_ACCUMULATION_LEVEL null \
    --METRIC_ACCUMULATION_LEVEL ALL_READS \
    --ASSUME_SORTED true \
    --TMP_DIR /<path>/<to>/results/2018-08-08/tmp
```