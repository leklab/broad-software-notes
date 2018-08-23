## GATK4 example using Cromwell

This short example was done in the docker version of GATK4. [See detailed notes on GATK4](gatk4_pipeline_detailed_notes.md)

### Part I. A toy example of how cromwell works

1. Running a simple cromwell workflow
	- The `myWorkflow.wdl`
```bash
workflow myWorkflow {

	call echo_string
	call cat_file {
		input:
			echo_file = echo_string.echo_file
	}
}

task echo_string {

	String out_string
	String out_filename
	command {
		echo ${out_string} > ${out_filename}.txt
	}
	output {
		#String out=read_string(stdout())
		File echo_file = "${out_filename}.txt"
	}
}

task cat_file {
	File echo_file
	command {
		cat ${echo_file}
	}
	output {
		String out=read_string(stdout())
	}

}
```

2. Setting up a config file (e.g. your.conf)
	- Add this to the top of the config file: `include required(classpath("application"))`
	- The default Cromwell configuration values are set via Cromwell's `application.conf`. To ensure that you always have the defaults from the application.conf, you must include it at the top of your new configuration file.
	- Now run the workflow again with the cromwell config file
```bash
java -jar -Dconfig.file=your.conf $CROM_JAR run myWorkflow.wdl
```

3. Validate your wdl script using wdltool.jar and generate template json file for input variables
	- [Quickstart tutorial page](https://software.broadinstitute.org/wdl/documentation/quickstart)
```bash
# Check for syntax errors
# $WOM_JAR = wdltool.jar
java -jar $WOM_JAR validate myWorkflow.wdl

# create template JSON file for workflow inputs
java -jar $WOM_JAR inputs myWorkflow.wdl > myWorkflow_inputs.json
```
	- fill out the JSON input template file with the literals
```json
{
  "myWorkflow.echo_string.out_string": "Monkol's Minions!",
  "myWorkflow.echo_string.out_filename": "testy"
}
```
4. Run with all inputs
```bash
java -jar -Dconfig.file=your.conf $CROM_JAR run myWorkflow.wdl --inputs myWorkflow_inputs.json
```

5.  The output is organized into structured directories, and its crazy deep. May need another task just to move the final result files to a specified folder.
```bash
/gatk/my_data/sandbox/cromwell-executions/myWorkflow/<hash_num>/<task_name>/execution/
/gatk/my_data/sandbox/cromwell-executions/myWorkflow/<hash_num>/<task_name>/input/<hash_num>/

# Actual example
/gatk/my_data/sandbox/cromwell-executions/myWorkflow/40b54d81-94e8-4a1c-adbc-f222765655e5/call-cat_file/execution
/gatk/my_data/sandbox/cromwell-executions/myWorkflow/40b54d81-94e8-4a1c-adbc-f222765655e5/call-cat_file/inputs/664275994/testy.txt
```

### Part II. Cromwell test using revertsam in docker

1. Important files and directories for the script.
	- Run WDL scripts here:
		`/gatk/my_data/results/2018-08-05_gatk4_cromwell`
	- Ref genome:
		```/gatk/my_data/references/hg19_decoy
				hs37d5.dict
				hs37d5.fa
				hs37d5.fa.amb
				hs37d5.fa.ann
				hs37d5.fa.bwt
				hs37d5.fa.fai
				hs37d5.fa.pac
				hs37d5.fa.sa
		```
	- Input Data:
		`/gatk/my_data/test_data/tutorial_6484_RevertSam/6484_snippet.bam`

2. Original `RevertSam` command for gatk4
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
```

3. Convered to wdl
```bash
workflow dockerGatk4 {
	call RevertSam
}

task RevertSam {
	File gatk
	File orig_aln_bam
	String java_opts
	String sample_name
	command {
		gatk --java-options $java_opts RevertSam \
		    -I ${orig_aln_bam} \
		    -O ${sample_name}_revertsam.unmapped.bam \
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
	}
	output {
		#String out=read_string(stdout())
		File revertsam_ubam = "${sample_name}_revertsam.unmapped.bam"
	}
}
```

4. Validate WDL
```bash
java -jar $WOM_JAR validate dockerGatk4.wdl
```

5. Generate input template json and fill it out
	- Generating the input file
```bash
java -jar $WOM_JAR inputs dockerGatk4.wdl > dockerGatk4_inputs.json
```
	- Original template
```json
{
  "dockerGatk4.RevertSam.gatk_path": "String",
  "dockerGatk4.RevertSam.orig_aln_bam": "File",
  "dockerGatk4.RevertSam.sample_name": "String",
  "dockerGatk4.RevertSam.java_opts": "String"
}
```
	- Filling out the variables
```json
{
  "dockerGatk4.RevertSam.gatk_path": "/gatk/gatk",
  "dockerGatk4.RevertSam.sample_name": "6484_snippet",
  "dockerGatk4.RevertSam.java_opts": "-Xmx8G",
  "dockerGatk4.RevertSam.orig_aln_bam": "/gatk/my_data/test_data/tutorial_6484_RevertSam/6484_snippet.bam"
}
```

6. Run the command to see if it will work!
```bash
cd /gatk/my_data/results/2018-08-05_gatk4_cromwell

java -jar -Dconfig.file=dockerGatk4.conf $CROM_JAR run dockerGatk4.wdl --inputs dockerGatk4_inputs.json
```
	- Looks like it works, but because our MAX_DISCARD_FRACTION is pretty low, RevertSam produced and error code 3 upon exit.  This was interpreted as a failure because we did not setup the pipeline to handle error code 3 (0 means success, 1 means failure) even though it is only an informational setting.  We could set the MAX_DISCARD_FRACTION to the default of 0.01 or change the `ContinueOnReturnCode` option.

7. Changing the default return codes so Cromwell doesn't quit on harmless exit codes
	- Exit codes can be found in the RC file:
```bash
/gatk/my_data/results/2018-08-05_gatk4_cromwell/cromwell-executions/dockerGatk4/<hash_num>/call-RevertSam/execution/rc
```
	- [Documentation on setting the runtime attributes](http://cromwell.readthedocs.io/en/develop/wf_options/Overview/#setting-default-runtime-attributes)
	- Create a `dockerGatk4_options.json` and add the acceptable return codes
```bash
{
    "default_runtime_attributes": {
        "continueOnReturnCode": [0, 3]
    }
}
```

8. Rerun with the supplied `options.json` file
```bash
cd /gatk/my_data/results/2018-08-05_gatk4_cromwell

java -jar -Dconfig.file=dockerGatk4.conf $CROM_JAR run dockerGatk4.wdl --inputs dockerGatk4_inputs.json --options dockerGatk4_options.json
```

9. Lets continue the pipeline by adding the next task: MarkIlluminaAdapters
	- Original command
```bash
gatk --java-options "-Xmx8G" MarkIlluminaAdapters \
    -I /gatk/my_data/results/revertsam/6484_snippet_revertsam.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    -M /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters_metrics.txt \
    --ADAPTERS PAIRED_END
    # --TMP_DIR <path>/tmp
```
	- new draft of WDL script. Reorganized to include workflow (global) level variables. Many variable names have changed in spelling and notice there are now input variables in each `call`
```bash
workflow dockerGatk4 {
    # workflow (global) variables
    File gatk_path
    String name

    call RevertSam {
        input:
            GATK = gatk_path,
            sampleName = name
    }

    call MarkIlluminaAdapters {
        input:
            GATK = gatk_path,
            sampleName = name,
            revertSam_uBam = RevertSam.revertSam_uBam
    }
}

task RevertSam {
    # Task (local) level variables
    File GATK
    File origAlnBam
    String javaOpts
    String sampleName
    command {
        ${GATK} --java-options ${javaOpts} RevertSam \
            -I ${origAlnBam} \
            -O ${sampleName}_revertsam.unmapped.bam \
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
    }
    output {
        File revertSam_uBam = "${sampleName}_revertsam.unmapped.bam"
    }
}
task MarkIlluminaAdapters {
    File GATK
    File revertSam_uBam
    String sampleName
    String javaOpts
    command {
        ${GATK} --java-options ${javaOpts} MarkIlluminaAdapters \
        -I ${revertSam_uBam} \
        -O ${sampleName}_markilluminaadapters.unmapped.bam \
        -M ${sampleName}_markilluminaadapters_metrics.txt \
        --ADAPTERS PAIRED_END
    }
    output {
        File markIlluminaAdapters_uBam = "${sampleName}_markilluminaadapters.unmapped.bam"
        File markIlluminaAdapters_metrics = "${sampleName}_markilluminaadapters_metrics.txt"
    }
}
```
	- New `input.json` file
```bash
{
  "dockerGatk4.RevertSam.origAlnBam": "File",
  "dockerGatk4.MarkIlluminaAdapters.javaOpts": "String",
  "dockerGatk4.gatk_path": "File",
  "dockerGatk4.RevertSam.javaOpts": "String",
  "dockerGatk4.name": "String"
}
```
	- Organize it
```bash
{
  "##_COMMENT1": "WORKFLOW LEVEL VARIABLES",
  "dockerGatk4.name": "String",
  "dockerGatk4.gatk_path": "File",

  "##_COMMENT2": "REVERTSAM VARIABLES",
  "dockerGatk4.RevertSam.origAlnBam": "File",
  "dockerGatk4.RevertSam.javaOpts": "String",
  
  "##_COMMENT3": "MARKILLUMINAADAPTERS VARIABLES",
  "dockerGatk4.MarkIlluminaAdapters.javaOpts": "String",
}
```
	- Fill it out
```bash
{
  "##_COMMENT1": "WORKFLOW LEVEL VARIABLES",
  "dockerGatk4.name": "6484_snippet",
  "dockerGatk4.gatk_path": "/gatk/gatk",

  "##_COMMENT2": "REVERTSAM VARIABLES",
  "dockerGatk4.RevertSam.origAlnBam": "/gatk/my_data/test_data/tutorial_6484_RevertSam/6484_snippet.bam",
  "dockerGatk4.RevertSam.javaOpts": "-Xmx8G",
  
  "##_COMMENT3": "MARKILLUMINAADAPTERS VARIABLES",
  "dockerGatk4.MarkIlluminaAdapters.javaOpts": "-Xmx8G"
}
```

10. Ok, lets try to add the rest of the steps that we can perform in the docker container (stop before Base Recalibration step)
```bash
gatk --java-options "-Xmx8G" SamToFastq \
    -I /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    --FASTQ /gatk/my_data/results/revertsam/6484_snippet_samtofastq_interleave.fastq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INTERLEAVE true \
    --NON_PF true
    # --TMP_DIR <path>/tmp

bwa mem -M -t 7 -p /gatk/my_data/references/hg19_decoy/hs37d5.fa /gatk/my_data/results/revertsam/6484_snippet_samtofastq_interleave.fastq > 6484_snippet_bwa_mem.sam

gatk --java-options "-Xmx8G" MergeBamAlignment \
    -R /gatk/my_data/references/hg19_decoy/hs37d5.fa \
    --UNMAPPED_BAM /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    --ALIGNED_BAM /gatk/my_data/results/revertsam/6484_snippet_bwa_mem.sam \
    -O /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment.bam \
    --CREATE_INDEX true \
    --ADD_MATE_CIGAR true \
    --CLIP_ADAPTERS false \
    --CLIP_OVERLAPPING_READS=true \
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
    --UNMAPPED_READ_STRATEGY COPY_TO_TAG
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --UNMAP_CONTAMINANT_READS true  
    #--TMP_DIR=/<path>/tmp # optional to process large files

gatk --java-options "-Xmx8G" CollectQualityYieldMetrics \
    -I /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_markilluminaadapters_qualityyieldmetrics.txt

gatk --java-options "-Xmx8G" CollectAlignmentSummaryMetrics \
    -R /gatk/my_data/references/hg19_decoy/hs37d5.fa \
    -I /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_mergebamalignment_alignmentsummarymetrics.txt 

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

gatk --java-options "-Xmx8G" SortSam \
    -I /gatk/my_data/results/revertsam/6484_snippet_markduplicates.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_markduplicates_sorted.bam \
    --SORT_ORDER coordinate

gatk --java-options "-Xmx8G" SetNmMdAndUqTags \
    -R /gatk/my_data/references/hg19_decoy/hs37d5.fa \
    -I /gatk/my_data/results/revertsam/6484_snippet_markduplicates_sorted.bam \
    -O /gatk/my_data/results/revertsam/6484_snippet_sorted_fixed.bam
```
	- Json input
```bash
{
  "dockerGatk4.ref": "/gatk/my_data/references/hg19_decoy/hs37d5.fa",
  "dockerGatk4.gatk_path": "/gatk/gatk",  
  "dockerGatk4.bwa_path": "/gatk/my_data/bwa-0.7.17/bwa",
  "dockerGatk4.java_options": "-Xmx8G",
  "dockerGatk4.name": "6484_snippet",
  "dockerGatk4.RevertSam.origAlnBam": "/gatk/my_data/test_data/tutorial_6484_RevertSam/6484_snippet.bam",
  "dockerGatk4.BwaAlignment.threads": "7",
}
```

11. Rerun with the supplied options.json file
```bash
cd /gatk/my_data/results/2018-08-05_gatk4_cromwell

java -jar -Dconfig.file=dockerGatk4.conf $CROM_JAR run dockerGatk4.wdl --inputs dockerGatk4_inputs.json --options dockerGatk4_options.json
```

12. I ran into a bug. So the example stops here.

  - Docker can't make hard-links to link to my reference data and soft-links are no good inside docker for whatever reason. So I narrowed down the issue to likely a permissions problem. I been running Docker as root, so maybe cromwell cannot access some of these directories because of that.

  - Turns out creating a new user for the docker container to avoid using the root user requires modifying the original GATK4 Dockerfile, a huge pain. So I will just stop the cromwell example here. I will implement the the whole thing in on the cluster later

  - See [the hard-link issue discussed here](https://gatkforums.broadinstitute.org/wdl/discussion/10347/localization-via-hard-link-has-failed) and [docker and ROOT issue here](https://medium.com/@mccode/processes-in-containers-should-not-run-as-root-2feae3f0df3b)
