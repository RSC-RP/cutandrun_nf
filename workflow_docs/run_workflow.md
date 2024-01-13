Run Cut&Run Peak Calling
================

-   <a href="#about-the-pipeline" id="toc-about-the-pipeline">About the
    Pipeline</a>
-   <a href="#activate-the-environment-on-hpc"
    id="toc-activate-the-environment-on-hpc">Activate the Environment on
    HPC</a>
    -   <a href="#1-interactive-session" id="toc-1-interactive-session">1)
        Interactive Session</a>
    -   <a href="#2-open-cutandrun-workflow-folder"
        id="toc-2-open-cutandrun-workflow-folder">2) Open CutandRun workflow
        folder</a>
    -   <a href="#3-activate-conda-environement"
        id="toc-3-activate-conda-environement">3) Activate conda
        environement</a>
-   <a href="#examine-the-sample-sheet"
    id="toc-examine-the-sample-sheet">Examine the Sample Sheet</a>
    -   <a href="#examples" id="toc-examples">Examples</a>
-   <a href="#run-the-example-data" id="toc-run-the-example-data">Run the
    Example Data</a>
-   <a href="#configure-pipeline-for-your-data"
    id="toc-configure-pipeline-for-your-data">Configure Pipeline for Your
    Data</a>
    -   <a href="#configuration-file" id="toc-configuration-file">Configuration
        file</a>
    -   <a href="#global-params" id="toc-global-params">Global Params</a>
    -   <a href="#genomic-references" id="toc-genomic-references">Genomic
        References</a>
    -   <a href="#seacr" id="toc-seacr">SEACR</a>
    -   <a href="#optional-macs2" id="toc-optional-macs2">Optional: MACS2</a>
    -   <a href="#advanced-options" id="toc-advanced-options">Advanced
        Options</a>
-   <a href="#run-script" id="toc-run-script">Run Script</a>
    -   <a href="#alignment-and-peak-calls"
        id="toc-alignment-and-peak-calls">Alignment and Peak Calls</a>
    -   <a href="#optional-build-the-index-and-exit-pipeline"
        id="toc-optional-build-the-index-and-exit-pipeline">Optional: Build the
        Index and Exit Pipeline</a>
-   <a href="#expected-outputs" id="toc-expected-outputs">Expected
    Outputs</a>
    -   <a href="#final-outputs" id="toc-final-outputs">Final Outputs</a>
    -   <a href="#file-structure" id="toc-file-structure">File Structure</a>
    -   <a href="#detailed-file-structure"
        id="toc-detailed-file-structure">Detailed File Structure</a>
    -   <a href="#pipeline-reports" id="toc-pipeline-reports">Pipeline
        Reports</a>

# About the Pipeline

The pipeline runs the
[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
alignment, quality trimming of reads with trimgalore,
[SEACR](https://github.com/FredHutch/SEACR) peak calling, and optionally
[MACS2](https://github.com/macs3-project/MACS) peak calling. MACS2
requires an effective genome size to call peaks, which you can provide
directly or call
[`unique-kmers.py`](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html)
to calculate the effective genome size on the fly. Coverage tracks are
produced for visualization in [IGV](https://igv.org/doc/desktop/).

It will also perform general QC statistics on the fastqs with
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/),
the alignment, peak calling, and sample similarity using
[deeptools](https://deeptools.readthedocs.io/en/develop/). Finally, the
QC reports are collected into a single file using
[multiQC](https://multiqc.info/).

A DAG (directed acyclic graph) of the default workflow is show below:

<img src="images/dag.png" width="4120" style="display: block; margin: auto;" />

# Activate the Environment on HPC

The directions to set-up the Nextflow workflow requirements are found in
the [README.md](../README.md). Ensure that you have followed the steps
to fork and clone the repository and created the conda nextflow
environment before starting with this document.

### 1) Interactive Session

Optional but recommended: use `tmux` on the cybertron login nodes. Name
the session nextflow and then request an interactive session, then
activate the nextflow conda environment. The project codes can be found
with `project info` command. Change the `$QUEUE` and `$NAME` variables
in the code chunk below to be accurate for your Cybertron projects.

``` bash
tmux new-session -s nextflow
project info
NAME="RSC_adhoc"
QUEUE="sceaq"
qsub -I -q $QUEUE -P $(project code $NAME) -l select=1:ncpus=1:mem=8g -l walltime=8:00:00
```

### 2) Open CutandRun workflow folder

Navigate to where you place the cloned (copied) cutandrun_nf directory,
and then checkout the latest release branch.

``` bash
cd /path/to/my/cutandrun_nf
git fetch
# this will list all branches of the repository. Find the release branch with the latest version in red text, which means you don't yet have a local copy of that branch. 
git branch -a
# then select the latest version, for example 2.0.0. This downloads the stable version of pipeline locally. 
git checkout release/2.0.0
# Now the * indicates that you're on the release branch and its no longer red text. 
git branch
```

### 3) Activate conda environement

Activate the Nextflow conda environment.

``` bash
conda env create -f env/nextflow.yaml
conda activate nextflow
```

# Examine the Sample Sheet

A sample sheet in csv (comma separated values) format is used as input
to the pipeline. This sample sheet **must** have the following column
names in any order:

-   “sample”
-   “sample_id”
-   “target_or_control”
-   “single_end”
-   “read1”
-   “read2”

| column_name       | column_description                                                                                                                                                                                                                                                                                                                                                                                |
|:------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| sample            | Any alphanumeric string for each biological sample in the dataset. Will have the same sample IDs for each antibody used. For example SAMPLE_1 has both H3K27me3 and IgG control CUT&RUN, and thus SAMPLE_1 has 1 row with the files for H3K27me3, and SAMPLE_1 has 2nd row with the files for IgG data.                                                                                           |
| sample_id         | Any alphanumeric string for each unique sample+condition. No duplicates allowed. For example SAMPLE_1 has both H3K27me3 and IgG control CUT&RUN. Thus, SAMPLE_1 is the value in `sample`, and “SAMPLE_1\_H3K27me3” is the value in `sample_id`. Again, SAMPLE_1 has 2nd row with the files for IgG data, where SAMPLE_1 is the value in `sample`, and “SAMPLE_1\_IgG” is the value in `sample_id` |
| target_or_control | Must contain the values \[target or control\] case-sensitive. Target is for the antibodies using the immunoprecipitation for the proteins of interest, such as transcription factors or histone modifications like H3K27me3, or the value control for the isotype control (eg IgG).                                                                                                               |
| read1             | Contain absolute filepaths to read 1 in paired-end fastqs.                                                                                                                                                                                                                                                                                                                                        |
| read2             | Contain absolute filepaths to read 2 in paired-end fastqs.                                                                                                                                                                                                                                                                                                                                        |
| single_end        | For CUT&RUN data it should always be \[false\] case-sensitive.                                                                                                                                                                                                                                                                                                                                    |

### Examples

**1)** Below is an example of a complete sample sheet for use in the
pipeline, which can be edited for your own samples in
`test_data/test_dataset_sample_sheet.csv`.

-   It contains IgG control samples for peak calling.
-   This sample sheets OK to use even if you elect to skip IgG
    normalization in SEACR or use IgG background in MACS2. The pipeline
    will simply not use the controls.
-   Use `threshold`, and `no_control_macs2` parameters in
    `nextflow.config` to change this. Details found in [Configure
    Pipeline for Your Data](#configure-pipeline-for-your-data)

| sample | sample_id   | single_end | target_or_control | read1                                                                                         | read2                                                                                         |
|:-------|:------------|:-----------|:------------------|:----------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------|
| M1     | M1_H3K27_NK | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K27_NK_chr17_R1_ecoli.fastq.gz | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K27_NK_chr17_R2_ecoli.fastq.gz |
| M1     | M1_H3K4_NK  | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K4_NK_chr17_R1_ecoli.fastq.gz  | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K4_NK_chr17_R2_ecoli.fastq.gz  |
| M1     | M1_IgG_NK   | false      | control           | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_IgG_NK_chr17_R1_ecoli.fastq.gz   | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_IgG_NK_chr17_R2_ecoli.fastq.gz   |
| M2     | M2_H3K27_NK | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K27_NK_chr17_R1_ecoli.fastq.gz | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K27_NK_chr17_R2_ecoli.fastq.gz |
| M2     | M2_H3K4_NK  | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K4_NK_chr17_R1_ecoli.fastq.gz  | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K4_NK_chr17_R2_ecoli.fastq.gz  |
| M2     | M2_IgG_NK   | false      | control           | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_IgG_NK_chr17_R1_ecoli.fastq.gz   | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_IgG_NK_chr17_R2_ecoli.fastq.gz   |

**2)** Below is another example of a complete sample sheet for use in
the pipeline.

-   It lacks IgG control samples for peak calling.
-   This sample sheets OK to use only if you modify the parameters to
    skip using IgG controls.
-   Use `threshold`, and `no_control_macs2` parameters in
    `nextflow.config` to modify this. Details found in [Configure
    Pipeline for Your Data](#configure-pipeline-for-your-data).

| sample | sample_id   | single_end | target_or_control | read1                                                                                         | read2                                                                                         |
|:-------|:------------|:-----------|:------------------|:----------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------|
| M1     | M1_H3K27_NK | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K27_NK_chr17_R1_ecoli.fastq.gz | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K27_NK_chr17_R2_ecoli.fastq.gz |
| M1     | M1_H3K4_NK  | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K4_NK_chr17_R1_ecoli.fastq.gz  | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_H3K4_NK_chr17_R2_ecoli.fastq.gz  |
| M1     | M1_IgG_NK   | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_IgG_NK_chr17_R1_ecoli.fastq.gz   | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M1_IgG_NK_chr17_R2_ecoli.fastq.gz   |
| M2     | M2_H3K27_NK | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K27_NK_chr17_R1_ecoli.fastq.gz | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K27_NK_chr17_R2_ecoli.fastq.gz |
| M2     | M2_H3K4_NK  | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K4_NK_chr17_R1_ecoli.fastq.gz  | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_H3K4_NK_chr17_R2_ecoli.fastq.gz  |
| M2     | M2_IgG_NK   | false      | target            | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_IgG_NK_chr17_R1_ecoli.fastq.gz   | /gpfs/shared_data/demo_data/mus_musculus/cutandrun/fastqs/M2_IgG_NK_chr17_R2_ecoli.fastq.gz   |

# Run the Example Data

To ensure that the pipeline works, first run the test data set. This
example will run using the data found in the `test_sample_sheet.csv`.

``` bash
./main_run.sh "test_dataset"
```

# Configure Pipeline for Your Data

## Configuration file

Open the configuration file `nextflow.config` and edit the necessary
parameters for building the index, and/or running the alignment or peak
calling steps.

    ## //working directory for temporary/intermediate files produced in the workflow processes
    ## workDir = "$HOME/temp"
    ## 
    ## //global parameters
    ## params {
    ##     // general options
    ##     sample_sheet                = "./test_data/test_dataset_sample_sheet.csv"
    ##     queue                       = 'paidq'
    ##     project                     = '207f23bf-acb6-4835-8bfe-142436acb58c'
    ##     outdir                      = "./results/mouse"
    ##     peaks_outdir                = "${params.outdir}/macs_no_igg"
    ##     publish_dir_mode            = 'copy'
    ## 
    ##     //Bowtie params for target genome
    ##     build_index                 = false
    ##     fasta                       = '/gpfs/shared_data/Bowtie2/mm39.fa' // required
    ##     index                       = '/gpfs/shared_data/Bowtie2/mm39_index/' // bowtie2 index path is required unless `build_index = true`
    ##     save_unaligned              = false
    ## 
    ##     // Bowtie params for spike-in genome
    ## <...>

## Global Params

Be sure to change the following lines for the global parameters:

-   sample_sheet
-   queue
-   project code
-   outdir
-   peak_outdir

<!-- -->

    ## Warning in params_lines[1]:end: numerical expression has 2 elements: only the
    ## first used

    ## //global parameters
    ## params {
    ##     // general options
    ##     sample_sheet                = "./test_data/test_dataset_sample_sheet.csv"
    ##     queue                       = 'paidq'
    ##     project                     = '207f23bf-acb6-4835-8bfe-142436acb58c'
    ##     outdir                      = "./results/mouse"
    ##     peaks_outdir                = "${params.outdir}/macs_no_igg"
    ##     publish_dir_mode            = 'copy'

## Genomic References

Additionally, determine if you require a new bowtie2 index to be build
for the target genome and/or the spike-in genome. The pipeline requires
either a fasta filepath OR Bowtie2 index filepath. This is also required
for the spike-in, with E. Coli provided as a default.

E. coli is the default since it that is a carry over DNA from the
Cut&Run library prep methodology and is expected to be present in all
Cut&Run experiments regardless if exogenous spike-in is used like Yeast.
Please see [here](https://doi.org/10.7554/eLife.46314) for more
information on spike-in normalization.

Change the following lines for alignment reference files when needed:

-   build_index
-   fasta
-   index
-   build_spike_index
-   spike_fasta
-   spike_index

<!-- -->

    ##     //Bowtie params for target genome
    ##     build_index                 = false
    ##     fasta                       = '/gpfs/shared_data/Bowtie2/mm39.fa' // required
    ##     index                       = '/gpfs/shared_data/Bowtie2/mm39_index/' // bowtie2 index path is required unless `build_index = true`
    ##     save_unaligned              = false
    ## 
    ##     // Bowtie params for spike-in genome

## SEACR

SEACR defaults to using IgG control normalization and stringent peak
calling, eg `SEACR_1.3.sh target_bedgraph igg_bedgraph norm stringent`.

If skipping the use of IgG control all together, set `threshold` to any
value \> 0.

If you would like to use a spike-in normalization, either E. Coli or an
exongenous spike-in like Drosophila, set `spike_norm` to true. You must
see [Advanced Options](#advanced-options) for details on turning off igg
normalization.

-   threshold
-   spike_norm
-   chrom_sizes
-   scale_factor_constant

``` r
c(config_file[grep("SEACR params", config_file):(grep("scale_factor_constant", config_file) +
    1)]) %>%
    cat(., sep = "\n")
```

## Optional: MACS2

Finally, decide whether to run MACS2 calls along with the SEACR peak
calling algorithm (default = true). MACS2 will use the effective genome
size value provided in `gsize` parameter.

If you are using a non-model organism or simply don’t want to use the
effective genome size provided in literature or MACS2 documentation, you
can set `calc_effective_gsize = true` to calculate an effective genome
size using the target genome fasta `fasta` filepath and read-length.

-   run_macs2
-   no_control_macs2
-   gsize
-   calc_effective_gsize
-   read_length

<!-- -->

    ##     //MACS2 params
    ##     run_macs2                   = true
    ##     no_control_macs2            = false   // if true, do not use IgG control in peak calling
    ## 
    ##     gsize                       = 1.87e9 //default effective genome size for mouse from MACS2 
    ##     calc_effective_gsize        = false  //if true, will override the value in gsize parameter
    ##     read_length                 = 150    //if calc_effective_gsize, provide illumina read-length in base pairs

## Advanced Options

In the `nextflow.config`, you can define additional command line
arguments to the scientific software under [process
scope](https://www.nextflow.io/docs/latest/process.html#). You may use
the advanced options to change computational resources requested for
different processes. The CPUs and memory parameters can updated to
request a larger amount of resources like CPUs or memory if files are
large. You may also edit the commandline parameters for processes in the
workflow using the `ext.arg`
[directive](https://www.nextflow.io/docs/latest/process.html#directives).
Please be aware the **default command line parameters for `Bowtie2`
processes are already provided for both target and spike-in alignment**,
but can be edited.

The most commonly modified and important process parameters are listed
toward the top of the process scope in the `nextflow.config` file.

You can edit the command line parameters for `SEACR` and `MACS2`
parameters that often need to be re-run multiple times when deciding on
the appropriate peak-set to use. For example, `MACS2` broad and narrow
peak calling parameters for different histone modifications which can be
modified using the `ext.args` parameter.

    ##  [1] // Computational resource allocation for the processes run in the workflow                                     
    ##  [2] process {                                                                                                      
    ##  [3]     //Bowtie2 aligner process specific parameters                                                              
    ##  [4]     withName: BOWTIE2_ALIGN {                                                                                  
    ##  [5]         cpus = { 2 * task.attempt }                                                                            
    ##  [6]         memory = { 32.GB * task.attempt }                                                                      
    ##  [7]         ext.prefix = { "${meta.id}.sort" }                                                                     
    ##  [8]         ext.args = '--local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700'
    ##  [9]         ext.args2 = ''      //command line arguments for `samtools sort`                                       
    ## [10]     }                                                                                                          
    ## [11]     //SEACR peak calling resources                                                                             
    ## [12]     withName: SEACR_CALLPEAK {                                                                                 
    ## [13]         cpus = { 1 * task.attempt }                                                                            
    ## [14]         memory = { 16.GB * task.attempt }                                                                      
    ## [15]         ext.version = '1.4' //version 1.3 and 1.4 supported                                                    
    ## [16]         ext.args = '--normalize norm --mode stringent --remove yes'                                            
    ## [17]         publishDir = [...]                                                                                     
    ## [18]                                                                                                                
    ## [19]     }                                                                                                          
    ## [20]     //MACS2 peak calling resources                                                                             
    ## [21]     withName: MACS2_CALLPEAK {                                                                                 
    ## [22]         cpus = { 1 * task.attempt }                                                                            
    ## [23]         memory = { 16.GB * task.attempt }                                                                      
    ## [24]         ext.args = '-q 0.01 --keep-dup all --bdg'                                                              
    ## [25]         publishDir = [...]                                                                                     
    ## [26]                                                                                                                
    ## [27]     }                                                                                                          
    ## [28]     //BAMCOVERAGE bigwig file  parameters                                                                      
    ## [29]     withName: DEEPTOOLS_BAMCOVERAGE {                                                                          
    ## [30]         cpus = { 4 * task.attempt }                                                                            
    ## [31]         memory = { 16.GB * task.attempt }                                                                      
    ## [32]         ext.args = '--normalizeUsing CPM --centerReads --verbose'                                              
    ## [33]     }

SEACR has the option to be set to SEACR v1.4 or SEACR v1.3 - which have
particularly different commandline interfaces, changes in the methods
for normalization to IgG, and v1.4 can optionally remove peaks found in
IgG. Please see [here](https://github.com/FredHutch/SEACR/releases) for
the full changelog.

For SEACR v1.3, Often, you will need to change `SEACR` from “non” to
“norm” for different normalization strategies whether you’re using IgG
normalization or spike-in normalization. The example below demonstrates
how to change the commandline params and version by editing
`ext.version` and `ext.args`.

    ##     //SEACR peak calling resources
    ##     withName: SEACR_CALLPEAK {
    ##         cpus = { 1 * task.attempt }
    ##         memory = { 16.GB * task.attempt }
    ##         ext.version = '1.3' //version 1.3 and 1.4 supported
    ##         ext.args = 'norm stringent'
    ##         publishDir = [...]
    ##                     
    ##     }

# Run Script

``` r
usethis::edit_file(here::here("main_run.sh"))
```

    ## • Edit '/active/taylor_s/people/jsmi26/RPDEV/cutandrun_nf/main_run.sh'

Decide on the `NFX_PROFILE`, which allows you to run the processes
either locally, or using the PBS job scheduler on Cybertron, and
determine if you’d like to use singularity containers or docker
containers.

2)  `PBS_singularity` \[DEFAULT, recommended\] \* you can submit a PBS
    job that will use singularity containers on Cybertron \* This takes
    care of requesting the appropriate resources using PBS

3)  `local_singularity` \* locally on an interactive session Cybertron
    with singularity \* requires appropriate computational resources be
    requested using
    `qsub -I -q <queue_name> -P <project_code> -l select=1:ncpus=4:mem=32GB`

Edit the script `main_run.sh` and change the values for the
`NFX_PROFILE` variable if desired.

    ## #Options: 'PBS_apptainer','local_apptainer','local_singularity', 'PBS_singularity'
    ## NFX_PROFILE='PBS_apptainer'

## Alignment and Peak Calls

Edit the variables in the `main_run.sh` script for entry-point of the
workflow. The default option *“align_call_peaks”* for the `NFX_ENTRY`
will run the full pipeline (QC, alignment, peak calling, coverage
tracks).

    ## #Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
    ## NFX_ENTRY='align_call_peaks'

If you already have aligned BAM files, see
`test_data/test_dataset_bams_sample_sheet.csv` for an example to call
peaks only using the entry `call_peaks`.

    ## #Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
    ## NFX_ENTRY='call_peaks'

Then, execute the `main_run.sh` script in order to complete the peak
calling on the samples. Provide a small descriptive prefix for the
pipeline run.

``` bash
./main_run.sh "my_analysis"
```

## Optional: Build the Index and Exit Pipeline

You can also change the entry-point of the workflow, which is
accomplished by setting the `NFX_ENTRY` variable in the `main_run.sh`
script to be `bowtie2_index_only`. This will allow the pipeline to run
only the Bowtie2 build process and exit upon completion of the index
building step.

    ## #Options: 'bowtie2_index_only', 'align_call_peaks', 'call_peaks'
    ## NFX_ENTRY='bowtie2_index_only'

``` bash
./main_run.sh "bowtie2_index"
```

# Expected Outputs

Under the path provided in the nextflow config for params “outdir”, you
will find directories named for each of the modules.

### Final Outputs

`results/{params.outdir}` - `samtools_view/` - aligned, coordinate
sorted, marked duplicates, and optionally quality filtered bam file

    -   {sample_id}.markedDup.filter.bam

-   `samtools_index/`
    -   {sample_id}.markedDup.filter.sort.bam.bai
-   `deeptools_bamcoverage/`
    -   Counts per million normalized coverage track (bigwig) file.
    -   {sample_id}\_CPM.bigWig
-   `{params.peaks_outdir}/seacr_callpeak/`
    -   If using IgG normalization, the {sample_id} of the IgG control
        used is appended to the target
    -   {sample_id}\_\[norm,non\].\[stringent,relaxed\].bed
-   `{params.peaks_outdir}/macs2_callpeak/`
    -   Optional output if `run_macs2 = true`
    -   {sample_id}\_peaks.\[narrowPeak,broadPeak\]

### File Structure

There will be the following file structure:

    ## ../results/mouse
    ## ├── bamtobedgraph
    ## │   ├── M1_H3K27_NK_aligned.bed
    ## │   ├── M1_H3K27_NK_aligned.clean.bed
    ## │   ├── M1_H3K27_NK_aligned.fragments.bed
    ## │   ├── M1_H3K27_NK_aligned_fragments.bg
    ## │   ├── M1_H3K4_NK_aligned.bed
    ## │   ├── M1_H3K4_NK_aligned.clean.bed
    ## │   ├── M1_H3K4_NK_aligned.fragments.bed
    ## │   ├── M1_H3K4_NK_aligned_fragments.bg
    ## │   ├── M1_IgG_NK_aligned.bed
    ## │   ├── M1_IgG_NK_aligned.clean.bed
    ## │   ├── M1_IgG_NK_aligned.fragments.bed
    ## │   ├── M1_IgG_NK_aligned_fragments.bg
    ## │   ├── M2_H3K27_NK_aligned.bed
    ## │   ├── M2_H3K27_NK_aligned.clean.bed
    ## │   ├── M2_H3K27_NK_aligned.fragments.bed
    ## │   ├── M2_H3K27_NK_aligned_fragments.bg
    ## │   ├── M2_H3K4_NK_aligned.bed
    ## │   ├── M2_H3K4_NK_aligned.clean.bed
    ## │   ├── M2_H3K4_NK_aligned.fragments.bed
    ## │   ├── M2_H3K4_NK_aligned_fragments.bg
    ## │   ├── M2_IgG_NK_aligned.bed
    ## │   ├── M2_IgG_NK_aligned.clean.bed
    ## │   ├── M2_IgG_NK_aligned.fragments.bed
    ## │   └── M2_IgG_NK_aligned_fragments.bg
    ## ├── bowtie2_align
    ## │   ├── M1_H3K27_NK.sort.bam
    ## │   ├── M1_H3K27_NK.sort.bowtie2.log
    ## │   ├── M1_H3K4_NK.sort.bam
    ## │   ├── M1_H3K4_NK.sort.bowtie2.log
    ## │   ├── M1_IgG_NK.sort.bam
    ## │   ├── M1_IgG_NK.sort.bowtie2.log
    ## │   ├── M2_H3K27_NK.sort.bam
    ## │   ├── M2_H3K27_NK.sort.bowtie2.log
    ## │   ├── M2_H3K4_NK.sort.bam
    ## │   ├── M2_H3K4_NK.sort.bowtie2.log
    ## │   ├── M2_IgG_NK.sort.bam
    ## │   └── M2_IgG_NK.sort.bowtie2.log
    ## ├── deeptools_bamcoverage
    ## │   ├── M1_H3K27_NK_CPM.bigWig
    ## │   ├── M1_H3K4_NK_CPM.bigWig
    ## │   ├── M1_IgG_NK_CPM.bigWig
    ## │   ├── M2_H3K27_NK_CPM.bigWig
    ## │   ├── M2_H3K4_NK_CPM.bigWig
    ## │   └── M2_IgG_NK_CPM.bigWig
    ## ├── deeptools_multibigwigsummary
    ## │   └── test_dataset_sample_sheet_scores_per_bin.npz
    ## ├── deeptools_plotcorrelation
    ## │   ├── test_dataset_sample_sheet.plotCorrelation.mat.tab
    ## │   └── test_dataset_sample_sheet.plotCorrelation.pdf
    ## ├── deeptools_plotfingerprint
    ## │   ├── M1_H3K27_NK.plotFingerprint.pdf
    ## │   ├── M1_H3K27_NK.plotFingerprint.qcmetrics.txt
    ## │   ├── M1_H3K27_NK.plotFingerprint.raw.txt
    ## │   ├── M1_H3K4_NK.plotFingerprint.pdf
    ## │   ├── M1_H3K4_NK.plotFingerprint.qcmetrics.txt
    ## │   ├── M1_H3K4_NK.plotFingerprint.raw.txt
    ## │   ├── M1_IgG_NK.plotFingerprint.pdf
    ## │   ├── M1_IgG_NK.plotFingerprint.qcmetrics.txt
    ## │   ├── M1_IgG_NK.plotFingerprint.raw.txt
    ## │   ├── M2_H3K27_NK.plotFingerprint.pdf
    ## │   ├── M2_H3K27_NK.plotFingerprint.qcmetrics.txt
    ## │   ├── M2_H3K27_NK.plotFingerprint.raw.txt
    ## │   ├── M2_H3K4_NK.plotFingerprint.pdf
    ## │   ├── M2_H3K4_NK.plotFingerprint.qcmetrics.txt
    ## │   ├── M2_H3K4_NK.plotFingerprint.raw.txt
    ## │   ├── M2_IgG_NK.plotFingerprint.pdf
    ## │   ├── M2_IgG_NK.plotFingerprint.qcmetrics.txt
    ## │   ├── M2_IgG_NK.plotFingerprint.raw.txt
    ## │   ├── test_dataset_sample_sheet.plotFingerprint.pdf
    ## │   ├── test_dataset_sample_sheet.plotFingerprint.qcmetrics.txt
    ## │   └── test_dataset_sample_sheet.plotFingerprint.raw.txt
    ## ├── deeptools_plotpca
    ## │   ├── test_dataset_sample_sheet.plotPCA.pdf
    ## │   └── test_dataset_sample_sheet.plotPCA.tab
    ## ├── fastqc
    ## │   ├── M1_H3K27_NK_FASTQC
    ## │   ├── M1_H3K4_NK_FASTQC
    ## │   ├── M1_IgG_NK_FASTQC
    ## │   ├── M2_H3K27_NK_FASTQC
    ## │   ├── M2_H3K4_NK_FASTQC
    ## │   └── M2_IgG_NK_FASTQC
    ## ├── fastqc_trim
    ## │   ├── M1_H3K27_NK_FASTQC_TRIM
    ## │   ├── M1_H3K4_NK_FASTQC_TRIM
    ## │   ├── M1_IgG_NK_FASTQC_TRIM
    ## │   ├── M2_H3K27_NK_FASTQC_TRIM
    ## │   ├── M2_H3K4_NK_FASTQC_TRIM
    ## │   └── M2_IgG_NK_FASTQC_TRIM
    ## ├── macs_igg
    ## │   ├── macs2_callpeak
    ## │   ├── macs2_plotenrichment
    ## │   ├── macspeakstobed
    ## │   ├── seacr_callpeak
    ## │   └── seacr_plotenrichment
    ## ├── macs_no_igg
    ## │   ├── macs2_callpeak
    ## │   ├── macs2_plotenrichment
    ## │   ├── macspeakstobed
    ## │   ├── seacr_callpeak
    ## │   └── seacr_plotenrichment
    ## ├── multiqc
    ## │   ├── test_dataset_sample_sheet_multiqc_report.html
    ## │   └── test_dataset_sample_sheet_multiqc_report_data
    ## ├── picard_markduplicates
    ## │   ├── M1_H3K27_NK.markedDup.MarkDuplicates.metrics.txt
    ## │   ├── M1_H3K27_NK.markedDup.bai
    ## │   ├── M1_H3K27_NK.markedDup.bam
    ## │   ├── M1_H3K27_NK.markedDup.bam.md5
    ## │   ├── M1_H3K4_NK.markedDup.MarkDuplicates.metrics.txt
    ## │   ├── M1_H3K4_NK.markedDup.bai
    ## │   ├── M1_H3K4_NK.markedDup.bam
    ## │   ├── M1_H3K4_NK.markedDup.bam.md5
    ## │   ├── M1_IgG_NK.markedDup.MarkDuplicates.metrics.txt
    ## │   ├── M1_IgG_NK.markedDup.bai
    ## │   ├── M1_IgG_NK.markedDup.bam
    ## │   ├── M1_IgG_NK.markedDup.bam.md5
    ## │   ├── M2_H3K27_NK.markedDup.MarkDuplicates.metrics.txt
    ## │   ├── M2_H3K27_NK.markedDup.bai
    ## │   ├── M2_H3K27_NK.markedDup.bam
    ## │   ├── M2_H3K27_NK.markedDup.bam.md5
    ## │   ├── M2_H3K4_NK.markedDup.MarkDuplicates.metrics.txt
    ## │   ├── M2_H3K4_NK.markedDup.bai
    ## │   ├── M2_H3K4_NK.markedDup.bam
    ## │   ├── M2_H3K4_NK.markedDup.bam.md5
    ## │   ├── M2_IgG_NK.markedDup.MarkDuplicates.metrics.txt
    ## │   ├── M2_IgG_NK.markedDup.bai
    ## │   ├── M2_IgG_NK.markedDup.bam
    ## │   └── M2_IgG_NK.markedDup.bam.md5
    ## ├── samtools_faidx
    ## │   └── mm39.fa.fai
    ## ├── samtools_index
    ## │   ├── M1_H3K27_NK.markedDup.filter.sort.bam.bai
    ## │   ├── M1_H3K4_NK.markedDup.filter.sort.bam.bai
    ## │   ├── M1_IgG_NK.markedDup.filter.sort.bam.bai
    ## │   ├── M2_H3K27_NK.markedDup.filter.sort.bam.bai
    ## │   ├── M2_H3K4_NK.markedDup.filter.sort.bam.bai
    ## │   └── M2_IgG_NK.markedDup.filter.sort.bam.bai
    ## ├── samtools_nsort
    ## │   ├── M1_H3K27_NK.markedDup.filter.nsort.bam
    ## │   ├── M1_H3K4_NK.markedDup.filter.nsort.bam
    ## │   ├── M1_IgG_NK.markedDup.filter.nsort.bam
    ## │   ├── M2_H3K27_NK.markedDup.filter.nsort.bam
    ## │   ├── M2_H3K4_NK.markedDup.filter.nsort.bam
    ## │   └── M2_IgG_NK.markedDup.filter.nsort.bam
    ## ├── samtools_sort
    ## │   ├── M1_H3K27_NK.markedDup.filter.sort.bam
    ## │   ├── M1_H3K4_NK.markedDup.filter.sort.bam
    ## │   ├── M1_IgG_NK.markedDup.filter.sort.bam
    ## │   ├── M2_H3K27_NK.markedDup.filter.sort.bam
    ## │   ├── M2_H3K4_NK.markedDup.filter.sort.bam
    ## │   └── M2_IgG_NK.markedDup.filter.sort.bam
    ## ├── samtools_stats
    ## │   ├── M1_H3K27_NK.markedDup.stats
    ## │   ├── M1_H3K4_NK.markedDup.stats
    ## │   ├── M1_IgG_NK.markedDup.stats
    ## │   ├── M2_H3K27_NK.markedDup.stats
    ## │   ├── M2_H3K4_NK.markedDup.stats
    ## │   └── M2_IgG_NK.markedDup.stats
    ## ├── samtools_view
    ## │   ├── M1_H3K27_NK.markedDup.filter.bam
    ## │   ├── M1_H3K4_NK.markedDup.filter.bam
    ## │   ├── M1_IgG_NK.markedDup.filter.bam
    ## │   ├── M2_H3K27_NK.markedDup.filter.bam
    ## │   ├── M2_H3K4_NK.markedDup.filter.bam
    ## │   └── M2_IgG_NK.markedDup.filter.bam
    ## └── trimgalore
    ##     ├── M1_H3K27_NK_1.fastq.gz_trimming_report.txt
    ##     ├── M1_H3K27_NK_1_val_1.fq.gz
    ##     ├── M1_H3K27_NK_2.fastq.gz_trimming_report.txt
    ##     ├── M1_H3K27_NK_2_val_2.fq.gz
    ##     ├── M1_H3K4_NK_1.fastq.gz_trimming_report.txt
    ##     ├── M1_H3K4_NK_1_val_1.fq.gz
    ##     ├── M1_H3K4_NK_2.fastq.gz_trimming_report.txt
    ##     ├── M1_H3K4_NK_2_val_2.fq.gz
    ##     ├── M1_IgG_NK_1.fastq.gz_trimming_report.txt
    ##     ├── M1_IgG_NK_1_val_1.fq.gz
    ##     ├── M1_IgG_NK_2.fastq.gz_trimming_report.txt
    ##     ├── M1_IgG_NK_2_val_2.fq.gz
    ##     ├── M2_H3K27_NK_1.fastq.gz_trimming_report.txt
    ##     ├── M2_H3K27_NK_1_val_1.fq.gz
    ##     ├── M2_H3K27_NK_2.fastq.gz_trimming_report.txt
    ##     ├── M2_H3K27_NK_2_val_2.fq.gz
    ##     ├── M2_H3K4_NK_1.fastq.gz_trimming_report.txt
    ##     ├── M2_H3K4_NK_1_val_1.fq.gz
    ##     ├── M2_H3K4_NK_2.fastq.gz_trimming_report.txt
    ##     ├── M2_H3K4_NK_2_val_2.fq.gz
    ##     ├── M2_IgG_NK_1.fastq.gz_trimming_report.txt
    ##     ├── M2_IgG_NK_1_val_1.fq.gz
    ##     ├── M2_IgG_NK_2.fastq.gz_trimming_report.txt
    ##     └── M2_IgG_NK_2_val_2.fq.gz

### Detailed File Structure

Within each directory you will find the following files (top 5 files per
directory are shown):

| path                                                                                         | type      | process                           | filename                                          |
|:---------------------------------------------------------------------------------------------|:----------|:----------------------------------|:--------------------------------------------------|
| ../results/mouse/bamtobedgraph                                                               | directory | /bamtobedgraph                    |                                                   |
| ../results/mouse/bamtobedgraph/M1_H3K27_NK_aligned.bed                                       | file      | /bamtobedgraph                    | M1_H3K27_NK_aligned.bed                           |
| ../results/mouse/bamtobedgraph/M1_H3K27_NK_aligned.clean.bed                                 | file      | /bamtobedgraph                    | M1_H3K27_NK_aligned.clean.bed                     |
| ../results/mouse/bamtobedgraph/M1_H3K27_NK_aligned.fragments.bed                             | file      | /bamtobedgraph                    | M1_H3K27_NK_aligned.fragments.bed                 |
| ../results/mouse/bamtobedgraph/M1_H3K27_NK_aligned_fragments.bg                              | file      | /bamtobedgraph                    | M1_H3K27_NK_aligned_fragments.bg                  |
| ../results/mouse/bowtie2_align                                                               | directory | /bowtie2_align                    |                                                   |
| ../results/mouse/bowtie2_align/M1_H3K27_NK.sort.bam                                          | file      | /bowtie2_align                    | M1_H3K27_NK.sort.bam                              |
| ../results/mouse/bowtie2_align/M1_H3K27_NK.sort.bowtie2.log                                  | file      | /bowtie2_align                    | M1_H3K27_NK.sort.bowtie2.log                      |
| ../results/mouse/bowtie2_align/M1_H3K4_NK.sort.bam                                           | file      | /bowtie2_align                    | M1_H3K4_NK.sort.bam                               |
| ../results/mouse/bowtie2_align/M1_H3K4_NK.sort.bowtie2.log                                   | file      | /bowtie2_align                    | M1_H3K4_NK.sort.bowtie2.log                       |
| ../results/mouse/deeptools_bamcoverage                                                       | directory | /deeptools_bamcoverage            |                                                   |
| ../results/mouse/deeptools_bamcoverage/M1_H3K27_NK_CPM.bigWig                                | file      | /deeptools_bamcoverage            | M1_H3K27_NK_CPM.bigWig                            |
| ../results/mouse/deeptools_bamcoverage/M1_H3K4_NK_CPM.bigWig                                 | file      | /deeptools_bamcoverage            | M1_H3K4_NK_CPM.bigWig                             |
| ../results/mouse/deeptools_bamcoverage/M1_IgG_NK_CPM.bigWig                                  | file      | /deeptools_bamcoverage            | M1_IgG_NK_CPM.bigWig                              |
| ../results/mouse/deeptools_bamcoverage/M2_H3K27_NK_CPM.bigWig                                | file      | /deeptools_bamcoverage            | M2_H3K27_NK_CPM.bigWig                            |
| ../results/mouse/deeptools_multibigwigsummary                                                | directory | /deeptools_multibigwigsummary     |                                                   |
| ../results/mouse/deeptools_multibigwigsummary/test_dataset_sample_sheet_scores_per_bin.npz   | file      | /deeptools_multibigwigsummary     | test_dataset_sample_sheet_scores_per_bin.npz      |
| ../results/mouse/deeptools_plotcorrelation                                                   | directory | /deeptools_plotcorrelation        |                                                   |
| ../results/mouse/deeptools_plotcorrelation/test_dataset_sample_sheet.plotCorrelation.mat.tab | file      | /deeptools_plotcorrelation        | test_dataset_sample_sheet.plotCorrelation.mat.tab |
| ../results/mouse/deeptools_plotcorrelation/test_dataset_sample_sheet.plotCorrelation.pdf     | file      | /deeptools_plotcorrelation        | test_dataset_sample_sheet.plotCorrelation.pdf     |
| ../results/mouse/deeptools_plotfingerprint                                                   | directory | /deeptools_plotfingerprint        |                                                   |
| ../results/mouse/deeptools_plotfingerprint/M1_H3K27_NK.plotFingerprint.pdf                   | file      | /deeptools_plotfingerprint        | M1_H3K27_NK.plotFingerprint.pdf                   |
| ../results/mouse/deeptools_plotfingerprint/M1_H3K27_NK.plotFingerprint.qcmetrics.txt         | file      | /deeptools_plotfingerprint        | M1_H3K27_NK.plotFingerprint.qcmetrics.txt         |
| ../results/mouse/deeptools_plotfingerprint/M1_H3K27_NK.plotFingerprint.raw.txt               | file      | /deeptools_plotfingerprint        | M1_H3K27_NK.plotFingerprint.raw.txt               |
| ../results/mouse/deeptools_plotfingerprint/M1_H3K4_NK.plotFingerprint.pdf                    | file      | /deeptools_plotfingerprint        | M1_H3K4_NK.plotFingerprint.pdf                    |
| ../results/mouse/deeptools_plotpca                                                           | directory | /deeptools_plotpca                |                                                   |
| ../results/mouse/deeptools_plotpca/test_dataset_sample_sheet.plotPCA.pdf                     | file      | /deeptools_plotpca                | test_dataset_sample_sheet.plotPCA.pdf             |
| ../results/mouse/deeptools_plotpca/test_dataset_sample_sheet.plotPCA.tab                     | file      | /deeptools_plotpca                | test_dataset_sample_sheet.plotPCA.tab             |
| ../results/mouse/fastqc                                                                      | directory | /fastqc                           |                                                   |
| ../results/mouse/fastqc/M1_H3K27_NK_FASTQC                                                   | directory | /fastqc                           |                                                   |
| ../results/mouse/fastqc/M1_H3K27_NK_FASTQC/M1_H3K27_NK_1\_fastqc.html                        | file      | /fastqc                           | M1_H3K27_NK_1\_fastqc.html                        |
| ../results/mouse/fastqc/M1_H3K27_NK_FASTQC/M1_H3K27_NK_1\_fastqc.zip                         | file      | /fastqc                           | M1_H3K27_NK_1\_fastqc.zip                         |
| ../results/mouse/fastqc/M1_H3K27_NK_FASTQC/M1_H3K27_NK_2\_fastqc.html                        | file      | /fastqc                           | M1_H3K27_NK_2\_fastqc.html                        |
| ../results/mouse/fastqc_trim                                                                 | directory | /fastqc_trim                      |                                                   |
| ../results/mouse/fastqc_trim/M1_H3K27_NK_FASTQC_TRIM                                         | directory | /fastqc_trim                      |                                                   |
| ../results/mouse/fastqc_trim/M1_H3K27_NK_FASTQC_TRIM/M1_H3K27_NK_1\_fastqc.html              | file      | /fastqc_trim                      | M1_H3K27_NK_1\_fastqc.html                        |
| ../results/mouse/fastqc_trim/M1_H3K27_NK_FASTQC_TRIM/M1_H3K27_NK_1\_fastqc.zip               | file      | /fastqc_trim                      | M1_H3K27_NK_1\_fastqc.zip                         |
| ../results/mouse/fastqc_trim/M1_H3K27_NK_FASTQC_TRIM/M1_H3K27_NK_2\_fastqc.html              | file      | /fastqc_trim                      | M1_H3K27_NK_2\_fastqc.html                        |
| ../results/mouse/macs_igg                                                                    | directory | /macs_igg                         |                                                   |
| ../results/mouse/macs_igg/macs2_callpeak                                                     | directory | /macs_igg/macs2_callpeak          |                                                   |
| ../results/mouse/macs_igg/macs2_callpeak/M1_H3K27_NK_control_lambda.bdg                      | file      | /macs_igg/macs2_callpeak          | M1_H3K27_NK_control_lambda.bdg                    |
| ../results/mouse/macs_igg/macs2_callpeak/M1_H3K27_NK_peaks.narrowPeak                        | file      | /macs_igg/macs2_callpeak          | M1_H3K27_NK_peaks.narrowPeak                      |
| ../results/mouse/macs_igg/macs2_callpeak/M1_H3K27_NK_peaks.xls                               | file      | /macs_igg/macs2_callpeak          | M1_H3K27_NK_peaks.xls                             |
| ../results/mouse/macs_igg/macs2_callpeak/M1_H3K27_NK_summits.bed                             | file      | /macs_igg/macs2_callpeak          | M1_H3K27_NK_summits.bed                           |
| ../results/mouse/macs_igg/macs2_plotenrichment                                               | directory | /macs_igg/macs2_plotenrichment    |                                                   |
| ../results/mouse/macs_igg/macs2_plotenrichment/M1_H3K27_NK.plotEnrichment.pdf                | file      | /macs_igg/macs2_plotenrichment    | M1_H3K27_NK.plotEnrichment.pdf                    |
| ../results/mouse/macs_igg/macs2_plotenrichment/M1_H3K27_NK.plotEnrichment.txt                | file      | /macs_igg/macs2_plotenrichment    | M1_H3K27_NK.plotEnrichment.txt                    |
| ../results/mouse/macs_igg/macs2_plotenrichment/M1_H3K4_NK.plotEnrichment.pdf                 | file      | /macs_igg/macs2_plotenrichment    | M1_H3K4_NK.plotEnrichment.pdf                     |
| ../results/mouse/macs_igg/macs2_plotenrichment/M1_H3K4_NK.plotEnrichment.txt                 | file      | /macs_igg/macs2_plotenrichment    | M1_H3K4_NK.plotEnrichment.txt                     |
| ../results/mouse/macs_igg/macspeakstobed                                                     | directory | /macs_igg/macspeakstobed          |                                                   |
| ../results/mouse/macs_igg/macspeakstobed/M1_H3K27_NK_peaks.bed                               | file      | /macs_igg/macspeakstobed          | M1_H3K27_NK_peaks.bed                             |
| ../results/mouse/macs_igg/macspeakstobed/M1_H3K4_NK_peaks.bed                                | file      | /macs_igg/macspeakstobed          | M1_H3K4_NK_peaks.bed                              |
| ../results/mouse/macs_igg/macspeakstobed/M2_H3K27_NK_peaks.bed                               | file      | /macs_igg/macspeakstobed          | M2_H3K27_NK_peaks.bed                             |
| ../results/mouse/macs_igg/macspeakstobed/M2_H3K4_NK_peaks.bed                                | file      | /macs_igg/macspeakstobed          | M2_H3K4_NK_peaks.bed                              |
| ../results/mouse/macs_igg/seacr_callpeak                                                     | directory | /macs_igg/seacr_callpeak          |                                                   |
| ../results/mouse/macs_igg/seacr_callpeak/M1_H3K27_NK_vs_M1_IgG_NK_norm.stringent.bed         | file      | /macs_igg/seacr_callpeak          | M1_H3K27_NK_vs_M1_IgG_NK_norm.stringent.bed       |
| ../results/mouse/macs_igg/seacr_callpeak/M1_H3K4_NK_vs_M1_IgG_NK_norm.stringent.bed          | file      | /macs_igg/seacr_callpeak          | M1_H3K4_NK_vs_M1_IgG_NK_norm.stringent.bed        |
| ../results/mouse/macs_igg/seacr_callpeak/M2_H3K27_NK_vs_M2_IgG_NK_norm.stringent.bed         | file      | /macs_igg/seacr_callpeak          | M2_H3K27_NK_vs_M2_IgG_NK_norm.stringent.bed       |
| ../results/mouse/macs_igg/seacr_callpeak/M2_H3K4_NK_vs_M2_IgG_NK_norm.stringent.bed          | file      | /macs_igg/seacr_callpeak          | M2_H3K4_NK_vs_M2_IgG_NK_norm.stringent.bed        |
| ../results/mouse/macs_igg/seacr_plotenrichment                                               | directory | /macs_igg/seacr_plotenrichment    |                                                   |
| ../results/mouse/macs_igg/seacr_plotenrichment/M1_H3K27_NK.plotEnrichment.pdf                | file      | /macs_igg/seacr_plotenrichment    | M1_H3K27_NK.plotEnrichment.pdf                    |
| ../results/mouse/macs_igg/seacr_plotenrichment/M1_H3K27_NK.plotEnrichment.txt                | file      | /macs_igg/seacr_plotenrichment    | M1_H3K27_NK.plotEnrichment.txt                    |
| ../results/mouse/macs_igg/seacr_plotenrichment/M1_H3K4_NK.plotEnrichment.pdf                 | file      | /macs_igg/seacr_plotenrichment    | M1_H3K4_NK.plotEnrichment.pdf                     |
| ../results/mouse/macs_igg/seacr_plotenrichment/M1_H3K4_NK.plotEnrichment.txt                 | file      | /macs_igg/seacr_plotenrichment    | M1_H3K4_NK.plotEnrichment.txt                     |
| ../results/mouse/macs_no_igg                                                                 | directory | /macs_no_igg                      |                                                   |
| ../results/mouse/macs_no_igg/macs2_callpeak                                                  | directory | /macs_no_igg/macs2_callpeak       |                                                   |
| ../results/mouse/macs_no_igg/macs2_callpeak/M1_H3K27_NK_control_lambda.bdg                   | file      | /macs_no_igg/macs2_callpeak       | M1_H3K27_NK_control_lambda.bdg                    |
| ../results/mouse/macs_no_igg/macs2_callpeak/M1_H3K27_NK_peaks.narrowPeak                     | file      | /macs_no_igg/macs2_callpeak       | M1_H3K27_NK_peaks.narrowPeak                      |
| ../results/mouse/macs_no_igg/macs2_callpeak/M1_H3K27_NK_peaks.xls                            | file      | /macs_no_igg/macs2_callpeak       | M1_H3K27_NK_peaks.xls                             |
| ../results/mouse/macs_no_igg/macs2_callpeak/M1_H3K27_NK_summits.bed                          | file      | /macs_no_igg/macs2_callpeak       | M1_H3K27_NK_summits.bed                           |
| ../results/mouse/macs_no_igg/macs2_plotenrichment                                            | directory | /macs_no_igg/macs2_plotenrichment |                                                   |
| ../results/mouse/macs_no_igg/macs2_plotenrichment/M1_H3K27_NK.plotEnrichment.pdf             | file      | /macs_no_igg/macs2_plotenrichment | M1_H3K27_NK.plotEnrichment.pdf                    |
| ../results/mouse/macs_no_igg/macs2_plotenrichment/M1_H3K27_NK.plotEnrichment.txt             | file      | /macs_no_igg/macs2_plotenrichment | M1_H3K27_NK.plotEnrichment.txt                    |
| ../results/mouse/macs_no_igg/macs2_plotenrichment/M1_H3K4_NK.plotEnrichment.pdf              | file      | /macs_no_igg/macs2_plotenrichment | M1_H3K4_NK.plotEnrichment.pdf                     |
| ../results/mouse/macs_no_igg/macs2_plotenrichment/M1_H3K4_NK.plotEnrichment.txt              | file      | /macs_no_igg/macs2_plotenrichment | M1_H3K4_NK.plotEnrichment.txt                     |
| ../results/mouse/macs_no_igg/macspeakstobed                                                  | directory | /macs_no_igg/macspeakstobed       |                                                   |
| ../results/mouse/macs_no_igg/macspeakstobed/M1_H3K27_NK_peaks.bed                            | file      | /macs_no_igg/macspeakstobed       | M1_H3K27_NK_peaks.bed                             |
| ../results/mouse/macs_no_igg/macspeakstobed/M1_H3K4_NK_peaks.bed                             | file      | /macs_no_igg/macspeakstobed       | M1_H3K4_NK_peaks.bed                              |
| ../results/mouse/macs_no_igg/macspeakstobed/M2_H3K27_NK_peaks.bed                            | file      | /macs_no_igg/macspeakstobed       | M2_H3K27_NK_peaks.bed                             |
| ../results/mouse/macs_no_igg/macspeakstobed/M2_H3K4_NK_peaks.bed                             | file      | /macs_no_igg/macspeakstobed       | M2_H3K4_NK_peaks.bed                              |
| ../results/mouse/macs_no_igg/seacr_callpeak                                                  | directory | /macs_no_igg/seacr_callpeak       |                                                   |
| ../results/mouse/macs_no_igg/seacr_callpeak/M1_H3K27_NK_vs_M1_IgG_NK_norm.stringent.bed      | file      | /macs_no_igg/seacr_callpeak       | M1_H3K27_NK_vs_M1_IgG_NK_norm.stringent.bed       |
| ../results/mouse/macs_no_igg/seacr_callpeak/M1_H3K4_NK_vs_M1_IgG_NK_norm.stringent.bed       | file      | /macs_no_igg/seacr_callpeak       | M1_H3K4_NK_vs_M1_IgG_NK_norm.stringent.bed        |
| ../results/mouse/macs_no_igg/seacr_callpeak/M2_H3K27_NK_vs_M2_IgG_NK_norm.stringent.bed      | file      | /macs_no_igg/seacr_callpeak       | M2_H3K27_NK_vs_M2_IgG_NK_norm.stringent.bed       |
| ../results/mouse/macs_no_igg/seacr_callpeak/M2_H3K4_NK_vs_M2_IgG_NK_norm.stringent.bed       | file      | /macs_no_igg/seacr_callpeak       | M2_H3K4_NK_vs_M2_IgG_NK_norm.stringent.bed        |
| ../results/mouse/macs_no_igg/seacr_plotenrichment                                            | directory | /macs_no_igg/seacr_plotenrichment |                                                   |
| ../results/mouse/macs_no_igg/seacr_plotenrichment/M1_H3K27_NK.plotEnrichment.pdf             | file      | /macs_no_igg/seacr_plotenrichment | M1_H3K27_NK.plotEnrichment.pdf                    |
| ../results/mouse/macs_no_igg/seacr_plotenrichment/M1_H3K27_NK.plotEnrichment.txt             | file      | /macs_no_igg/seacr_plotenrichment | M1_H3K27_NK.plotEnrichment.txt                    |
| ../results/mouse/macs_no_igg/seacr_plotenrichment/M1_H3K4_NK.plotEnrichment.pdf              | file      | /macs_no_igg/seacr_plotenrichment | M1_H3K4_NK.plotEnrichment.pdf                     |
| ../results/mouse/macs_no_igg/seacr_plotenrichment/M1_H3K4_NK.plotEnrichment.txt              | file      | /macs_no_igg/seacr_plotenrichment | M1_H3K4_NK.plotEnrichment.txt                     |
| ../results/mouse/multiqc                                                                     | directory | /multiqc                          |                                                   |
| ../results/mouse/picard_markduplicates                                                       | directory | /picard_markduplicates            |                                                   |
| ../results/mouse/picard_markduplicates/M1_H3K27_NK.markedDup.MarkDuplicates.metrics.txt      | file      | /picard_markduplicates            | M1_H3K27_NK.markedDup.MarkDuplicates.metrics.txt  |
| ../results/mouse/picard_markduplicates/M1_H3K27_NK.markedDup.bai                             | file      | /picard_markduplicates            | M1_H3K27_NK.markedDup.bai                         |
| ../results/mouse/picard_markduplicates/M1_H3K27_NK.markedDup.bam                             | file      | /picard_markduplicates            | M1_H3K27_NK.markedDup.bam                         |
| ../results/mouse/picard_markduplicates/M1_H3K27_NK.markedDup.bam.md5                         | file      | /picard_markduplicates            | M1_H3K27_NK.markedDup.bam.md5                     |
| ../results/mouse/samtools_faidx                                                              | directory | /samtools_faidx                   |                                                   |
| ../results/mouse/samtools_faidx/mm39.fa.fai                                                  | file      | /samtools_faidx                   | mm39.fa.fai                                       |
| ../results/mouse/samtools_index                                                              | directory | /samtools_index                   |                                                   |
| ../results/mouse/samtools_index/M1_H3K27_NK.markedDup.filter.sort.bam.bai                    | file      | /samtools_index                   | M1_H3K27_NK.markedDup.filter.sort.bam.bai         |
| ../results/mouse/samtools_index/M1_H3K4_NK.markedDup.filter.sort.bam.bai                     | file      | /samtools_index                   | M1_H3K4_NK.markedDup.filter.sort.bam.bai          |
| ../results/mouse/samtools_index/M1_IgG_NK.markedDup.filter.sort.bam.bai                      | file      | /samtools_index                   | M1_IgG_NK.markedDup.filter.sort.bam.bai           |
| ../results/mouse/samtools_index/M2_H3K27_NK.markedDup.filter.sort.bam.bai                    | file      | /samtools_index                   | M2_H3K27_NK.markedDup.filter.sort.bam.bai         |
| ../results/mouse/samtools_nsort                                                              | directory | /samtools_nsort                   |                                                   |
| ../results/mouse/samtools_nsort/M1_H3K27_NK.markedDup.filter.nsort.bam                       | file      | /samtools_nsort                   | M1_H3K27_NK.markedDup.filter.nsort.bam            |
| ../results/mouse/samtools_nsort/M1_H3K4_NK.markedDup.filter.nsort.bam                        | file      | /samtools_nsort                   | M1_H3K4_NK.markedDup.filter.nsort.bam             |
| ../results/mouse/samtools_nsort/M1_IgG_NK.markedDup.filter.nsort.bam                         | file      | /samtools_nsort                   | M1_IgG_NK.markedDup.filter.nsort.bam              |
| ../results/mouse/samtools_nsort/M2_H3K27_NK.markedDup.filter.nsort.bam                       | file      | /samtools_nsort                   | M2_H3K27_NK.markedDup.filter.nsort.bam            |
| ../results/mouse/samtools_sort                                                               | directory | /samtools_sort                    |                                                   |
| ../results/mouse/samtools_sort/M1_H3K27_NK.markedDup.filter.sort.bam                         | file      | /samtools_sort                    | M1_H3K27_NK.markedDup.filter.sort.bam             |
| ../results/mouse/samtools_sort/M1_H3K4_NK.markedDup.filter.sort.bam                          | file      | /samtools_sort                    | M1_H3K4_NK.markedDup.filter.sort.bam              |
| ../results/mouse/samtools_sort/M1_IgG_NK.markedDup.filter.sort.bam                           | file      | /samtools_sort                    | M1_IgG_NK.markedDup.filter.sort.bam               |
| ../results/mouse/samtools_sort/M2_H3K27_NK.markedDup.filter.sort.bam                         | file      | /samtools_sort                    | M2_H3K27_NK.markedDup.filter.sort.bam             |
| ../results/mouse/samtools_stats                                                              | directory | /samtools_stats                   |                                                   |
| ../results/mouse/samtools_stats/M1_H3K27_NK.markedDup.stats                                  | file      | /samtools_stats                   | M1_H3K27_NK.markedDup.stats                       |
| ../results/mouse/samtools_stats/M1_H3K4_NK.markedDup.stats                                   | file      | /samtools_stats                   | M1_H3K4_NK.markedDup.stats                        |
| ../results/mouse/samtools_stats/M1_IgG_NK.markedDup.stats                                    | file      | /samtools_stats                   | M1_IgG_NK.markedDup.stats                         |
| ../results/mouse/samtools_stats/M2_H3K27_NK.markedDup.stats                                  | file      | /samtools_stats                   | M2_H3K27_NK.markedDup.stats                       |
| ../results/mouse/samtools_view                                                               | directory | /samtools_view                    |                                                   |
| ../results/mouse/samtools_view/M1_H3K27_NK.markedDup.filter.bam                              | file      | /samtools_view                    | M1_H3K27_NK.markedDup.filter.bam                  |
| ../results/mouse/samtools_view/M1_H3K4_NK.markedDup.filter.bam                               | file      | /samtools_view                    | M1_H3K4_NK.markedDup.filter.bam                   |
| ../results/mouse/samtools_view/M1_IgG_NK.markedDup.filter.bam                                | file      | /samtools_view                    | M1_IgG_NK.markedDup.filter.bam                    |
| ../results/mouse/samtools_view/M2_H3K27_NK.markedDup.filter.bam                              | file      | /samtools_view                    | M2_H3K27_NK.markedDup.filter.bam                  |
| ../results/mouse/trimgalore                                                                  | directory | /trimgalore                       |                                                   |
| ../results/mouse/trimgalore/M1_H3K27_NK_1.fastq.gz_trimming_report.txt                       | file      | /trimgalore                       | M1_H3K27_NK_1.fastq.gz_trimming_report.txt        |
| ../results/mouse/trimgalore/M1_H3K27_NK_1\_val_1.fq.gz                                       | file      | /trimgalore                       | M1_H3K27_NK_1\_val_1.fq.gz                        |
| ../results/mouse/trimgalore/M1_H3K27_NK_2.fastq.gz_trimming_report.txt                       | file      | /trimgalore                       | M1_H3K27_NK_2.fastq.gz_trimming_report.txt        |
| ../results/mouse/trimgalore/M1_H3K27_NK_2\_val_2.fq.gz                                       | file      | /trimgalore                       | M1_H3K27_NK_2\_val_2.fq.gz                        |

## Pipeline Reports

In addition, there will be an HTML report with information on where the
temp data is stored in the `workDir` path, and general run statistics
such as resource utilized versus requested, which helps with
optimization. It will also provide information on how much walltime was
used per sample, total CPU hours, etc.

The HTML file is found in `reports` directory and will have the prefix
defined on the command line when the `./main_run.sh "my_analysis"` was
invoked, so in this example it would be named
“my_analysis\_{DATE}.html”.

There will also be a detailed nextflow log file that is useful for
de-bugging which will also be named in this example,
“my_analysis\_{DATE}\_nextflow.log”.

Finally, the pipeline will produce a DAG - Directed acyclic graph, which
describes the workflow channels (inputs) and the modules. The DAG image
will be saved under `dag/` directory with the name
“my_analysis\_{DATE}\_dag.pdf”.

<img src="images/dag.png" width="4120" style="display: block; margin: auto;" />
