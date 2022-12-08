# Cut&Run Alignement, QC, and Peak Calling Pipeline

This pipeline uses publically available modules from [nf-core](https://nf-co.re/) with some locally created modules.

The primary functionality is to run a workflow on 10s - 1000s of samples in parallel on the Cybertron HPC using the PBS job scheduler and containerized scientific software.

The step-by-step instructions to run the workflow can be found in workflow_docs/workflow_run.Rmd