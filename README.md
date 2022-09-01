

# To Do:
1) Glenn meeting - overview the Cut&Run technology
2) create a test dataset of mouse samples that runs in < 5 min
   1) SEACR BED output --> chr17 subset bed --> subset the bam --> grab read IDs --> filter fastq for read IDs
3) Add samtools index module
4) test khmer module for effective genome size calculation
   *  https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
   *  https://khmer.readthedocs.io/en/stable/user/scripts.html#unique-kmers-py
5) Glenn will add fastq and multiqc modules from NF-core
6) Determine logic on whether to have broad or narrow peak calls
   1) should it be in the sample manifest or this that too much?
7) Add Module and test for picard Markduplicates and samtools stats
8) Add module and test module for the spike-in index building and alignement
9)  Add module and test spike-in normalization factor calculation 
10) Add module and test for the spike-in normalization script https://github.com/Henikoff/Cut-and-Run
