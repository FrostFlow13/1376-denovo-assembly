# 1376-denovo-assembly

This repository contains the full code/processes I used for generation of de novo assemblies using deep long- and short-read sequencing. I cannot guarantee it will work for every organism or dataset, especially for organisms/strains that have low heterozygosity of their genomes or long (~75-100+ kb) tracts of homozygosity between heterozygous regions, but it worked for a CRISPR-competent strain of _Candida albicans_ SC5314 (specifically a Chr5AB disomic derivative of AHY940 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5422035/).

The code will be in a .sh file, but is not intended to be run as the file itself. I highly suggest downloading it and opening it in a software such as Notepad++ for nicer visualization.

Additionally, as a general rule of thumb, if a section of code DOESN'T start with "source /[HOMEPATH]/miniconda3/bin/activate", that means it's very likely something that needs to be done via the terminal! This is because the jobs (which is how I ran these scripts) need conda activated first.
