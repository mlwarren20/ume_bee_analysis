## *Prunus mume* nectar, and *Apis mellifera* and *Apis cerana* mouth and crop bacteria

This repository includes all the scripts necessary to process and analyze the next-generation sequences and isolate sequences associated with the Japanese apricot nectar and honey bee mouth and crop samples. These scripts will reproduce the results in the corresponding manuscript titled "Bacteria in honeybee crops are decoupled from those in floral nectar and bee mouths", doi:... . In this study, we compared the bacteria found in the nectar, mouth, and crop samples collected summer 2018 and winter 2019 from apricot orchards in Japan. We found that the honey bee crop contains a bacterial community that is more consistent across season and location and distinct from environmental bacteria than previously thought.  

### To run these scripts you will need the following: 
1. bash and python packages: cutadapt v , mafft v , trimAI v , SuperCRUNCH v , and iqtree v2.2.2.6
2. All R packages can be installed by running `renv::restore()` using the renv package.
3. The amplicon sequences available for download from NCBI (PRJNA1076572).
4. The Silva 138.1 prokaryotic SSU taxonomic training data can be downloaded from <https://zenodo.org/records/4587955>.
5. The strain sequences are available to download from NCBI and their accession numbers are available in table S4. 

### To reproduce the processing and analysis: 
1. The raw amplicon sequences are processed first with the bash script "ume_cutadapt.sh" and then cleaned up with the scripts "dada2_processing.Rmd" followed by "Data_cleanup.Rmd". 
2. The phylogeny is reproduced using the script "Ak_phylogeny.sh".
3. The analyses and figures are reproduced in their corresponding scripts: "Figure_1.Rmd", "Figure_2_and_S2.Rmd", "Figures_3_and_S4.Rmd", "Figure_4.Rmd", "Figures_5_S1_and_S3.Rmd", "Figure_6.Rmd", "Figures_7_and_S7.Rmd", "Figures_S5_and_S6.Rmd".

### The extras:
- The script, "helper_functions.R", contains custom functions that were used in the analysis and visualization of the data. 
- There are few steps in the analysis that take a while, so there are saved data objects that can be loaded in order to skip this long step: "seqtab_nochim.rds", "ume_hb_fitGTR.rds", "cmn_dsq_fitGTR.rds".
- The phyloseq object created using iqtree is loaded with "partition.nex.treefile". 
- The "renv" folder and "renv.lock" file load the versions of the packages used in this analysis.

Updated February 29, 2024
