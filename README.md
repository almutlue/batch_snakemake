# Batch snakemake

## General

### Aim

This pipeline is build to explore batch effects in single cell RNAseq data sets and simulate realistic batch effects from them.
It also includes a comparison of the batch effects and other dataset features between simulation and real data.

### Structure

The pipeline consists of 3 major steps 
1. Characterize batch effects in real data
2. Simulate single cell data with a corresponding batch effect
3. Validate simulation using CountsimQC and batch characterization

## Batch effects

### Definition
As batch effect we consider all kinds of unwanted variation. Thus a batch effect is a signal caused by something that is not the biological signal of interest, but conflicts with this signal. So we need to adjust and/or understand the batch effect in order to use the full potential of your biological signal. In this definition a batch effect could be caused by patient differences or media differences in the one case, while in other cases this is the signal of interest. So it is very variable and always depends on the question asked.

### Analysis
We analysed single cell RNAseq dataset with batch effects from different sources
+ Technical batch effects (e.g. different sequencing protocols)
+ Biological batch effects (e.g. different patients)
+ Conditional batch effects (e.g. different media)

## Results
View results [here](https://almutlue.github.io/batch_snakemake/index.html).

## Setup

### Preparations

To setup this pipeline follow these instructions (Step 1 -2 explain one possible way to setup and run snakemake):

1. Set up and activate an Anaconda enviroment with __Snakemake >= v.5.6.0__ (or sth. eqivalent)  
2. Make sure your path to __R__ is exported within snakemake  
  * e.g. adding `*export PATH="/your/prefered/R/bin:$PATH"*` in your `*~/.bashrc*`
3. Install all required R packages using **packrat**
4. Clone this repository 
  * Caution: If you don't want to get all analysis that came with this repo you need to clean the `docs` directory from all files except `_site.yaml`
5. Create `**log**` and `**out**` directories.
6. Run: `*snakemake dir_setup*` to set up the neccessary directory structure to make all rules work.
7. If you want to view or share your analysis as website, activate github pages within your corresponding repo and specify the `*/docs*` as source directory. 

### Run
To run the entire pipeline:
1. Copy your preprocessed `*SingleCellExperiment*` dataset into `*/src/data/*`
2. Generate a corresponding metadata file and save it at `*/src/meta_files/*`
3. Run **snakemake**
4. Push results to github and refresh it's web deployment.
