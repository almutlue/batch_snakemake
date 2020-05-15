# Main Workflow - Characterize batch effects in single cell RNAseq data
#
# Contributors: @almutlue @zjanna

configfile: "config.yaml"

import glob
import os
import logging

sample = glob_wildcards(config["src_data"] + "{sample_name}.rds").sample_name
logger.info("=====")
logger.info(sample)
logger.info("=====")

sim_name = glob_wildcards(config["src_meta_sim"] + "{sim_name}.rds").sim_name
logger.info("=====")
logger.info(sim_name)
logger.info("=====")

sample_sim = sample.copy()
sample_sim.remove("pancreas")
logger.info("=====")
logger.info(sample_sim)
logger.info("=====")

# --- Subworkflows --- #

## -------------------------------------------------------------------------- ##
## Characterize batch effects
## -------------------------------------------------------------------------- ##

subworkflow characterization:
   workdir: config["ROOT"]
   snakefile:  config["src_characterization"] + "snakefile_char"


   
## -------------------------------------------------------------------------- ##
## Simulation
## -------------------------------------------------------------------------- ##

subworkflow simulation:
   workdir: config["ROOT"]
   snakefile:  config["src_sim"] + "snakefile_sim"

## -------------------------------------------------------------------------- ##
## QC Simualtion
## -------------------------------------------------------------------------- ##




# --- Build --- #

## -------------------------------------------------------------------------- ##
## All
## -------------------------------------------------------------------------- ##
rule all:
	input:
	  characterization(expand(config["out_set"] + "settings_{sample}.rds", sample = sample)),
	  characterization(expand(config["out_norm"] + "norm_{sample}_sce.rds", sample = sample)),
	  characterization(expand(config["out_vp"] + "vp_{sample}_sce.rds", sample = sample)),
	  characterization(expand(config["out_abund"] + "abundance_{sample}.rds", sample = sample)),
	  characterization(expand(config["out_cms"] + "cms_{sample}_sce.rds", sample = sample)),
	  characterization(expand(config["out_de"] + "de_{sample}.rds", sample = sample)),
	  characterization(expand(config["out_de"] + "de_{sample}_sce.rds", sample = sample)),
	  characterization(expand(config["out_type"] + "type_{sample}_sce.rds", sample = sample)), 
	  characterization(expand(config["out_summary"] + "summary_{sample}.rds", sample = sample)),
	  characterization(expand(config["docs"] + "batch_effect_{sample}.html", sample = sample)),
	  characterization(config["docs"] + "overall_batch_characteristics.html"),
	  characterization(config["docs"] + "batch_index.html"),
	  characterization(config["docs"] + "dataset_index.html"),
	  simulation(expand(config["out_edger"] + "edgeR_{sample}.rds", sample = sample)),
	  simulation(expand(config["out_sim"] + "{sample}/sim_{sample}_{sim_name}_sce.rds", sample = sample, sim_name = sim_name)),
	  simulation(expand(config["out_sim_char"] + "{sample}/sim_{sample}_{sim_name}_sce.rds", sample = sample, sim_name = sim_name)),
	  simulation(expand(config["docs"] + "countSimQC_{sample_sim}.html", sample_sim = sample_sim)),
	  simulation(expand(config["out_adjust_params"] + "adjust_params_{sample_sim}.rds", sample_sim = sample_sim)),
	  simulation(expand(config["out_sim_vp"] + "vp_sim_{sample_s}_sce.rds", sample_s = sample_sim)),
	  simulation(expand(config["out_sim_abund"] + "abundance_sim_{sample_sim}.rds", sample_sim = sample_sim)),
	  simulation(expand(config["out_sim_de"] + "de_sim_{sample_sim}_sce.rds", sample_sim = sample_sim)),
	  simulation(expand(config["out_sim_de"] + "de_sim_{sample_sim}.rds", sample_sim = sample_sim)),
	  simulation(expand(config["out_sim_type"] + "type_sim_{sample_sim}_sce.rds", sample_sim = sample_sim)),
	  simulation(expand(config["out_sim_sumcms"] + "summarize_cms_sim_{sample_sim}_sce.rds", sample_sim = sample_sim)),
	  simulation(expand(config["out_sim_summary"] + "summary_sim_{sample_sim}.rds", sample_sim = sample_sim)),
	  simulation(expand(config["docs"] + "simulation_{sample_sim}.html", sample_sim = sample_sim)),
	  simulation(expand(config["docs"] + "vis_sim_{sample_sim}.html", sample_sim = sample_sim)),
	  simulation(config["docs"] + "overall_sim_batch_characteristics.html"),
	  simulation(config["docs"] + "simulation_index.html"),
	  simulation(config["docs"] + "sim_char_index.html"),
	  simulation(config["docs"] + "countSimQC_index.html"),
	  simulation(config["docs"] + "overview_index.html")
	  


# --- Optional Rules  --- #

## -------------------------------------------------------------------------- ##
## Directory setup
## -------------------------------------------------------------------------- ##
rule dir_setup:
    input: 
        script = config["src_data_mgt"] + "dir_setup.R"
    params:
        dir_names = config["dir_names"].replace(" ", "")
    log:
        config["log_dir"] + "dir_setup.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args dir_names='{params.dir_names}'" {input.script} {log}'''
        
        
## -------------------------------------------------------------------------- ##
## generate sim_vars
## -------------------------------------------------------------------------- ##
rule generate_sim_vars:
    input: 
        script = config["src_data_mgt"] + "generate_sim_vars.R"
    params:
        out_path = config["src_meta_sim"]
    log:
        config["log_sim_vars"] + "generate_sim_vars.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args out_path='{params.out_path}'" {input.script} {log}'''
        
        
        
