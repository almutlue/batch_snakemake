# Main Workflow - Characterize batch effects in single cell RNAseq data
#
# Contributors: @almutlue @zjanna

configfile: "config.yaml"

import glob
import os
import logging

sample = glob_wildcards(config["src_data"] + "{sample_name}.rds").sample_name
print("=====")
print(sample)
print("=====")

sim_name = glob_wildcards(config["src_meta_sim"] + "{sim_name}.rds").sim_name
print("=====")
print(sim_name)
print("=====")


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
	  simulation(expand(config["out_edger"] + "edgeR_{sample}.rds", sample = sample)),
	  simulation(expand(config["out_sim"] + "{sample}/sim_{sample}_{sim_name}_sce.rds", sample = sample, sim_name = sim_name)),
	  simulation(expand(config["out_sim_char"] + "{sample}/sim_{sample}_{sim_name}_sce.rds", sample = sample, sim_name = sim_name)),
	  simulation(expand(config["docs"] + "countSimQC_{sample}.html", sample = sample)),
	  simulation(expand(config["out_adjust_params"] + "adjust_params_{sample}.rds", sample = sample)),
	  simulation(expand(config["out_sim_vp"] + "vp_sim_{sample}_sce.rds", sample = sample)),
	  simulation(expand(config["out_sim_abund"] + "abundance_sim_{sample}.rds", sample = sample)),
	  simulation(expand(config["out_sim_de"] + "de_sim_{sample}_sce.rds", sample = sample)),
	  simulation(expand(config["out_sim_de"] + "de_sim_{sample}.rds", sample = sample)),
	  simulation(expand(config["out_sim_type"] + "type_sim_{sample}_sce.rds", sample = sample)),
	  simulation(expand(config["out_sim_sumcms"] + "summarize_cms_sim_{sample}_sce.rds", sample = sample)),
	  simulation(expand(config["out_sim_summary"] + "summary_sim_{sample}.rds", sample = sample)),
	  simulation(expand(config["docs"] + "simulation_{sample}.html", sample = sample)),
	  characterization(config["docs"] + "index.md")
	  


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
        
        
        
