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



# --- Build --- #

## -------------------------------------------------------------------------- ##
## All
## -------------------------------------------------------------------------- ##
rule all:
	input:
	  expand(config["out_set"] + "settings_{sample}.rds", sample = sample),
	  expand(config["out_norm"] + "norm_{sample}_sce.rds", sample = sample),
	  expand(config["out_vp"] + "vp_{sample}_sce.rds", sample = sample),
	  expand(config["out_abund"] + "abundance_{sample}.rds", sample = sample),
	  expand(config["out_cms"] + "cms_{sample}_sce.rds", sample = sample),
	  expand(config["out_de"] + "de_{sample}.rds", sample = sample),
	  expand(config["out_type"] + "type_{sample}_sce.rds", sample = sample), 
	  expand(config["out_summary"] + "summary_{sample}.rds", sample = sample),
	  expand(config["docs"] + "batch_effect_{sample}.html", sample = sample),
	  config["docs"] + "index.md"

# --- Main Build Rules --- #

## -------------------------------------------------------------------------- ##
## Settings
## -------------------------------------------------------------------------- ##
 
rule settings:
    input: 
        data = config["src_data"] + "{sample}.rds",
        meta = config["src_meta"] + "{sample}_meta.rds",
        script = config["src_data_mgt"] + "settings.R"
    output:
        out = config["out_set"] + "settings_{sample}.rds"
    log:
        config["log_set"] + "prepare_settings_{sample}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' meta='{input.meta}' outputfile='{output.out}'" {input.script} {log}'''


## -------------------------------------------------------------------------- ##
## Normalization and input checks
## -------------------------------------------------------------------------- ##

rule normalize:
    input: 
        data = config["src_data"] + "{sample}.rds",
        param = config["out_set"] + "settings_{sample}.rds",
        script = config["src_data_mgt"] + "normalization.R"
    output:
        out = config["out_norm"] + "norm_{sample}_sce.rds"
    log:
        config["log_norm"] + "normalization_{sample}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' params='{input.param}' outputfile='{output.out}'" {input.script} {log}'''
        
        
## -------------------------------------------------------------------------- ##
## Variance partitioning
## -------------------------------------------------------------------------- ##

rule variance_part:
    input: 
        data = config["out_norm"] + "norm_{sample}_sce.rds",
        param = config["out_set"] + "settings_{sample}.rds",
        script = config["src_analysis"] + "variance_part.R"
    output:
        out = config["out_vp"] + "vp_{sample}_sce.rds"
    log:
        config["log_vp"] + "vp_{sample}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' params='{input.param}' outputfile='{output.out}'" {input.script} {log}'''
        
        
## -------------------------------------------------------------------------- ##
## Differential abundance
## -------------------------------------------------------------------------- ##

rule abundance:
    input: 
        data = config["out_vp"] + "vp_{sample}_sce.rds",
        param = config["out_set"] + "settings_{sample}.rds",
        script = config["src_analysis"] + "abundance.R"
    output:
        out = config["out_abund"] + "abundance_{sample}.rds"
    log:
        config["log_abund"] + "abundance_{sample}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' params='{input.param}' outputfile='{output.out}'" {input.script} {log}'''


## -------------------------------------------------------------------------- ##
## Cellspecific Mixing Score
## -------------------------------------------------------------------------- ##

rule cms:
    input: 
        data = config["out_vp"] + "vp_{sample}_sce.rds",
        param = config["out_set"] + "settings_{sample}.rds",
        script = config["src_analysis"] + "cms.R"
    output:
        out = config["out_cms"] + "cms_{sample}_sce.rds"
    log:
        config["log_cms"] + "cms_{sample}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' params='{input.param}' outputfile='{output.out}'" {input.script} {log}'''
        
        
        
## -------------------------------------------------------------------------- ##
## Differential expression
## -------------------------------------------------------------------------- ##

rule de:
    input: 
        data = config["out_cms"] + "cms_{sample}_sce.rds",
        param = config["out_set"] + "settings_{sample}.rds",
        script = config["src_analysis"] + "de.R"
    output:
        out = config["out_de"] + "de_{sample}.rds"
    log:
        config["log_de"] + "de_{sample}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' params='{input.param}' outputfile='{output.out}'" {input.script} {log}'''
        
        

## -------------------------------------------------------------------------- ##
## Batch type
## -------------------------------------------------------------------------- ##

rule batch_type:
    input: 
        data = config["out_cms"] + "cms_{sample}_sce.rds",
        param = config["out_set"] + "settings_{sample}.rds",
        script = config["src_analysis"] + "batch_type.R"
    output:
        out = config["out_type"] + "type_{sample}_sce.rds"
    log:
        config["log_type"] + "type_{sample}.Rout"
    shell:
        '''R CMD BATCH --no-restore --no-save "--args data='{input.data}' params='{input.param}' outputfile='{output.out}'" {input.script} {log}'''
    
    
    
## -------------------------------------------------------------------------- ##
## Summarize batch characterization
## -------------------------------------------------------------------------- ##

rule summary:
    input: 
        data = config["out_type"] + "type_{sample}_sce.rds",
        param = config["out_set"] + "settings_{sample}.rds",
        de = config["out_de"] + "de_{sample}.rds",
        gs = config["src_data"] + "c5.mf.v6.2.symbols.gmt",
        abund = config["out_abund"] + "abundance_{sample}.rds",
        script = config["src_summary"] + "summary_batch_characterization.Rmd",
        outdir = config["docs"]
    output:
        out = config["docs"] + "batch_effect_{sample}.html",
        out_file = config["out_summary"] + "summary_{sample}.rds"
    shell:
        '''R -e "rmarkdown::render(input ='{input.script}', output_file=basename('{output.out}'), output_dir='{input.outdir}', knit_root_dir=getwd(), params = list(data='{input.data}', param ='{input.param}', de ='{input.de}', gs = '{input.gs}', abund = '{input.abund}', out_file = '{output.out_file}'))"'''
        

## -------------------------------------------------------------------------- ##
## Update index file
## -------------------------------------------------------------------------- ##

rule index:
    input:
        script = config["src_summary"] + "write_index.sh",
        ind_template = config["src_summary"] + "index_template.md"
    output:
        out = config["docs"] + "index.md"
    shell:
        "./{input.script} {input.ind_template} {output.out} {sample}"


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
        
        
        
