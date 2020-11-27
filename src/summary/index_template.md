---
title: "index"
author: "Almut Lütge"
date: "17 October, 2019"
output: 
  html_document:  
    self_contained: no
    lib_dir: '../../docs/site_libs'
    code_folding: show
    theme: slate
    highlight: tango
  number_sections: no
  toc: yes
  toc_depth: 3
  toc_float:
    collapsed: no
    smooth_scroll: yes
---



A systematic analysis of different kinds of unwanted varaition in single cell RNAseq datasets.

# Analysis

## Batch effects

![](batch_effect_pbmc2_media_files/figure-html/cms celltypes-2.png)

## Simulation

![](simulation_csf_patient_files/figure-html/compare lfcs-4.png)

## CountSimQC

![](external/countsimQC.png)

## Visualize Simulations

![](vis_sim_csf_media_files/figure-html/tsne1__1-1.png)



# Summarized batch characteristics

## Real data
[Summary real data](overall_batch_characteristics.html)
[Variance partition real data](var_part_ct.html)

![](overall_batch_characteristics_files/figure-html/densitycellbench-1.png)


## Simulations
[Summary simulated data](overall_sim_batch_characteristics.html)
[Variance partition simulations](var_part_ct_sim.html)

![](overall_sim_batch_characteristics_files/figure-html/plot batch size-1.png)

# Datasets

Detailed information on reference datasets and preprocessing can be found [here](https://almutlue.github.io/batch_dataset/).

